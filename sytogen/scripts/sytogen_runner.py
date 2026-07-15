"""
sytogen_runner.py
=================
Top-level pipeline entry point called by api.py:

    result = run_sytogen_pipeline(seq_record, codon_df, motif_df, params)

Returns a dict with:
    altered_fasta   : str   — the RM-silent sequence in FASTA format
    original_fasta  : str   — the input sequence in FASTA format
    decision_matrix : list  — one row per candidate considered (see _MatrixRow)
    summary         : dict  — counts of motifs resolved, edits applied, etc.

The decision matrix is the main output the user sees. Every synonymous
candidate considered for every motif position is recorded, with a 'chosen'
flag on the winner, so the user can reconstruct exactly why each edit was made.
"""

import io
import csv
import pandas as pd
import Bio.Seq
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

from sytogen.scripts.genome_model import (
    GenomeModel,
    Gene,
    Motif,
    ProtectedRegion,
    compile_iupac,
    motif_destroyed,
    is_gc_preserving_swap,
)
from sytogen.scripts.assembly_planner import (
    fragment_sequence,
    DEFAULT_FRAGMENT_SIZE,
    DEFAULT_OVERLAP_LENGTH,
)


# ============================================================
# PUBLIC ENTRY POINT
# ============================================================

def _final_new_motif_check(original_sequence, final_sequence, motifs, topology):
    """
    Whole-construct safety net, run once after ALL edits are applied.

    evaluate_mutation()'s "creates new motif" check only ever looks in a
    small (~50bp) window around each individual edit, evaluated one edit
    at a time as candidates are considered. This function is deliberately
    independent of that — a full re-scan of the entire final sequence for
    every distinct target pattern, compared against the same re-scan of
    the original sequence. Anything present in the final scan that wasn't
    in the original is a newly introduced site: whether that's from a
    combined effect of two nearby edits neither window caught alone, or
    just a second, unrelated verification pass on the incremental
    per-edit bookkeeping.

    Returns a list of {'motif','start','end','strand'} dicts for any
    newly introduced site — empty if none.
    """
    unique_patterns = sorted({m.motif for m in motifs})
    if not unique_patterns:
        return []

    pattern_df = pd.DataFrame({"motif": unique_patterns})
    # No 'start'/'end' columns -> _parse_motifs does a full, exhaustive,
    # topology-aware self-search rather than trusting known positions.
    original_hits = _parse_motifs(pattern_df, original_sequence, topology)
    final_hits    = _parse_motifs(pattern_df, final_sequence, topology)

    def _hit_key(m):
        return (m.motif, m.start, m.end, m.strand)

    original_keys = {_hit_key(m) for m in original_hits}
    new_hits = [m for m in final_hits if _hit_key(m) not in original_keys]

    return [
        {"motif": m.motif, "start": m.start, "end": m.end, "strand": m.strand}
        for m in new_hits
    ]


def run_sytogen_pipeline(seq_record, codon_df, motif_df, params=None):
    """
    Parameters
    ----------
    seq_record : Bio.SeqRecord
        MotifFinder-annotated GenBank record.
    codon_df : pd.DataFrame
        CodonBias output. Expected columns: 'codon', 'fraction' (or 'frequency').
    motif_df : pd.DataFrame
        Restriction motif table. Expected columns: 'motif', 'start', 'end',
        optionally 'strand'.
    params : dict, optional
        'topology'             : 'circular' | 'linear'  (default 'circular')
        'preserve_gc'          : bool                    (default False, reserved)
        'include_assembly_plan': bool                    (default False) — if
            True, also runs Gibson Assembly fragment/overlap planning
            (assembly_planner.fragment_sequence) on the final RM-silent
            sequence and returns it under the 'assembly_plan' key.
        'fragment_size'        : int  (default assembly_planner.DEFAULT_FRAGMENT_SIZE)
        'overlap_length'       : int  (default assembly_planner.DEFAULT_OVERLAP_LENGTH)

    Returns
    -------
    dict with keys: altered_fasta, original_fasta, decision_matrix, summary,
    assembly_plan (None unless 'include_assembly_plan' was set)
    """
    params = params or {}
    topology = params.get("topology", "circular")

    sequence = str(seq_record.seq).upper()
    seqid    = seq_record.id or "sequence"

    # ----------------------------------------------------------
    # 1. Validate inputs, then parse
    # ----------------------------------------------------------
    # Without this, a motif table missing its 'motif' column (wrong file
    # uploaded, wrong delimiter, a header row that didn't parse the way
    # the person expected, etc.) makes _parse_motifs silently skip every
    # row and return an empty list. The pipeline then "succeeds" with
    # motifs_input: 0 and an unmodified sequence — indistinguishable from
    # a construct that genuinely has zero restriction sites, which is a
    # very different situation to hand back silently. Same story for the
    # codon table and its 'codon' column: _parse_codon_usage would
    # silently return {}, and every synonymous choice downstream would
    # become an arbitrary tie instead of an actual codon-usage decision.
    _validate_motif_table(motif_df)
    _validate_codon_table(codon_df)

    genes             = _parse_genes(seq_record)
    motifs            = _parse_motifs(motif_df, sequence, topology)
    protected_regions = _parse_protected_regions(seq_record)
    codon_usage       = _parse_codon_usage(codon_df)

    # ----------------------------------------------------------
    # 2. Build genome model
    # ----------------------------------------------------------
    genome = GenomeModel(
        sequence=sequence,
        topology=topology,
        genes=genes,
        motifs=motifs,
        protected_regions=protected_regions,
        codon_usage=codon_usage,
    )

    # ----------------------------------------------------------
    # 3. Run motif-by-motif optimization, collecting full matrix
    # ----------------------------------------------------------
    decision_matrix   = []   # every candidate row (including losers)
    applied_mutations = []   # only the chosen edits
    motifs_resolved   = 0
    motifs_unresolved = 0
    # Explicit resolution status per motif, keyed the same way decision
    # matrix rows are grouped (motif, start, end, strand). Tracked
    # directly here rather than re-derived from decision_matrix rows,
    # because a motif resolved as a side effect of a DIFFERENT edit (the
    # motif_destroyed() check below, hit before any candidate is even
    # generated for it) never gets a row of its own — there'd be nothing
    # to re-derive its status from otherwise.
    def _motif_key(m):
        return (m.motif, m.start, m.end, m.strand)
    resolved_motif_keys = set()

    for motif in motifs:
        if motif_destroyed(genome, motif):
            # Already gone — a prior edit may have cleared it
            motifs_resolved += 1
            resolved_motif_keys.add(_motif_key(motif))
            continue

        # Generate every valid candidate for this motif: synonymous
        # codon substitutions inside genes, plus destructive-but-safe
        # single-base substitutions at unprotected non-coding positions.
        candidates = (
            genome.generate_synonymous_candidates(motif)
            + genome.generate_neutral_candidates(motif)
        )

        if not candidates:
            motifs_unresolved += 1
            diagnostic = genome.explain_no_candidates(motif)
            _record_unresolvable(decision_matrix, motif, diagnostic)
            continue

        # Score all candidates. GC-class preservation (is_gc_preserving_swap)
        # is a secondary key here, deliberately — it only decides between
        # candidates whose primary score (destroyed/usage/edits) is
        # exactly tied. A real usage_score difference always wins; GC
        # preservation never overrides it, only breaks a tie between
        # equally-good options.
        scored = [
            (c, genome.score_candidate(c), is_gc_preserving_swap(c.mutation.old, c.mutation.new))
            for c in candidates
        ]
        scored.sort(key=lambda x: (x[1], x[2]), reverse=True)
        best_candidate, best_score, _ = scored[0]

        # Record every candidate in the matrix, mark the winner
        for candidate, score, gc_preserving in scored:
            chosen = (candidate is best_candidate)
            decision_matrix.append(
                _make_matrix_row(motif, candidate, score, chosen, genome,
                                  best_candidate=best_candidate, best_score=best_score,
                                  gc_preserving=gc_preserving)
            )

        # Apply the winning edit to the live genome
        try:
            genome.sequence = genome.apply_mutation(best_candidate.mutation)
            genome.topology_engine = genome.build_topology(genome.sequence)
            applied_mutations.append(best_candidate.mutation)
        except ValueError as e:
            # Sequence has drifted from what the candidate expected —
            # a prior overlapping edit changed this region. Skip and note it.
            motifs_unresolved += 1
            _mark_last_rows_as_skipped(decision_matrix, motif)
            continue

        if motif_destroyed(genome, motif):
            motifs_resolved += 1
            resolved_motif_keys.add(_motif_key(motif))
        else:
            # Edit applied but motif persists — needs another pass
            motifs_unresolved += 1

    # ----------------------------------------------------------
    # 3b. Final whole-construct check: any new motifs introduced anywhere?
    # ----------------------------------------------------------
    # Independent of evaluate_mutation()'s per-edit window check — see
    # _final_new_motif_check's docstring for why that's not redundant.
    new_motifs = _final_new_motif_check(sequence, genome.sequence, motifs, topology)

    # ----------------------------------------------------------
    # 4. Optional: Gibson Assembly fragment/overlap planning
    # ----------------------------------------------------------
    # Runs against the *final* edited genome/sequence, using the decision
    # matrix so fragment boundaries stay clear of the positions SyToGen
    # just edited (collect_edit_positions() reads the same 'chosen' /
    # 'edit_position' columns _make_matrix_row() already produces).
    assembly_plan = None
    if params.get("include_assembly_plan"):
        assembly_plan = fragment_sequence(
            genome.sequence,
            genome,
            decision_matrix,
            fragment_size=params.get("fragment_size", DEFAULT_FRAGMENT_SIZE),
            overlap_length=params.get("overlap_length", DEFAULT_OVERLAP_LENGTH),
            topology=topology,
        )

    # ----------------------------------------------------------
    # 5. Serialise outputs
    # ----------------------------------------------------------
    altered_fasta  = _to_fasta(genome.sequence, seqid + "_sytogen")
    original_fasta = _to_fasta(sequence, seqid + "_original")

    summary = {
        "sequence_id":       seqid,
        "topology":          topology,
        "original_length":   len(sequence),
        "altered_length":    len(genome.sequence),
        "motifs_input":      len(motifs),
        "motifs_resolved":   motifs_resolved,
        "motifs_unresolved": motifs_unresolved,
        "edits_applied":     len(applied_mutations),
        "candidates_total":  len(decision_matrix),
        "new_motifs_introduced": len(new_motifs),
    }

    return {
        "altered_fasta":   altered_fasta,
        "original_fasta":  original_fasta,
        "altered_sequence": genome.sequence,
        "applied_mutations": applied_mutations,
        "motifs":          motifs,
        "resolved_motif_keys": resolved_motif_keys,
        "decision_matrix": decision_matrix,
        "summary":         summary,
        "assembly_plan":   assembly_plan,
        "new_motifs":      new_motifs,
    }


# ============================================================
# SERIALISATION HELPERS
# ============================================================

def decision_matrix_to_tsv(matrix):
    """Serialise the decision matrix to a TSV string."""
    if not matrix:
        return ""
    fieldnames = list(matrix[0].keys())
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(matrix)
    return buf.getvalue()


def decision_matrix_to_json(matrix):
    """Serialise the decision matrix to a JSON-serialisable list of dicts."""
    # Already a list of plain dicts — nothing to transform.
    return matrix


def assembly_plan_to_tsv(plan):
    """Serialise an AssemblyPlan's fragments + overlaps + primers to a TSV string."""
    if not plan or not plan.fragments:
        return ""

    fieldnames = [
        "fragment", "start", "end", "length",
        "overlap_left_start", "overlap_left_end", "overlap_left_seq",
        "overlap_left_tm", "overlap_left_gc", "overlap_left_score",
        "overlap_right_start", "overlap_right_end", "overlap_right_seq",
        "overlap_right_tm", "overlap_right_gc", "overlap_right_score",
        "forward_primer_seq", "forward_primer_tm",
        "reverse_primer_seq", "reverse_primer_tm",
    ]
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for frag in plan.fragments:
        row = {
            "fragment": frag.name,
            "start":    frag.start,
            "end":      frag.end,
            "length":   len(frag.sequence),
        }
        for prefix, overlap in (
            ("overlap_left",  frag.overlap_left),
            ("overlap_right", frag.overlap_right),
        ):
            if overlap:
                row[f"{prefix}_start"] = overlap.start
                row[f"{prefix}_end"]   = overlap.end
                row[f"{prefix}_seq"]   = overlap.sequence
                row[f"{prefix}_tm"]    = round(overlap.tm, 2)
                row[f"{prefix}_gc"]    = round(overlap.gc_percent, 2)
                row[f"{prefix}_score"] = round(overlap.score, 2)
            else:
                row[f"{prefix}_start"] = ""
                row[f"{prefix}_end"]   = ""
                row[f"{prefix}_seq"]   = ""
                row[f"{prefix}_tm"]    = ""
                row[f"{prefix}_gc"]    = ""
                row[f"{prefix}_score"] = ""

        row["forward_primer_seq"] = frag.forward_primer.sequence if frag.forward_primer else ""
        row["forward_primer_tm"]  = round(frag.forward_primer.tm, 2) if frag.forward_primer else ""
        row["reverse_primer_seq"] = frag.reverse_primer.sequence if frag.reverse_primer else ""
        row["reverse_primer_tm"]  = round(frag.reverse_primer.tm, 2) if frag.reverse_primer else ""

        writer.writerow(row)

    return buf.getvalue()


def assembly_primers_to_tsv(plan):
    """
    A flat, order-sheet-ready TSV: one row per primer (two per fragment),
    with the full primer sequence (Gibson homology tail + annealing
    region already concatenated) and its annealing-region Tm.
    """
    if not plan or not plan.fragments:
        return ""

    fieldnames = ["primer_name", "fragment", "role", "sequence", "length", "anneal_tm"]
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for frag in plan.fragments:
        for role, primer in (("forward", frag.forward_primer), ("reverse", frag.reverse_primer)):
            if not primer:
                continue
            writer.writerow({
                "primer_name": primer.name,
                "fragment":    frag.name,
                "role":        role,
                "sequence":    primer.sequence,
                "length":      len(primer.sequence),
                "anneal_tm":   round(primer.tm, 2),
            })

    return buf.getvalue()


def assembly_plan_fragments_fasta(plan):
    """Serialise each fragment's sequence to a multi-record FASTA string."""
    if not plan or not plan.fragments:
        return ""
    return "".join(_to_fasta(frag.sequence, frag.name) for frag in plan.fragments)


def assembly_plan_summary(plan):
    """Small JSON-serialisable summary of the assembly plan as a whole."""
    if not plan:
        return {}
    return {
        "fragment_count":       len(plan.fragments),
        "assembly_score":       round(plan.assembly_score, 4),
        "overlap_length":       plan.overlap_length,
        "target_fragment_size": plan.target_fragment_size,
        "warnings":             plan.warnings,
    }


# ============================================================
# PARSING HELPERS
# ============================================================

GENE_FEATURE_TYPES = {"CDS", "ORF", "Marker"}


def _parse_genes(seq_record):
    """
    Extract Gene objects from protein-coding features in the SeqRecord.

    Not every annotation source uses the strict GenBank 'CDS' type — plasmid
    editors like SnapGene/ApE frequently label coding regions 'ORF', and
    selection-marker genes sometimes show up as 'Marker' (e.g. an AmpR/bla
    beta-lactamase gene). All three are treated as genes here so codon-level
    synonymous editing works anywhere a real protein-coding region exists,
    not just where it happens to be typed exactly 'CDS'.

    A gene that spans the circular origin is written in GenBank as a
    CompoundLocation, e.g. join(9950..10000,1..50). BioPython's
    CompoundLocation.start/.end return the min/max across all parts — for
    an origin-spanning join that's practically the whole molecule (0 to
    10000 here), not the small ~100bp region the gene actually occupies.
    Left uncorrected, that falsely blocks almost every position in the
    plasmid from being edited. Detected below by comparing that naive
    span against the gene's actual declared length (sum of its parts);
    a real gap between them is the signature of a wrap, and the raw,
    un-clamped start/end convention already used for Motif/Overlap
    wraparound is reused here (see Gene._resolve / Gene.codon_start).
    """
    sequence_length = len(seq_record.seq)
    genes = []
    for i, feature in enumerate(seq_record.features):
        if feature.type not in GENE_FEATURE_TYPES:
            continue

        parts = getattr(feature.location, "parts", [feature.location])
        naive_start = int(feature.location.start)
        naive_end   = int(feature.location.end) - 1
        declared_length = sum(len(part) for part in parts)

        if len(parts) > 1 and (naive_end - naive_start + 1) > declared_length:
            # Origin-spanning join(): trust the parts' declared order
            # (GenBank convention lists them in biological reading order —
            # the segment right before the origin, then the segment right
            # after it) rather than re-deriving it from position.
            start = int(parts[0].start)
            end   = sequence_length + int(parts[-1].end) - 1  # raw, >= sequence_length
        else:
            start = naive_start
            end   = naive_end

        strand = "+" if feature.location.strand >= 0 else "-"
        gene_id = (
            feature.qualifiers.get("gene",     [None])[0]
            or feature.qualifiers.get("locus_tag", [None])[0]
            or feature.qualifiers.get("sequence",  [None])[0]  # ORF-style label used by some tools
            or f"{feature.type}_{i}"
        )
        genes.append(Gene(gene_id=gene_id, start=start, end=end, strand=strand))
    return genes


def _validate_motif_table(motif_df):
    """
    Raise a clear, actionable error for a malformed motif table.

    An empty table (0 rows, but the right column present) is a legitimate
    input — it just means "nothing to silence" — and is NOT flagged here.
    What IS flagged is the table missing the 'motif' column entirely,
    which otherwise causes _parse_motifs to silently skip every row.
    """
    if motif_df is None or not hasattr(motif_df, "columns"):
        raise ValueError("Motif table could not be read as a table at all.")

    columns = {str(c).strip().lower() for c in motif_df.columns}
    if "motif" not in columns:
        raise ValueError(
            "Motif table is missing a required 'motif' column (the "
            "recognition sequence). Optional columns: 'start'/'end' or "
            "'position_1based'/'hit_seq', and 'strand'. "
            f"Found columns: {list(motif_df.columns)!r}."
        )


def strip_backbone(seq_record, backbone_record, topology="circular"):
    """
    Remove a known vector backbone from a full construct, returning just
    the insert. Mirrors legacy_sytogen.sequence_preprocess(), adapted to
    this codebase's validate-via-ValueError convention and made
    topology-aware — legacy always used a doubled-sequence search
    regardless of topology, which only makes sense for a circular
    molecule (a linear insert has no origin for a backbone to wrap
    around, so a plain substring search is correct and sufficient there).

    Both the construct and the backbone must contain only A/C/G/T, and
    the backbone must be found in the construct EXACTLY once — zero
    matches means it isn't actually there, more than one is ambiguous
    (which copy is the real backbone?). Same validation legacy applied.

    Returns a new SeqRecord for the insert only. Feature coordinates are
    adjusted (or the feature dropped, if it fell entirely inside the
    removed backbone region) by BioPython's own SeqRecord slicing —
    nothing custom needed there.
    """
    sequence = str(seq_record.seq).upper()
    backbone = str(backbone_record.seq).upper()

    if not backbone:
        raise ValueError("Backbone sequence is empty.")

    non_canonical_seq = set(sequence) - set("ACGT")
    if non_canonical_seq:
        raise ValueError(
            f"Input sequence contains non-canonical bases: {sorted(non_canonical_seq)}."
        )
    non_canonical_backbone = set(backbone) - set("ACGT")
    if non_canonical_backbone:
        raise ValueError(
            f"Backbone sequence contains non-canonical bases: {sorted(non_canonical_backbone)}."
        )

    length = len(sequence)

    if topology == "circular":
        # Double the sequence (enough to catch a backbone that spans the
        # origin as one contiguous match) — same technique _parse_motifs
        # uses for circular motif search.
        search_space = sequence + sequence[:len(backbone) - 1]
        match_count = search_space.count(backbone)
        if match_count == 0:
            raise ValueError("Backbone not found in the input construct.")
        if match_count > 1:
            raise ValueError("Backbone found more than once in the input construct.")

        match_start = search_space.find(backbone)
        start = match_start % length
        end = (match_start + len(backbone)) % length

        if start < end:
            insert_record = seq_record[:start] + seq_record[end:]
        else:
            # Backbone wraps the origin — the insert is the single
            # contiguous piece between where the backbone's tail ends
            # and where its head begins.
            insert_record = seq_record[end:start]
    else:
        # Linear: no wraparound possible, plain substring search.
        match_count = sequence.count(backbone)
        if match_count == 0:
            raise ValueError("Backbone not found in the input construct.")
        if match_count > 1:
            raise ValueError("Backbone found more than once in the input construct.")

        start = sequence.find(backbone)
        end = start + len(backbone)
        insert_record = seq_record[:start] + seq_record[end:]

    if len(insert_record.seq) == 0:
        raise ValueError("Removing the backbone leaves an empty insert — nothing left to process.")

    insert_record.annotations["molecule_type"] = "DNA"
    insert_record.id = seq_record.id
    insert_record.name = seq_record.name
    return insert_record


def _validate_codon_table(codon_df):
    """
    Same reasoning as _validate_motif_table: without this, a codon table
    missing its 'codon' column makes _parse_codon_usage silently return
    an empty usage dict, and every synonymous codon choice downstream
    becomes an arbitrary tie (usage_score=0 for everything) instead of an
    actual strain-specific preference — with nothing telling the person
    their codon-usage file didn't parse the way they expected.
    """
    if codon_df is None or not hasattr(codon_df, "columns"):
        raise ValueError("Codon usage table could not be read as a table at all.")

    columns = {str(c).strip().lower() for c in codon_df.columns}
    if "codon" not in columns:
        raise ValueError(
            "Codon usage table is missing a required 'codon' column. "
            "Also expected one usage-score column: 'fraction', 'frequency', "
            "'value', 'usage', 'proportion', 'ranking_ratio', 'ranking', "
            f"or 'count'. Found columns: {list(codon_df.columns)!r}."
        )

    score_columns = {
        "fraction", "frequency", "value", "usage",
        "proportion", "ranking_ratio", "ranking", "count",
    }
    if not (columns & score_columns):
        raise ValueError(
            "Codon usage table has a 'codon' column but none of the "
            "recognized usage-score columns (expected one of: "
            f"{sorted(score_columns)}). Found columns: {list(codon_df.columns)!r}."
        )


def _parse_motifs(motif_df, sequence, topology="circular"):
    """
    Build Motif objects from the motif table.

    The table must have a 'motif' column with the recognition sequence.
    'start' / 'end' columns are optional — if absent we search the sequence
    ourselves so every occurrence is covered.

    For circular topology, the self-search also has to catch motifs that
    span the origin (e.g. the last 4 bases of the sequence followed by the
    first 6). Those are represented with `end >= len(sequence)` — an
    intentionally un-clamped, "raw" coordinate — rather than dropped or
    wrapped to `end < start`. GenomeModel's topology engine already knows
    how to resolve a raw position like that back onto the circular
    sequence (see `normalize_position` / `get_interval`), so every
    downstream consumer (candidate generation, position indexing,
    motif_destroyed) handles it correctly without further special-casing.
    """
    motif_df = motif_df.rename(
        columns={column: str(column).strip().lower() for column in motif_df.columns}
    )

    motifs = []
    seen   = set()   # deduplicate (motif, start) pairs — shared across
                      # strands so a palindromic hit found on both the
                      # forward and reverse search isn't recorded twice.

    has_coords = (
        "start" in motif_df.columns
        and "end"   in motif_df.columns
    )
    has_motiffinder_coords = "position_1based" in motif_df.columns
    is_circular = topology == "circular"

    for _, row in motif_df.iterrows():
        if "motif" not in motif_df.columns:
            continue

        motif_seq = str(row["motif"]).strip().upper()
        if not motif_seq:
            continue

        strand = str(row.get("strand", "+")).strip() if "strand" in motif_df.columns else "+"

        if has_motiffinder_coords and not _is_empty(row.get("position_1based")):
            start = int(row["position_1based"]) - 1
            hit_seq = str(row.get("hit_seq", "")).strip().upper()
            motif_len = len(hit_seq) if hit_seq else len(motif_seq)
            end = start + motif_len - 1
            key = (motif_seq, start)
            if key not in seen:
                seen.add(key)
                motifs.append(Motif(motif=motif_seq, start=start, end=end, strand=strand))
        elif has_coords and not _is_empty(row.get("start")) and not _is_empty(row.get("end")):
            start = int(row["start"])
            end   = int(row["end"])
            key   = (motif_seq, start)
            if key not in seen:
                seen.add(key)
                motifs.append(Motif(motif=motif_seq, start=start, end=end, strand=strand))
        else:
            # No coordinates — search the sequence for all occurrences.
            regex = compile_iupac(motif_seq)
            L = len(sequence)

            if is_circular:
                # Search a doubled sequence so matches spanning the origin
                # are found; keep only matches that *start* in the first
                # copy (start < L) so each circular occurrence is counted
                # once. The kept end coordinate is left un-clamped (may be
                # >= L) to signal a wrap.
                doubled = sequence + sequence
                for m in regex.finditer(doubled):
                    if m.start() >= L:
                        continue
                    key = (motif_seq, m.start())
                    if key not in seen:
                        seen.add(key)
                        motifs.append(
                            Motif(
                                motif=motif_seq,
                                start=m.start(),
                                end=m.end() - 1,
                                strand="+",
                            )
                        )

                # Reverse complement: search the doubled RC sequence, then
                # map each hit back to a forward-strand genomic start/end.
                # reverse_complement(sequence + sequence) == rc_seq + rc_seq
                # (since both halves are identical), so we can build that
                # directly rather than reverse-complementing the doubled
                # forward sequence.
                rc_seq = _reverse_complement(sequence)
                doubled_rc = rc_seq + rc_seq
                for m in regex.finditer(doubled_rc):
                    match_len = m.end() - m.start()
                    fwd_start = 2 * L - m.start() - match_len
                    if not (0 <= fwd_start < L):
                        continue  # duplicate representative from the second copy
                    fwd_end = fwd_start + match_len - 1
                    key = (motif_seq, fwd_start)
                    if key not in seen:
                        seen.add(key)
                        motifs.append(
                            Motif(
                                motif=motif_seq,
                                start=fwd_start,
                                end=fwd_end,
                                strand="-",
                            )
                        )
            else:
                for m in regex.finditer(sequence):
                    key = (motif_seq, m.start())
                    if key not in seen:
                        seen.add(key)
                        motifs.append(
                            Motif(
                                motif=motif_seq,
                                start=m.start(),
                                end=m.end() - 1,
                                strand="+",
                            )
                        )
                # Also search reverse complement
                rc_seq = _reverse_complement(sequence)
                for m in regex.finditer(rc_seq):
                    # Convert RC coordinate back to forward-strand genomic coordinate
                    fwd_end   = L - m.start() - 1
                    fwd_start = L - m.end()
                    key = (motif_seq, fwd_start)
                    if key not in seen:
                        seen.add(key)
                        motifs.append(
                            Motif(
                                motif=motif_seq,
                                start=fwd_start,
                                end=fwd_end,
                                strand="-",
                            )
                        )

    return motifs


def _is_motiffinder_hit_marker(feature):
    """
    MotifFinder writes each motif it finds back into the GenBank as its own
    misc_feature entry (see api.py /motiffinder/run: qualifiers include
    'ID'='motif_hit_NNNN', 'motif', and 'hit_seq'). These mark exactly the
    sites SyToGen is meant to edit and must never be treated as protected
    regulatory regions — otherwise every motif site would be locked from
    editing by the very annotation that identified it.
    """
    qualifiers = feature.qualifiers
    feature_id = str(qualifiers.get("ID", [""])[0])
    return (
        "motif" in qualifiers
        or "hit_seq" in qualifiers
        or feature_id.startswith("motif_hit_")
    )


def _parse_protected_regions(seq_record, max_protected_length=100):
    """
    Extract ProtectedRegion objects from regulatory / misc_feature annotations.

    Only features no longer than max_protected_length (bp) are treated as
    protected. Narrow regulatory elements — a promoter box, an RBS, an
    operator sequence, a single protein-binding/recognition site — are
    exactly the kind of thing that genuinely needs to stay untouched at the
    nucleotide level, and are typically well under 100 bp.

    Broad, whole-region descriptive annotations (e.g. "Staphylococcus aureus
    plasmid pC194 region" spanning 2.6 kb, or an entire inducible-promoter
    expression cassette spanning 1.4 kb) are informational labels about
    provenance/function, not narrow no-edit zones — and since they can span
    entire genes, treating them as absolute editing blocks can silently
    veto every synonymous candidate in a gene without any biological reason
    to. Anything over the cutoff is excluded from protection so it doesn't
    block synonymous codon edits it was never really meant to guard against.
    """
    protected = []
    PROTECTED_TYPES = {"regulatory", "misc_feature", "rep_origin", "promoter", "RBS"}
    for feature in seq_record.features:
        if feature.type not in PROTECTED_TYPES:
            continue
        if _is_motiffinder_hit_marker(feature):
            continue
        start = int(feature.location.start)
        end   = int(feature.location.end) - 1
        if (end - start + 1) > max_protected_length:
            continue
        protected.append(ProtectedRegion(start=start, end=end))
    return protected


def _parse_codon_usage(codon_df):
    """
    Convert a CodonBias DataFrame into the {codon: float} dict GenomeModel expects.

    Accepts column names 'fraction', 'frequency', or 'value' for the usage score.

    Every consumer of this dict (Gene.ranked_synonymous_codons,
    GenomeModel.score_candidate) assumes a HIGHER usage_score means a MORE
    preferred codon. That's true for 'fraction'/'frequency'/'value'/'usage'/
    'proportion'/'count', but CodonBias's own output (see
    legacy_sytogen.codon_usage()) also emits 'Ranking' and 'Ranking_ratio',
    where LOWER is better (rank 1 = most-used codon). Those two must be
    inverted here, or SyToGen would silently prefer the rarest codon instead
    of the most strain-preferred one whenever a CodonBias export only
    includes Ranking-based columns.
    """
    codon_df = codon_df.rename(
        columns={column: str(column).strip().lower() for column in codon_df.columns}
    )

    usage = {}
    score_col = None
    for candidate_col in (
        "fraction",
        "frequency",
        "value",
        "usage",
        "proportion",
        "ranking_ratio",
        "ranking",
        "count",
    ):
        if candidate_col in codon_df.columns:
            score_col = candidate_col
            break

    if "codon" not in codon_df.columns or score_col is None:
        # Graceful fallback — equal usage for all codons
        return {}

    # 'ranking' / 'ranking_ratio' are lower-is-better; every other supported
    # column is higher-is-better. Negate the rank-based columns so a bigger
    # usage_score always means a more preferred codon, consistently.
    invert = score_col in ("ranking", "ranking_ratio")

    for _, row in codon_df.iterrows():
        codon = str(row["codon"]).strip().upper()
        try:
            value = float(row[score_col])
        except (ValueError, TypeError):
            continue
        usage[codon] = -value if invert else value

    return usage


# ============================================================
# DECISION MATRIX ROW BUILDERS
# ============================================================

def _make_matrix_row(motif, candidate, score, chosen, genome, best_candidate=None, best_score=None, gc_preserving=None):
    """Build one row of the decision matrix for a candidate edit."""
    gene = genome.find_gene(candidate.mutation.position)
    is_coding = gene is not None
    aa = _translate(candidate.codon) if is_coding else ""
    destroyed = candidate.result.get("destroyed", 0)
    created   = candidate.result.get("created",   0)
    if gc_preserving is None:
        gc_preserving = is_gc_preserving_swap(candidate.mutation.old, candidate.mutation.new)

    if is_coding:
        change_desc = f"{candidate.codon}\u2192{candidate.replacement}"
        context_desc = f"codon-usage score {candidate.usage_score:.3f}"
    else:
        change_desc = f"{candidate.codon}\u2192{candidate.replacement} (non-coding position)"
        context_desc = "no codon-usage concept outside a gene"

    # destroyed/created counts already live in their own columns below, so
    # the prose here doesn't restate them — avoids the two ever silently
    # disagreeing and keeps this column to context that isn't already
    # structured (which codon/base, why it beat the alternatives).
    if chosen:
        reasoning = (
            f"Chosen: {change_desc}, {context_desc} "
            f"(highest-scoring valid option for this motif)."
        )
    else:
        reasoning = (
            f"Valid alternative ({change_desc}), but scored lower than the chosen candidate."
        )

    return {
        # Motif context
        "motif":             motif.motif,
        "motif_start":       motif.start,
        "motif_end":         motif.end,
        "motif_strand":      motif.strand,
        # Genomic position of the edit
        "edit_position":     candidate.mutation.position,
        "gene_id":           gene.id if gene else "",
        "gene_strand":       gene.strand if gene else "",
        # Codon change (or single-base change, for non-coding positions)
        "original_codon":    candidate.codon,
        "replacement_codon": candidate.replacement,
        "AA_LetterCode":     aa,
        "synonymous":        True if is_coding else "",   # not applicable outside a gene
        # Scoring breakdown — the columns the user cares about
        "motifs_destroyed":  destroyed,
        "reasoning":         reasoning,
        "motifs_created":    created,
        "usage_score":       round(candidate.usage_score, 6),
        "gc_preserving":     gc_preserving,
        "total_score":       round(score, 4),
        # Decision
        "chosen":            chosen,
        "skip_reason":       "",
        # Diagnostic columns — only populated on the sentinel "no valid
        # candidate" rows built by _record_unresolvable(); blank here so
        # every row in the matrix shares the same columns.
        "attempted_count":       "",
        "rejected_count":        "",
        "top_rejection_reason":  "",
        "top_rejection_count":   "",
    }


def _record_unresolvable(matrix, motif, diagnostic):
    """
    Record a sentinel row for a motif where no valid candidate existed.
    `diagnostic` is the dict returned by GenomeModel.explain_no_candidates():
    reason_code, reasoning, and (when applicable) attempted_count/
    rejected_count/top_rejection_reason/top_rejection_count.
    """
    matrix.append({
        "motif":             motif.motif,
        "motif_start":       motif.start,
        "motif_end":         motif.end,
        "motif_strand":      motif.strand,
        "edit_position":     "",
        "gene_id":           "",
        "gene_strand":       "",
        "original_codon":    "",
        "replacement_codon": "",
        "AA_LetterCode":     "",
        "synonymous":        "",
        "motifs_destroyed":  0,
        "reasoning":         diagnostic.get("reasoning") or "No valid candidate could be constructed for this motif.",
        "motifs_created":    0,
        "usage_score":       0,
        "gc_preserving":     "",
        "total_score":       0,
        "chosen":            False,
        "skip_reason":       diagnostic.get("reason_code", "no_valid_candidate"),
        "attempted_count":       diagnostic.get("attempted_count") if diagnostic.get("attempted_count") is not None else "",
        "rejected_count":        diagnostic.get("rejected_count") if diagnostic.get("rejected_count") is not None else "",
        "top_rejection_reason":  diagnostic.get("top_rejection_reason") or "",
        "top_rejection_count":   diagnostic.get("top_rejection_count") if diagnostic.get("top_rejection_count") is not None else "",
    })


def _mark_last_rows_as_skipped(matrix, motif):
    """
    After a sequence-drift error, mark the rows we just appended for this
    motif as not applied so the user knows the edit didn't go through.
    """
    for row in reversed(matrix):
        if row["motif"] == motif.motif and row["motif_start"] == motif.start:
            row["chosen"]      = False
            row["skip_reason"] = "sequence_drift"
            row["reasoning"]   = (
                "This candidate could not be applied: an earlier motif's edit "
                "changed the sequence at or near this position, so the base "
                "this edit expected to find is no longer there."
            )
        else:
            break


# ============================================================
# UTILITY
# ============================================================

_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")

def _reverse_complement(seq):
    return seq.translate(_RC_TABLE)[::-1]


def _translate(codon):
    try:
        return str(Bio.Seq.Seq(codon).translate())
    except Exception:
        return "?"


def _to_fasta(sequence, seqid):
    lines = [f">{seqid}"]
    for i in range(0, len(sequence), 60):
        lines.append(sequence[i:i + 60])
    return "\n".join(lines) + "\n"


def _is_empty(val):
    if val is None:
        return True
    try:
        import math
        return math.isnan(float(val))
    except (TypeError, ValueError):
        return str(val).strip() == ""
