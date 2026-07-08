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
)


# ============================================================
# PUBLIC ENTRY POINT
# ============================================================

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
        'topology'    : 'circular' | 'linear'  (default 'circular')
        'preserve_gc' : bool                    (default False, reserved)

    Returns
    -------
    dict with keys: altered_fasta, original_fasta, decision_matrix, summary
    """
    params = params or {}
    topology = params.get("topology", "circular")

    sequence = str(seq_record.seq).upper()
    seqid    = seq_record.id or "sequence"

    # ----------------------------------------------------------
    # 1. Parse inputs
    # ----------------------------------------------------------
    genes             = _parse_genes(seq_record)
    motifs            = _parse_motifs(motif_df, sequence)
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

    for motif in motifs:
        if motif_destroyed(genome, motif):
            # Already gone — a prior edit may have cleared it
            motifs_resolved += 1
            continue

        # Generate every valid synonymous candidate for this motif
        candidates = genome.generate_synonymous_candidates(motif)

        if not candidates:
            motifs_unresolved += 1
            _record_unresolvable(decision_matrix, motif)
            continue

        # Score all candidates
        scored = [
            (c, genome.score_candidate(c))
            for c in candidates
        ]
        scored.sort(key=lambda x: x[1], reverse=True)
        best_candidate, best_score = scored[0]

        # Record every candidate in the matrix, mark the winner
        for candidate, score in scored:
            chosen = (candidate is best_candidate)
            decision_matrix.append(
                _make_matrix_row(motif, candidate, score, chosen, genome)
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
        else:
            # Edit applied but motif persists — needs another pass
            motifs_unresolved += 1

    # ----------------------------------------------------------
    # 4. Serialise outputs
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
    }

    return {
        "altered_fasta":   altered_fasta,
        "original_fasta":  original_fasta,
        "altered_sequence": genome.sequence,
        "applied_mutations": applied_mutations,
        "decision_matrix": decision_matrix,
        "summary":         summary,
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
    """
    genes = []
    for i, feature in enumerate(seq_record.features):
        if feature.type not in GENE_FEATURE_TYPES:
            continue
        start  = int(feature.location.start)
        end    = int(feature.location.end) - 1   # convert to inclusive
        strand = "+" if feature.location.strand >= 0 else "-"
        gene_id = (
            feature.qualifiers.get("gene",     [None])[0]
            or feature.qualifiers.get("locus_tag", [None])[0]
            or feature.qualifiers.get("sequence",  [None])[0]  # ORF-style label used by some tools
            or f"{feature.type}_{i}"
        )
        genes.append(Gene(gene_id=gene_id, start=start, end=end, strand=strand))
    return genes


def _parse_motifs(motif_df, sequence):
    """
    Build Motif objects from the motif table.

    The table must have a 'motif' column with the recognition sequence.
    'start' / 'end' columns are optional — if absent we search the sequence
    ourselves so every occurrence is covered.
    """
    motif_df = motif_df.rename(
        columns={column: str(column).strip().lower() for column in motif_df.columns}
    )

    motifs = []
    seen   = set()   # deduplicate (motif, start) pairs

    has_coords = (
        "start" in motif_df.columns
        and "end"   in motif_df.columns
    )
    has_motiffinder_coords = "position_1based" in motif_df.columns

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
            # No coordinates — search the sequence for all occurrences
            regex = compile_iupac(motif_seq)
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
                fwd_end   = len(sequence) - m.start() - 1
                fwd_start = len(sequence) - m.end()
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


def _parse_protected_regions(seq_record):
    """Extract ProtectedRegion objects from regulatory / misc_feature annotations."""
    protected = []
    PROTECTED_TYPES = {"regulatory", "misc_feature", "rep_origin", "promoter", "RBS"}
    for feature in seq_record.features:
        if feature.type in PROTECTED_TYPES and not _is_motiffinder_hit_marker(feature):
            start = int(feature.location.start)
            end   = int(feature.location.end) - 1
            protected.append(ProtectedRegion(start=start, end=end))
    return protected


def _parse_codon_usage(codon_df):
    """
    Convert a CodonBias DataFrame into the {codon: float} dict GenomeModel expects.

    Accepts column names 'fraction', 'frequency', or 'value' for the usage score.
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

    for _, row in codon_df.iterrows():
        codon = str(row["codon"]).strip().upper()
        try:
            usage[codon] = float(row[score_col])
        except (ValueError, TypeError):
            continue

    return usage


# ============================================================
# DECISION MATRIX ROW BUILDERS
# ============================================================

def _make_matrix_row(motif, candidate, score, chosen, genome):
    """Build one row of the decision matrix for a candidate edit."""
    gene = genome.find_gene(candidate.mutation.position)
    aa   = _translate(candidate.codon)

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
        # Codon change
        "original_codon":    candidate.codon,
        "replacement_codon": candidate.replacement,
        "amino_acid":        aa,
        "synonymous":        True,           # all candidates are synonymous by construction
        # Scoring breakdown — the columns the user cares about
        "motifs_destroyed":  candidate.result.get("destroyed", 0),
        "motifs_created":    candidate.result.get("created",   0),
        "usage_score":       round(candidate.usage_score, 6),
        "total_score":       round(score, 4),
        # Decision
        "chosen":            chosen,
        "skip_reason":       "",
    }


def _record_unresolvable(matrix, motif):
    """Record a sentinel row for a motif where no valid candidate existed."""
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
        "amino_acid":        "",
        "synonymous":        "",
        "motifs_destroyed":  0,
        "motifs_created":    0,
        "usage_score":       0,
        "total_score":       0,
        "chosen":            False,
        "skip_reason":       "no_valid_candidate",
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
