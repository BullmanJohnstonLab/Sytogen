DEBUG = False

def debug(msg):
    if DEBUG: print(msg)

# build the genome model from the genbank file and codon usage table
# Some oop to encapsulate the genome and its features, and provide methods for editing and scoring
import re
from enum import Enum
import Bio.Seq
from sytogen.scripts.legacy_sytogen import reverse_complement
from Bio.Data import CodonTable
from collections import defaultdict

STANDARD_TABLE = CodonTable.unambiguous_dna_by_name["Standard"]
SYNONYMOUS_CODONS = defaultdict(list)
for codon, amino_acid in STANDARD_TABLE.forward_table.items():
    SYNONYMOUS_CODONS[amino_acid].append(codon)
for stop_codon in STANDARD_TABLE.stop_codons:
    SYNONYMOUS_CODONS["*"].append(stop_codon)


# ============================================================
# IUPAC SUPPORT
# ============================================================

IUPAC_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]"}

def compile_iupac(motif):
    pattern: str = "".join(IUPAC_MAP[b] for b in motif.upper())
    return re.compile(pattern)


class RegionType(Enum):
    CDS = "CDS"
    REGULATORY = "REGULATORY"
    NEUTRAL = "NEUTRAL"


class Gene:
    def __init__(self, gene_id, start, end, strand="+"):
        self.id = gene_id
        self.start = start
        self.end = end          # raw / un-clamped — end >= genome_length
                                 # signals a gene that spans the circular
                                 # origin (see _parse_genes), same
                                 # convention already used for Motif/Overlap
        self.strand = strand

    def _resolve(self, position, genome_length):
        """
        Map a normalized position (always given in [0, genome_length)) into
        this gene's own coordinate frame. For a gene that wraps the origin
        (self.end >= genome_length), a position numerically "before"
        self.start actually belongs to the wrapped-around tail of the gene
        (the part that continues past position 0) and needs genome_length
        added to line up with self.end. A no-op for a normal gene, or when
        genome_length isn't known.
        """
        if genome_length and self.end >= genome_length and position < self.start:
            return position + genome_length
        return position

    def contains(self, position, genome_length=None):
        pos = self._resolve(position, genome_length)
        return self.start <= pos <= self.end

    def codon_start(self, position, genome_length=None):
        pos = self._resolve(position, genome_length)
        if self.strand == "+":
            offset = pos - self.start
            return self.start + (offset // 3) * 3
        else:
            offset = self.end - pos
            return self.end - ((offset // 3) + 1) * 3 + 1
        # The returned codon_start is itself left raw/un-clamped — for the
        # rare codon that straddles the origin it can be negative or
        # >= genome_length. Callers read the 3 bases via
        # topology_engine.get_interval(), which already resolves that
        # wraparound (modulo) the same way it does for Motif windows.

    def get_codon(self, genome, position):
        codon_start = self.codon_start(position, genome.length)
        genomic_codon = genome.topology_engine.get_interval(codon_start, codon_start + 3)
        debug(f"[get_codon]"
              f"pos={position}A"
              f"start={codon_start}"
              f"raw={genomic_codon}"
              f"strand={self.strand}")
        if len(genomic_codon) != 3:
            debug(f"[get_codon] INVALID codon length at {codon_start}")
            return genomic_codon
        if self.strand == "+":
            return genomic_codon
        rc = reverse_complement(genomic_codon)
        debug(f"[get_codon] reverse complement -> {rc}")
        return rc

    def mutate_codon(self, genome, mutation):
        codon_start = self.codon_start(mutation.position, genome.length)
        genomic_codon = genome.topology_engine.get_interval(codon_start, codon_start + 3)
        debug(f"[mutate_codon]"
              f"original={genomic_codon} at {codon_start} "
              f"mutation={mutation.position}:{mutation.old}->{mutation.new}")
        if len(genomic_codon) != 3:
            raise ValueError("Invalid codon length")
        codon_list = list(genomic_codon)
        # mutation.position is always a normalized, valid genomic index
        # (codon_mutations() already resolved it that way); codon_start
        # may be raw/negative/>=length for a straddling codon, so index
        # by the normalized *distance* between the two rather than
        # subtracting raw values directly.
        within_genomic = (mutation.position - codon_start) % genome.length
        codon_list[within_genomic] = mutation.new
        mutated_genomic = "".join(codon_list)
        debug(f"[mutate_codon] mutated genomic={mutated_genomic}")
        if self.strand == "+":
            return mutated_genomic
        rc = reverse_complement(mutated_genomic)
        debug(f"[mutate_codon] reverse complement mutated={rc}")
        return rc

    def get_codon_start(self, position):
        return self.codon_start(position)

    def affected_codons(self, start, end):
        codon_start = set()
        interval_start = max(start, self.start)
        interval_end = min(end, self.end)
        for pos in range(interval_start, interval_end + 1):
            codon_start.add(self.codon_start(pos))
        return sorted(codon_start)

    def get_codon_by_start(self, genome, codon_start):
        genomic_codon = genome.topology_engine.get_interval(codon_start, codon_start + 3)
        if len(genomic_codon) != 3:
            return None
        if self.strand == "+":
            return genomic_codon
        return reverse_complement(genomic_codon)

    def synonymous_codons(self, codon):
        codon = codon.upper()
        try:
            aa = str(Bio.Seq.Seq(codon).translate())
        except Exception:
            debug(f"[synonymous_codons] failed translation for {codon}")
            return []
        syn = [c for c in SYNONYMOUS_CODONS.get(aa, []) if c != codon]
        debug(f"[synonymous_codons] codon={codon} aa={aa} synonyms={syn}")
        return syn

    def ranked_synonymous_codons(self, codon, codon_usage):
        candidates = self.synonymous_codons(codon)
        return sorted(
            candidates,
            key=lambda c: codon_usage.get(c, 0),
            reverse=True)

    def codon_mutations(self, codon_start, original_codon, replacement_codon, genome_length=None):
        diffs = []
        debug(f"[codon_mutations] {original_codon} -> {replacement_codon} at {codon_start}")
        for i in range(3):
            if original_codon[i] != replacement_codon[i]:
                if self.strand == "+":
                    genomic_position = codon_start + i
                    old_base = original_codon[i]
                    new_base = replacement_codon[i]
                else:
                    genomic_position = codon_start + (2 - i)
                    genomic_original = reverse_complement(original_codon)
                    genomic_replacement = reverse_complement(replacement_codon)
                    # FIX: genomic_position uses the mirrored index (2 - i)
                    # since coding-frame base i sits at genomic position
                    # codon_start + (2 - i) on a minus-strand gene. old_base/
                    # new_base must be read from that same mirrored index,
                    # not i, or they describe the base at a different genomic
                    # position than the one just computed — causing
                    # apply_mutation()'s sequence check to fail with
                    # "Mutation mismatch: expected X, got Y" whenever the
                    # differing base isn't the middle codon position.
                    old_base = genomic_original[2 - i]
                    new_base = genomic_replacement[2 - i]
                    debug(f"[codon_mutations] diff at pos={genomic_position} {old_base}->{new_base}")
                # codon_start can be raw (negative or >= genome_length) for
                # the rare codon that straddles the circular origin — see
                # Gene.codon_start(). A single base's own position is
                # always meaningful once wrapped back into [0, genome_length),
                # which is all apply_mutation() needs; it never has to
                # handle a wrapping *range* itself, just this one valid index.
                if genome_length:
                    genomic_position %= genome_length
                diffs.append(
                    Mutation(
                        position=genomic_position,
                        old=old_base,
                        new=new_base))
        debug(f"[codon_mutations] total diffs={len(diffs)}")
        if len(diffs) == 1:
            return diffs
        return []


class ProtectedRegion:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def contains(self, pos):
        return self.start <= pos <= self.end


class Motif:
    def __init__(self, motif, start, end, strand="+"):
        self.motif = motif
        self.start = start
        self.end = end
        self.strand = strand
        self.length = len(motif)
        self.regex = compile_iupac(motif)

    def overlaps(self, start, end):
        return not (
            self.end < start or
            self.start > end)


class Mutation:
    def __init__(self, position, old, new):
        self.position = position
        self.old = old
        self.new = new
        self.start = position
        self.end = position + len(old) - 1


GC_BASES = frozenset("GC")
AT_BASES = frozenset("AT")


def is_gc_preserving_swap(old_base, new_base):
    """
    True if old_base -> new_base stays within the same GC-class (a G<->C
    or A<->T swap) rather than crossing between them. A same-class swap
    never changes local GC content at all.

    This ranks BELOW codon-usage preference (see run_sytogen_pipeline's
    candidate sort in sytogen_runner.py) — it's a tiebreaker for candidates
    that are otherwise equally good, not a reason to pick a worse-usage
    codon over a better one.
    """
    return (
        (old_base in GC_BASES and new_base in GC_BASES)
        or (old_base in AT_BASES and new_base in AT_BASES)
    )


class Candidate:
    def __init__(self, mutation, result, codon, replacement, usage_score=0):
        self.mutation = mutation
        self.result = result
        self.codon = codon
        self.replacement = replacement
        self.usage_score = usage_score


# Define circular and linear topology classes to handle sequence indexing and motif counting

class LinearTopology:
    name = "linear"

    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)

    def get_interval(self, start, end):
        """Return the subsequence from start to end, clamped to sequence bounds."""
        start = max(0, start)
        end = min(self.length, end)
        return self.sequence[start:end]

    def normalize_position(self, pos):
        """Clamp pos to [0, length-1]. Out-of-range positions are invalid on a linear molecule."""
        return max(0, min(self.length - 1, pos))

    def count_motif_hits(self, regex):
        """Count non-overlapping forward-strand hits across the linear sequence."""
        return len(regex.findall(self.sequence))


class CircularTopology:
    name = "circular"

    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)

    def get_interval(self, start, end):
        """Return the subsequence from start to end, wrapping around the origin if needed.
        Both start and end are taken modulo length so out-of-range indices always resolve."""
        L = self.length
        start %= L
        end %= L
        if start < end:
            return self.sequence[start:end]
        # Wraps around the origin
        return self.sequence[start:] + self.sequence[:end]

    def normalize_position(self, pos):
        """Wrap pos into [0, length-1] — all positions are valid on a circular molecule."""
        return pos % self.length

    def count_motif_hits(self, regex):
        """Count non-overlapping forward-strand hits, including those that span the origin.
        We scan a doubled sequence and deduplicate hits that fall in the first copy."""
        doubled = self.sequence + self.sequence
        hits = [m.start() for m in regex.finditer(doubled)]
        # Keep only hits whose start falls in the first copy
        return len([h for h in hits if h < self.length])


# ============================================================
# GENOME MODEL
# ============================================================

class GenomeModel:
    def __init__(
        self,
        sequence,
        topology="circular",
        genes=None,
        motifs=None,
        protected_regions=None,
        codon_usage=None):

        self.sequence = sequence.upper()
        self.length = len(sequence)
        self.topology = topology
        self.genes = genes or []
        self.motifs = motifs or []
        self.protected_regions = protected_regions or []
        self.codon_usage = codon_usage or {}
        self.position_index = {}
        self.topology_engine = self.build_topology(self.sequence)
        self.build_position_index()

    def build_topology(self, sequence):
        if self.topology == "circular":
            return CircularTopology(sequence)
        return LinearTopology(sequence)

    def set_topology(self, topology):
        """Toggle between 'circular' and 'linear' with a full re-parse from scratch.
        Rebuilds the topology engine, resets the position index, and updates length."""
        if topology not in ("circular", "linear"):
            raise ValueError(f"topology must be 'circular' or 'linear', got {topology!r}")
        if topology == self.topology:
            debug(f"[set_topology] already {topology}, no-op")
            return
        debug(f"[set_topology] switching {self.topology} → {topology}")
        self.topology = topology
        self.length = len(self.sequence)          # re-derive in case sequence was mutated
        self.topology_engine = self.build_topology(self.sequence)
        self.build_position_index()               # full reindex under new topology
        debug(f"[set_topology] done — engine={self.topology_engine.name}")

    def generate_synonymous_candidates(self, motif):
        candidates = []
        # A motif almost always spans more than one base of the same codon
        # (a 6bp restriction site inside a gene typically covers 2+
        # positions of the same 3bp codon). Without this, every one of
        # those positions independently re-expands the identical codon's
        # synonymous options, producing exact duplicate candidates/rows —
        # same mutation position, same old/new base — for every extra
        # position that happens to fall in that codon. Track codons
        # already expanded for this motif and only do it once each.
        seen_codons = set()
        for raw_pos in range(motif.start, motif.end + 1):
            # raw_pos may fall outside [0, length) when motif.end was left
            # un-clamped to represent a motif that spans the circular
            # origin (see _parse_motifs / Motif). normalize_position()
            # wraps it correctly for circular topology (and is a no-op —
            # positions here are always already in range — for linear).
            pos = self.topology_engine.normalize_position(raw_pos)
            debug(f"\n[position] {pos}")
            gene = self.find_gene(pos)
            # FIX: these lines were un-indented out of the for-loop body
            if gene is None:
                debug(f"[position] {pos} not in gene → skipping")
                continue
            debug(f"[position] {pos} in gene {gene.id} strand={gene.strand}")
            codon_start = gene.codon_start(pos, self.length)
            codon_key = (gene.id, codon_start)
            if codon_key in seen_codons:
                debug(f"[codon] {codon_key} already expanded for this motif → skipping duplicate")
                continue
            seen_codons.add(codon_key)
            codon = gene.get_codon(self, pos)
            if codon is None:
                debug(f"[codon] None at pos {pos} → skipping")
                continue
            debug(f"[codon] original codon={codon} start={codon_start}")
            synonymous = gene.synonymous_codons(codon)
            if not synonymous:
                debug(f"[synonymous] no alternatives for {codon}")
                continue
            debug(f"[synonymous] candidates={synonymous}")
            for replacement in synonymous:
                debug(f"\n[replacement] trying {codon} -> {replacement}")
                mutations = gene.codon_mutations(codon_start, codon, replacement, self.length)
                if not mutations:
                    debug(f"[replacement] rejected (multi-base or invalid)")
                    continue
                mutation = mutations[0]
                debug(f"[mutation] pos={mutation.position} {mutation.old}->{mutation.new}")

                result = self.evaluate_mutation(mutation)
                debug(f"[simulate] result={result}")

                if result is None or not result.get("valid"):
                    debug(f"[simulate] returned None or invalid → skipping")
                    continue
                usage_score = self.codon_usage.get(replacement, 0)
                debug(f"[usage] replacement={replacement} usage_score={usage_score}")
                candidate = Candidate(
                    mutation=mutation,
                    result=result,
                    codon=codon,
                    replacement=replacement,
                    usage_score=usage_score)
                debug(f"[candidate] ACCEPTED -> destroyed={result.get('destroyed')} edits={result.get('edits')}")
                candidates.append(candidate)
        debug(f"\n[generate_candidates] TOTAL candidates={len(candidates)}")
        return candidates

    def generate_neutral_candidates(self, motif):
        """
        For motif positions that fall outside any annotated gene AND outside
        any protected region (RegionType.NEUTRAL), there's no codon/amino-acid
        constraint to preserve — the position isn't coding for anything, so a
        single-base substitution is fine on its own merits. The only
        requirement is the same one evaluate_mutation() already enforces for
        every edit: it must not create a new motif site elsewhere in the
        local window (and if the substitution window happens to touch a
        protected region, evaluate_mutation() will still reject it there —
        this method doesn't bypass that, it just doesn't require synonymy).

        Only substitutions that actually destroy this motif are kept —
        undirected single-base changes that don't help silence anything
        aren't useful candidates.
        """
        candidates = []
        for raw_pos in range(motif.start, motif.end + 1):
            pos = self.topology_engine.normalize_position(raw_pos)
            if self.find_gene(pos) is not None:
                continue  # coding position — handled by generate_synonymous_candidates
            if self.get_region(pos) != RegionType.NEUTRAL:
                continue  # protected — leave untouched

            original_base = self.sequence[pos]
            for new_base in "ACGT":
                if new_base == original_base:
                    continue
                mutation = Mutation(position=pos, old=original_base, new=new_base)
                result = self.evaluate_mutation(mutation)
                if result is None or not result.get("valid"):
                    continue
                if result.get("destroyed", 0) <= 0:
                    continue  # only care about edits that actually silence this motif

                candidate = Candidate(
                    mutation=mutation,
                    result=result,
                    codon=original_base,
                    replacement=new_base,
                    usage_score=0)  # no codon-usage concept outside a gene
                candidates.append(candidate)
        return candidates

    def explain_no_candidates(self, motif):
        """
        Diagnostic companion to generate_synonymous_candidates() +
        generate_neutral_candidates(). Walks the same positions and gates
        both paths use, but instead of silently discarding a motif that
        produced zero candidates, returns why — as structured data
        (attempted/rejected counts, most common rejection reason) rather
        than only a prose sentence, so the decision matrix can carry it in
        dedicated columns instead of parsing it back out of text.
        """
        saw_gene_position = False
        saw_editable_gene_position = False
        saw_synonymous_alternative = False
        saw_neutral_position = False
        gene_rejection_reasons = []
        neutral_rejection_reasons = []

        seen_codons = set()
        for raw_pos in range(motif.start, motif.end + 1):
            pos = self.topology_engine.normalize_position(raw_pos)
            gene = self.find_gene(pos)

            if gene is not None:
                saw_gene_position = True
                if self.get_region(pos) == RegionType.REGULATORY:
                    continue
                saw_editable_gene_position = True

                codon_start = gene.codon_start(pos, self.length)
                codon_key = (gene.id, codon_start)
                if codon_key in seen_codons:
                    continue  # already evaluated this codon via another motif position
                seen_codons.add(codon_key)

                codon = gene.get_codon(self, pos)
                if not codon or len(codon) != 3:
                    continue
                synonymous = gene.synonymous_codons(codon)
                if not synonymous:
                    continue
                saw_synonymous_alternative = True

                for replacement in synonymous:
                    mutations = gene.codon_mutations(codon_start, codon, replacement, self.length)
                    if not mutations:
                        continue
                    result = self.evaluate_mutation(mutations[0])
                    if result and not result.get("valid"):
                        gene_rejection_reasons.append(result.get("reason", "invalid"))

            elif self.get_region(pos) == RegionType.NEUTRAL:
                saw_neutral_position = True
                original_base = self.sequence[pos]
                for new_base in "ACGT":
                    if new_base == original_base:
                        continue
                    mutation = Mutation(position=pos, old=original_base, new=new_base)
                    result = self.evaluate_mutation(mutation)
                    if result is None or not result.get("valid"):
                        reason = (result or {}).get("reason", "invalid")
                        neutral_rejection_reasons.append(reason)
                    elif result.get("destroyed", 0) <= 0:
                        neutral_rejection_reasons.append("does not destroy this motif")

        def _no_attempt(reason_code, reasoning):
            return {
                "reason_code": reason_code,
                "reasoning": reasoning,
                "attempted_count": None,
                "rejected_count": None,
                "top_rejection_reason": None,
                "top_rejection_count": None,
            }

        # Nothing editable at all — every position is either protected, or a
        # gene position whose region check still failed for some other reason.
        if not saw_gene_position and not saw_neutral_position:
            return _no_attempt(
                "blocked_by_protected_region",
                "Every position this motif spans falls inside a protected "
                "regulatory annotation, so no edit is allowed anywhere in it.")

        from collections import Counter

        # Prefer explaining the coding-side outcome if there was a gene here,
        # since that's usually the more informative story (protection vs.
        # amino-acid constraints vs. rejected edits).
        if saw_gene_position:
            if not saw_editable_gene_position and not saw_neutral_position:
                return _no_attempt(
                    "blocked_by_protected_region",
                    "Every gene position this motif overlaps falls inside a "
                    "protected regulatory annotation, so no edit is allowed here.")
            if saw_editable_gene_position and not saw_synonymous_alternative and not saw_neutral_position:
                return _no_attempt(
                    "no_synonymous_codon",
                    "The amino acid encoded here has no synonymous codon "
                    "alternative (e.g. Met/Trp), so no silent edit is possible.")
            if gene_rejection_reasons and not saw_neutral_position:
                top_reason, top_count = Counter(gene_rejection_reasons).most_common(1)[0]
                return {
                    "reason_code": "all_candidates_rejected",
                    "reasoning": "Every synonymous codon option at this position was rejected "
                                 "(see attempted/rejected columns for the tally).",
                    "attempted_count": len(gene_rejection_reasons),
                    "rejected_count": len(gene_rejection_reasons),
                    "top_rejection_reason": top_reason,
                    "top_rejection_count": top_count,
                }

        if saw_neutral_position:
            if neutral_rejection_reasons:
                top_reason, top_count = Counter(neutral_rejection_reasons).most_common(1)[0]
                return {
                    "reason_code": "all_candidates_rejected",
                    "reasoning": "Every non-coding substitution attempted at this position was "
                                 "rejected (see attempted/rejected columns for the tally).",
                    "attempted_count": len(neutral_rejection_reasons),
                    "rejected_count": len(neutral_rejection_reasons),
                    "top_rejection_reason": top_reason,
                    "top_rejection_count": top_count,
                }

        return _no_attempt(
            "no_valid_edit",
            "No valid single-base substitution could be constructed here.")

    def score_candidate(self, candidate):
        """
        The single source of truth for how a candidate edit is ranked —
        this is what run_sytogen_pipeline actually sorts candidates by.
        Codon choice affects the score entirely through usage_score
        (looked up from the codon-usage table when the candidate was
        built in generate_synonymous_candidates/generate_neutral_candidates),
        so a different codon choice for the same position naturally scores
        differently here without needing any separate mechanism for it.

        GC-class preservation (is_gc_preserving_swap) is deliberately NOT
        folded into this score — it's a lower priority than codon usage
        preference, so it should only decide between candidates that are
        otherwise exactly tied, not outrank a real usage_score difference.
        run_sytogen_pipeline applies it as a secondary sort key on top of
        this score for exactly that reason.
        """
        score = 0
        # Prioritize candidates that destroy more motifs
        score += (candidate.result["destroyed"] * 1000)
        # codon preference bonus
        score += (candidate.usage_score * 100)
        # Penalize edits
        score -= (candidate.result["edits"] * 10)
        debug(f"[score_candidate] "
              f"destroyed={candidate.result['destroyed']} "
              f"usage={candidate.usage_score} "
              f"edits={candidate.result['edits']} "
              f"score={score}")
        return score

    # REGION LOOKUPS
    def get_region(self, pos):
        for region in self.protected_regions:
            if region.contains(pos):
                return RegionType.REGULATORY
        for gene in self.genes:
            if gene.contains(pos, self.length):
                return RegionType.CDS
        return RegionType.NEUTRAL

    def is_protected(self, pos):
        for region in self.protected_regions:
            if region.contains(pos):
                return True
        return False

    def find_gene(self, position):
        for gene in self.genes:
            if gene.contains(position, self.length):
                return gene
        return None

    # MOTIF INDEX
    def build_position_index(self):
        self.position_index = {
            i: [] for i in range(self.length)}
        for motif in self.motifs:
            for raw_pos in range(motif.start, motif.end + 1):
                # motif.end can exceed self.length - 1 for a motif that
                # spans the circular origin — normalize instead of
                # dropping those positions, or the wrapped-around tail of
                # the motif would never get indexed.
                pos = self.topology_engine.normalize_position(raw_pos)
                self.position_index[pos].append(motif)

    def get_local_motifs(self, start, end, radius=25):
        motif_set = set()
        for raw_pos in range(start - radius, end + radius + 1):
            # Same reasoning as build_position_index: wrap the window
            # around the origin for circular topology instead of clamping
            # it away, so a mutation near one side of the origin still
            # sees motifs indexed just past it on the other side.
            pos = self.topology_engine.normalize_position(raw_pos)
            motif_set.update(self.position_index.get(pos, []))
        return list(motif_set)

    def get_overlapping_genes(self, start, end):
        hits = []
        for gene in self.genes:
            if gene.end >= self.length:
                # Wraps the origin — its true span is two real segments,
                # not one contiguous [start, end]. Check both.
                seg1 = (gene.start, self.length - 1)
                seg2 = (0, gene.end - self.length)
                overlaps = (
                    not (seg1[1] < start or seg1[0] > end)
                    or not (seg2[1] < start or seg2[0] > end)
                )
                if overlaps:
                    hits.append(gene)
            elif not (gene.end < start or gene.start > end):
                hits.append(gene)
        return hits

    # MUTATION HELPERS
    def apply_mutation(self, mutation):
        seq = list(self.sequence)
        expected = self.sequence[mutation.start:mutation.end + 1]
        if expected != mutation.old:
            raise ValueError(
                f"Mutation mismatch: expected {expected}, got {mutation.old}")
        seq[mutation.start:mutation.end + 1] = list(mutation.new)
        return "".join(seq)

    def debug_window(self, pos, window=10):
        start = max(0, pos - window)
        end = pos + window
        debug(f"[window] {start}:{end} -> {self.sequence[start:end]}")

    # EVALUATION ENGINE
    def evaluate_mutation(self, mutation, window_radius=25):
        region = self.get_region(mutation.position)
        # Regulatory regions are locked
        if region == RegionType.REGULATORY:
            return {
                "valid": False,
                "reason": "Protected region"}
        # CDS synonymous check
        if region == RegionType.CDS:
            if not is_synonymous(self, mutation):
                return {
                    "valid": False,
                    "reason": "Not synonymous"}

        # Build mutated sequence and topology
        mutated_sequence = self.apply_mutation(mutation)
        mutated_topology = self.build_topology(mutated_sequence)

        # Build windows around the mutation for motif checking
        original_window = self.topology_engine.get_interval(
            mutation.start - window_radius,
            mutation.end + window_radius)
        mutated_window = mutated_topology.get_interval(
            mutation.start - window_radius,
            mutation.end + window_radius)
        original_rc = reverse_complement(original_window)
        mutated_rc = reverse_complement(mutated_window)

        # Local motifs only
        local_motifs = self.get_local_motifs(
            mutation.start,
            mutation.end,
            radius=window_radius)
        destroyed = 0
        created = 0

        for motif in local_motifs:
            before = (
                motif.regex.search(original_window) is not None
                or motif.regex.search(original_rc) is not None)

            after = (
                motif.regex.search(mutated_window) is not None
                or motif.regex.search(mutated_rc) is not None)

            if before and not after:
                destroyed += 1
            elif not before and after:
                created += 1

            if created > 0:
                return {
                    "valid": False,
                    "reason": "Creates new motif"}

        return {
            "valid": True,
            "destroyed": destroyed,
            "created": created,
            "edits": 1,
            "score": destroyed - (created * 10)}


# ============================================================
# UTILITY FUNCTIONS (module-level, not methods)
# ============================================================

def is_synonymous(genome, mutation):
    gene = genome.find_gene(mutation.position)
    if gene is None:
        return False
    original_codon = gene.get_codon(genome, mutation.position)
    mutated_codon = gene.mutate_codon(genome, mutation)
    if original_codon is None or mutated_codon is None:
        return False
    return (Bio.Seq.Seq(original_codon).translate()
            == Bio.Seq.Seq(mutated_codon).translate())


def motif_destroyed(genome, motif):
    """FIX: was ignoring the reverse-complement check entirely."""
    window = genome.topology_engine.get_interval(
        motif.start,
        motif.end + motif.length)
    rc = reverse_complement(window)
    return (
        motif.regex.search(window) is None
        and motif.regex.search(rc) is None)
