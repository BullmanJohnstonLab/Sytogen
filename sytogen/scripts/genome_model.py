DEBUG = True

def debug(msg):
    if DEBUG: print(msg)

# build the genome model from the genbank file and codon usage table
# Some oop to encapsulate the genome and its features, and provide methods for editing and scoring
import copy
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
        self.end = end
        self.strand = strand

    def contains(self, pos):
        return self.start <= pos <= self.end

    def codon_start(self, position):
        if self.strand == "+":
            offset = position - self.start
            return self.start + (offset // 3) * 3
        else:
            offset = self.end - position
            return self.end - ((offset // 3) + 1) * 3 + 1

    def get_codon(self, genome, position):
        codon_start = self.codon_start(position)
        genomic_codon = genome.sequence[codon_start:codon_start + 3]
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
        codon_start = self.codon_start(mutation.position)
        genomic_codon = genome.sequence[codon_start:codon_start + 3]
        debug(f"[mutate_codon]"
              f"original={genomic_codon} at {codon_start} "
              f"mutation={mutation.position}:{mutation.old}->{mutation.new}")
        if len(genomic_codon) != 3:
            raise ValueError("Invalid codon length")
        codon_list = list(genomic_codon)
        within_genomic = mutation.position - codon_start
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
        genomic_codon = genome.sequence[codon_start:codon_start + 3]
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

    def codon_mutations(self, codon_start, original_codon, replacement_codon):
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
                    old_base = genomic_original[i]
                    new_base = genomic_replacement[i]
                    debug(f"[codon_mutations] diff at pos={genomic_position} {old_base}->{new_base}")
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


class Candidate:
    def __init__(self, mutation, result, codon, replacement, usage_score=0):
        self.mutation = mutation
        self.result = result
        self.codon = codon
        self.replacement = replacement
        self.usage_score = usage_score

    @property
    def score(self):
        return (self.result["destroyed"] * 100
                + self.usage_score)


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
        for pos in range(motif.start, motif.end + 1):
            debug(f"\n[position] {pos}")
            gene = self.find_gene(pos)
            # FIX: these lines were un-indented out of the for-loop body
            if gene is None:
                debug(f"[position] {pos} not in gene → skipping")
                continue
            debug(f"[position] {pos} in gene {gene.id} strand={gene.strand}")
            codon = gene.get_codon(self, pos)
            if codon is None:
                debug(f"[codon] None at pos {pos} → skipping")
                continue
            debug(f"[codon] original codon={codon}")
            codon_start = gene.codon_start(pos)
            debug(f"[codon] start={codon_start}")
            synonymous = gene.synonymous_codons(codon)
            if not synonymous:
                debug(f"[synonymous] no alternatives for {codon}")
                continue
            debug(f"[synonymous] candidates={synonymous}")
            for replacement in synonymous:
                debug(f"\n[replacement] trying {codon} -> {replacement}")
                mutations = gene.codon_mutations(codon_start, codon, replacement)
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

    def score_candidate(self, candidate):
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

    def best_candidate(self, motif):
        candidates = self.generate_synonymous_candidates(motif)
        debug(f"[best_candidate] {len(candidates)} candidates generated")
        if not candidates:
            return None
        best = max(candidates, key=self.score_candidate)
        debug(f"[best_candidate] SELECTED codon={best.codon} -> {best.replacement}")
        return best

    def optimize_motif(self, motif, max_iterations=10):
        edits = []
        for i in range(max_iterations):
            debug(f"\n[optimize] iteration {i}")
            candidate = self.best_candidate(motif)
            if candidate is None:
                debug("[optimize] no candidate found, stopping")
                break
            edits.append(candidate)
            debug(f"[optimize] applying mutation at {candidate.mutation.position}")
            self.sequence = self.apply_mutation(candidate.mutation)
            self.topology_engine = self.build_topology(self.sequence)
            if motif_destroyed(self, motif):
                debug("[optimize] motif destroyed, stopping")
                break
        return edits

    def lookahead_score(self, motif, depth):
        candidates = self.generate_synonymous_candidates(motif)
        if not candidates:
            return 0
        if depth == 0:
            return max(self.score_candidate(c) for c in candidates)
        best = float("-inf")
        for candidate in candidates:
            immediate = self.score_candidate(candidate)
            future_genome = self.simulate_candidate(candidate)
            future = future_genome.lookahead_score(motif, depth - 1)
            total = immediate + future
            best = max(best, total)
        return best

    # REGION LOOKUPS
    def get_region(self, pos):
        for region in self.protected_regions:
            if region.contains(pos):
                return RegionType.REGULATORY
        for gene in self.genes:
            if gene.contains(pos):
                return RegionType.CDS
        return RegionType.NEUTRAL

    def is_protected(self, pos):
        for region in self.protected_regions:
            if region.contains(pos):
                return True
        return False

    def find_gene(self, position):
        for gene in self.genes:
            if gene.contains(position):
                return gene
        return None

    # MOTIF INDEX
    def build_position_index(self):
        self.position_index = {
            i: [] for i in range(self.length)}
        for motif in self.motifs:
            for pos in range(motif.start, motif.end + 1):
                if pos < self.length:
                    self.position_index[pos].append(motif)

    def get_local_motifs(self, start, end, radius=25):
        motif_set = set()
        window_start = max(0, start - radius)
        window_end = min(self.length - 1, end + radius)
        for pos in range(window_start, window_end + 1):
            motif_set.update(self.position_index.get(pos, []))
        return list(motif_set)

    def get_overlapping_genes(self, start, end):
        hits = []
        for gene in self.genes:
            if not (gene.end < start or gene.start > end):
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

    def apply_best_candidate(self, motif):
        candidate = self.best_candidate(motif)
        if candidate is None:
            return None
        self.sequence = self.apply_mutation(candidate.mutation)
        self.topology_engine = self.build_topology(self.sequence)
        return candidate

    def clone(self):
        return copy.deepcopy(self)

    def simulate_candidate(self, candidate):
        temp = self.clone()
        temp.sequence = temp.apply_mutation(candidate.mutation)
        temp.topology_engine = temp.build_topology(temp.sequence)
        return temp

    def lookahead_best_candidate(self, motif, depth=2):
        candidates = self.generate_synonymous_candidates(motif)
        if not candidates:
            return None
        best_candidate = None
        best_score = float("-inf")
        for candidate in candidates:
            immediate = self.score_candidate(candidate)
            future_genome = self.simulate_candidate(candidate)
            future = future_genome.lookahead_score(motif, depth - 1)
            total = immediate + future
            if total > best_score:
                best_score = total
                best_candidate = candidate
        return best_candidate

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
