# build the genome model from the genbank file and codon usage table
# Some oop to encapsulate the genome and its features, and provide methods for editing and scoring
import os
import re
from enum import Enum
from tracemalloc import start
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
        genomic_codon = genome.sequence[
            codon_start:codon_start + 3]
        if len(genomic_codon) != 3:
            return None
        if self.strand == "+":
            return genomic_codon
        return reverse_complement(genomic_codon)
    def mutate_codon(self, genome, mutation):
        codon_start = self.codon_start(mutation.position)
        genomic_codon = genome.sequence[
            codon_start:codon_start + 3]
        if len(genomic_codon) != 3:
            raise ValueError("Invalid codon length")
        codon_list = list(genomic_codon)
        within_genomic = mutation.position - codon_start
        codon_list[within_genomic] = mutation.new
        mutated_genomic = "".join(codon_list)
        if self.strand == "+":
            return mutated_genomic
        return reverse_complement(mutated_genomic)
    def get_codon_start(self, position):
        return self.codon_start(position)
    def affected_codon_positions(self, mutation):
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
        try: aa = str(Bio.Seq.Seq(codon).translate())
        except Exception:
            return []
        return [
            c
            for c in SYNONYMOUS_CODONS.get(aa, [])
            if c != codon]
    def codon_mutations(self, codon_start, original_codon, replacement_codon):
        diffs = []
        for i in range(3):
            if original_codon[i] != replacement_codon[i]:
                if self.strand == "+":
                    genomic_position = codon_start + i
                    old_base = original_codon[i]
                    new_base = replacement_codon[i]
                else:
                    genomic_position = codon_start + i
                    genomic_original = reverse_complement(
                        original_codon)
                    genomic_replacement = reverse_complement(
                        replacement_codon)
                    old_base = genomic_original[i]
                    new_base = genomic_replacement[i]
                diffs.append(
                    Mutation(
                        position=genomic_position,
                        old=old_base,
                        new=new_base))
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

# Define circular and linear topology classes to handle sequence indexing and motif counting

class LinearTopology:
    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)

    def get_interval(self, start, end):
        start = max(0, start)
        end = min(self.length, end)
        return self.sequence[start:end]

class CircularTopology:
    def __init__(self, sequence):
        self.sequence = sequence
        self.length = len(sequence)

    def get_interval(self, start, end):
        L = self.length
        start %= L
        end %= L
        if start < end:
            return self.sequence[start:end]
        return (
            self.sequence[start:] +
            self.sequence[:end])
    
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
        self.topology_engine = self.build_topology(
            self.sequence)
        self.build_position_index()
    def build_topology(self, sequence):
        if self.topology == "circular":
            return CircularTopology(sequence)
        return LinearTopology(sequence)
    def ranked_synonymous_codons(self, codon, codon_usage):
        candidates = self.synonymous_codons(codon)
        return sorted(
            candidates,
            key=lambda c: codon_usage.get(c, 0),
            reverse=True)

    def generate_synonymous_candidates(self, motif):
        candidates = []
        genes = self.get_overlapping_genes(motif.start, motif.end)
        for gene in genes:
            codon_starts = gene.affected_codons(
                motif.start,
                motif.end)
            for codon_start in codon_starts:
                codon = gene.get_codon_by_start(
                    self,
                    codon_start)
                if codon is None:
                    continue
                synonyms = gene.ranked_synonymous_codons(
                    codon,
                    self.codon_usage)
                for synonym in synonyms:
                    mutations = gene.codon_mutations(
                        codon_start,
                        codon,
                        synonym)
                    for mutation in mutations:
                        result = self.evaluate_mutation(
                            mutation)
                        if result["valid"]:
                            candidates.append({
                                "mutation": mutation,
                                "result": result,
                                "codon": codon,
                                "replacement": synonym})
        return candidates

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
            for pos in range(
                motif.start,
                motif.end + 1):
                if pos < self.length:
                    self.position_index[pos].append(
                        motif)

    def get_local_motifs(self, start, end, radius=25):
        motif_set = set()
        window_start = max(0, start - radius)
        window_end = min(self.length - 1, end + radius)
        for pos in range(window_start, window_end + 1):
            motif_set.update(
                self.position_index.get(pos, []))
        return list(motif_set)
    
    def get_overlapping_genes(self, start, end):
        hits = []
        for gene in self.genes:
            if not (gene.end < start
                or
                gene.start > end):
                hits.append(gene)
        return hits

    # MUTATION HELPERS
    def apply_mutation(self, mutation):
        seq = list(self.sequence)
        expected = self.sequence[
            mutation.start:mutation.end + 1]
        if expected != mutation.old:
            raise ValueError(
                f"Mutation mismatch: expected {expected}, got {mutation.old}")
        seq[mutation.start:mutation.end + 1] = list(mutation.new)
        return "".join(seq)

    # EVALUATION ENGINE
    def evaluate_mutation(
        self,
        mutation,
        window_radius=25):
        region = self.get_region(
            mutation.position)
        # Regulatory regions are locked
        if region == RegionType.REGULATORY:
            return {
                "valid": False,
                "reason": "Protected region"}
        # CDS synonymous check
        if region == RegionType.CDS:
            if not is_synonymous(
                self,
                mutation):
                return {
                    "valid": False,
                    "reason": "Not synonymous"}

        # Build mutated sequence and topology
        mutated_sequence = self.apply_mutation(mutation)
        mutated_topology = self.build_topology(mutated_sequence)

        # Build windows around the mutation for motif checking
        original_window = (
            self.topology_engine.get_interval(
                mutation.start - window_radius,
                mutation.end + window_radius))
        mutated_window = (
            mutated_topology.get_interval(
                mutation.start - window_radius,
                mutation.end + window_radius))
        original_rc = reverse_complement(
            original_window)
        mutated_rc = reverse_complement(
            mutated_window)
        
        # Local motifs only
        local_motifs = self.get_local_motifs(
            mutation.start,
            mutation.end,
            radius=window_radius)
        destroyed = 0
        created = 0

        for motif in local_motifs:
            before = (
                motif.regex.search(
                    original_window
                ) is not None
                or
                motif.regex.search(
                    original_rc
                ) is not None)

            after = (
                motif.regex.search(
                    mutated_window
                ) is not None
                or
                motif.regex.search(
                    mutated_rc
                ) is not None)

            if before and not after:
                destroyed += 1

            elif not before and after:
                created += 1

            if created > 0:
                return {
                    "valid": False,
                    "reason": "Creates new motif"
                    }
        return {
            "valid": True,
            "destroyed": destroyed,
            "created": created,
            "edits": 1,
            "score": destroyed - (created * 10)}

# UTILITY FUNCTIONS
def is_synonymous(genome, mutation):
    gene = genome.find_gene(
        mutation.position)
    if gene is None:
        return False
    original_codon = gene.get_codon(
        genome,
        mutation.position)
    mutated_codon = gene.mutate_codon(
        genome,
        mutation)
    if (original_codon is None
        or
        mutated_codon is None):
        return False
    return (
        Bio.Seq.Seq(original_codon).translate()
        ==
        Bio.Seq.Seq(mutated_codon).translate()
    )