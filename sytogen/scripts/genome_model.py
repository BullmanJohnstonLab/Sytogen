# build the genome model from the genbank file and codon usage table
# Some oop to encapsulate the genome and its features, and provide methods for editing and scoring
import os
import re
from enum import Enum

from sytogen.scripts.legacy_sytogen import reverse_complement


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
    pattern = "".join(IUPAC_MAP[b] for b in motif.upper())
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


class ProtectedRegion:

    def __init__(self, name, start, end):
        self.name = name
        self.start = start
        self.end = end

    def contains(self, pos):
        return self.start <= pos <= self.end

    def __init__(
        self,
        topology="circular",
        sequence=None,
        genes=None,
        motifs=None,
        protected_regions=None,
        codon_usage=None):

        self.sequence = sequence
        self.length = len(sequence)
        self.topology = topology  # "circular" or "linear"
        self.genes = genes or []
        self.motifs = motifs or []
        self.protected_regions = protected_regions or []
        self.codon_usage = codon_usage or {}

        self.topology_engine = (
            CircularTopology(sequence)
            if topology == "circular"
            else LinearTopology(sequence))

class Gene:
    def __init__(self, id, start, end, strand):
        self.id = id
        self.start = start
        self.end = end
        self.strand = strand

class Motif:
    def __init__(
        self,
        motif,
        start,
        end,
        strand="+"):

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

    # --------------------------------------------------
    # REGION LOOKUPS
    # --------------------------------------------------

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

    # --------------------------------------------------
    # MOTIF INDEX
    # --------------------------------------------------

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

    def get_local_motifs(
        self,
        start,
        end,
        radius=25
    ):

        motif_set = set()

        window_start = max(
            0,
            start - radius
        )

        window_end = min(
            self.length - 1,
            end + radius
        )

        for pos in range(
            window_start,
            window_end + 1
        ):

            motif_set.update(
                self.position_index.get(pos, [])
            )

        return list(motif_set)

    # --------------------------------------------------
    # MUTATION HELPERS
    # --------------------------------------------------

    def apply_mutation(self, mutation):

        seq = list(self.sequence)

        seq[
            mutation.start:
            mutation.end + 1
        ] = list(mutation.new)

        return "".join(seq)


# ============================================================
# EVALUATION ENGINE
# ============================================================

def evaluate_mutation(
    genome,
    mutation,
    window_radius=25
):

    region = genome.get_region(
        mutation.position
    )

    # ----------------------------------------
    # Regulatory regions are locked
    # ----------------------------------------

    if region == RegionType.REGULATORY:

        return {
            "valid": False,
            "reason": "Protected region"
        }

    # ----------------------------------------
    # CDS synonymous check
    # ----------------------------------------

    if region == RegionType.CDS:

        if not is_synonymous(
            genome,
            mutation
        ):
            return {
                "valid": False,
                "reason": "Not synonymous"
            }

    # ----------------------------------------
    # Build mutated sequence
    # ----------------------------------------

    mutated_sequence = genome.apply_mutation(
        mutation
    )

    mutated_topology = genome.build_topology(
        mutated_sequence
    )

    # ----------------------------------------
    # Build windows
    # ----------------------------------------

    original_window = (
        genome.topology_engine.get_interval(
            mutation.start - window_radius,
            mutation.end + window_radius
        )
    )

    mutated_window = (
        mutated_topology.get_interval(
            mutation.start - window_radius,
            mutation.end + window_radius
        )
    )

    original_rc = reverse_complement(
        original_window
    )

    mutated_rc = reverse_complement(
        mutated_window
    )

    # ----------------------------------------
    # Local motifs only
    # ----------------------------------------

    local_motifs = genome.get_local_motifs(
        mutation.start,
        mutation.end,
        radius=window_radius
    )

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
            ) is not None
        )

        after = (
            motif.regex.search(
                mutated_window
            ) is not None
            or
            motif.regex.search(
                mutated_rc
            ) is not None
        )

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
        "edits": 1
    }


# ============================================================
# PLACEHOLDERS
# ============================================================

def is_synonymous(genome, mutation):
    """
    TODO

    Determine:
        - gene
        - codon affected
        - translated amino acid

    Return True if synonymous.
    """
    return True