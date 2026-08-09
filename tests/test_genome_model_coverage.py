"""
Additional coverage for genome_model.py.

The existing test_genome_model.py covers motif-overlap candidate ranking.
This file targets the other delicate areas of the module - the ones
called out by "FIX:" comments in the source as places that have already
had real bugs (circular-origin wraparound, minus-strand codon indexing) -
plus the basic building blocks (region priority, mutation validation)
that the rest of the pipeline depends on.
"""

import Bio.Seq
import pytest

from sytogen.scripts.genome_model import (
    Gene,
    GenomeModel,
    LinearTopology,
    CircularTopology,
    Mutation,
    ProtectedRegion,
    RegionType,
)
from sytogen.scripts.sequence_utils import compile_iupac


# =========================================================
# Topology: circular-origin wraparound
# =========================================================

def test_circular_topology_get_interval_wraps_across_origin():
    sequence = "ACGTACGTACGTACGTACGT"  # 20 bases
    topo = CircularTopology(sequence)

    result = topo.get_interval(18, 22)

    assert result == sequence[18:] + sequence[:2]
    assert len(result) == 4


def test_circular_topology_normalize_position_wraps():
    topo = CircularTopology("ACGT" * 5)  # length 20

    assert topo.normalize_position(-1) == 19
    assert topo.normalize_position(20) == 0
    assert topo.normalize_position(25) == 5


def test_linear_topology_normalize_position_clamps():
    topo = LinearTopology("ACGT" * 5)  # length 20

    assert topo.normalize_position(-5) == 0
    assert topo.normalize_position(100) == 19


def test_circular_topology_count_motif_hits_finds_origin_spanning_site():
    # "GAATTC" (EcoRI) is split as "GA" at the very end and "ATTC" at the
    # very start - only detectable if hit-counting actually wraps.
    sequence = "ATTC" + "A" * 10 + "GA"
    topo = CircularTopology(sequence)
    regex = compile_iupac("GAATTC")

    assert topo.count_motif_hits(regex) == 1


# =========================================================
# Gene: circular-origin wraparound
# =========================================================

def test_gene_contains_handles_circular_origin_wrap():
    # Raw start/end of 18-25 on a 20bp genome represents a gene that
    # wraps: genomic positions 18, 19, 0, 1, 2, 3, 4, 5 (8 bases).
    gene = Gene("wrap", start=18, end=25, strand="+")
    genome_length = 20

    assert gene.contains(19, genome_length) is True   # just before origin
    assert gene.contains(1, genome_length) is True    # just after origin
    assert gene.contains(5, genome_length) is True    # last base (boundary)
    assert gene.contains(6, genome_length) is False   # one past the end
    assert gene.contains(10, genome_length) is False  # nowhere near the gene


# =========================================================
# Gene: minus-strand codon handling
# =========================================================
# This is the historically bug-prone area (see the "FIX:" comment on
# Gene.codon_mutations for the specific bug this used to have).

@pytest.fixture
def minus_strand_gene_fixture():
    """
    A 9bp plus-strand sequence encoding, on the MINUS strand read 5'->3',
    the protein Met-Phe-Stop (ATG-TTT-TAA). Verified independently against
    Gene.get_codon() before being used as a fixture.
    """
    plus_seq = "TTA" + "AAA" + "CAT"
    gene = Gene("g1", start=0, end=8, strand="-")
    genome = GenomeModel(
        sequence=plus_seq,
        topology="linear",
        genes=[gene],
        motifs=[],
        protected_regions=[],
        codon_usage={},
    )
    return genome, gene


def test_minus_strand_get_codon_reads_coding_strand_5_to_3(minus_strand_gene_fixture):
    genome, gene = minus_strand_gene_fixture

    # Gene position (nearest 'end') = first mRNA codon = Met
    assert gene.get_codon(genome, 8) == "ATG"
    # Middle codon = Phe
    assert gene.get_codon(genome, 5) == "TTT"
    # Last codon (nearest 'start') = Stop
    assert gene.get_codon(genome, 2) == "TAA"


def test_minus_strand_codon_mutations_diff_at_non_middle_position():
    """
    Regression test for the bug documented in Gene.codon_mutations: for a
    minus-strand gene, a coding-frame base difference that ISN'T the
    middle position must still map to the correct genomic base and
    position. TTT -> TTC (both Phe) differs at coding-frame index 2 (the
    last base), which is exactly the case the fix addresses.
    """
    gene = Gene("g1", start=0, end=8, strand="-")

    diffs = gene.codon_mutations(
        codon_start=3,
        original_codon="TTT",
        replacement_codon="TTC",
        genome_length=9,
    )

    assert len(diffs) == 1
    mutation = diffs[0]
    assert mutation.position == 3
    assert mutation.old == "A"
    assert mutation.new == "G"


def test_minus_strand_codon_mutations_round_trips_through_apply_mutation():
    """
    The genomic Mutation produced for a minus-strand codon swap must be
    directly usable by GenomeModel.apply_mutation() - i.e. mutation.old
    must match what's actually at that position in the plus-strand
    sequence, and applying it must produce the expected new codon.
    """
    plus_seq = "TTA" + "AAA" + "CAT"
    gene = Gene("g1", start=0, end=8, strand="-")
    genome = GenomeModel(
        sequence=plus_seq,
        topology="linear",
        genes=[gene],
        motifs=[],
        protected_regions=[],
        codon_usage={},
    )

    diffs = gene.codon_mutations(
        codon_start=3,
        original_codon="TTT",
        replacement_codon="TTC",
        genome_length=len(plus_seq),
    )
    assert len(diffs) == 1

    mutated_sequence = genome.apply_mutation(diffs[0])
    mutated_gene_codon = gene.get_codon(
        GenomeModel(
            sequence=mutated_sequence,
            topology="linear",
            genes=[gene],
            motifs=[],
            protected_regions=[],
            codon_usage={},
        ),
        5,
    )
    assert mutated_gene_codon == "TTC"
    assert str(Bio.Seq.Seq(mutated_gene_codon).translate()) == "F"


def test_synonymous_codons_excludes_the_input_codon():
    gene = Gene("g", start=0, end=8, strand="+")

    synonyms = gene.synonymous_codons("TTT")  # Phe

    assert "TTT" not in synonyms
    assert "TTC" in synonyms  # the other Phe codon


def test_ranked_synonymous_codons_orders_by_usage():
    gene = Gene("g", start=0, end=8, strand="+")
    codon_usage = {"TTC": 0.9, "TTT": 0.1}  # TTT excluded anyway (it's the input)

    ranked = gene.ranked_synonymous_codons("TTT", codon_usage)

    assert ranked[0] == "TTC"


# =========================================================
# GenomeModel.apply_mutation
# =========================================================

def test_apply_mutation_replaces_the_correct_base():
    genome = GenomeModel(
        sequence="AAAACCCCGGGGTTTT",
        topology="linear",
        genes=[],
        motifs=[],
        protected_regions=[],
        codon_usage={},
    )

    mutated = genome.apply_mutation(Mutation(position=4, old="C", new="T"))

    assert mutated == "AAAATCCCGGGGTTTT"
    # Original sequence is untouched (apply_mutation must not mutate in place).
    assert genome.sequence == "AAAACCCCGGGGTTTT"


def test_apply_mutation_raises_on_mismatched_expected_base():
    genome = GenomeModel(
        sequence="AAAACCCCGGGGTTTT",
        topology="linear",
        genes=[],
        motifs=[],
        protected_regions=[],
        codon_usage={},
    )

    with pytest.raises(ValueError):
        genome.apply_mutation(Mutation(position=4, old="G", new="T"))


# =========================================================
# GenomeModel.get_region priority: protected > CDS > neutral
# =========================================================

def test_get_region_priority_protected_over_cds_over_neutral():
    protected = ProtectedRegion(start=5, end=9, source="test")
    gene = Gene("g", start=5, end=13, strand="+")
    genome = GenomeModel(
        sequence="A" * 20,
        topology="linear",
        genes=[gene],
        motifs=[],
        protected_regions=[protected],
        codon_usage={},
    )

    # Overlap of protected region and CDS -> protected wins.
    assert genome.get_region(6) == RegionType.REGULATORY
    # CDS only.
    assert genome.get_region(12) == RegionType.CDS
    # Neither.
    assert genome.get_region(17) == RegionType.NEUTRAL


def test_is_protected_matches_get_region():
    protected = ProtectedRegion(start=5, end=9, source="test")
    genome = GenomeModel(
        sequence="A" * 20,
        topology="linear",
        genes=[],
        motifs=[],
        protected_regions=[protected],
        codon_usage={},
    )

    assert genome.is_protected(6) is True
    assert genome.is_protected(15) is False


# =========================================================
# GenomeModel.evaluate_mutation
# =========================================================

def test_evaluate_mutation_rejects_edit_inside_protected_region():
    protected = ProtectedRegion(start=5, end=9, source="test")
    genome = GenomeModel(
        sequence="A" * 20,
        topology="linear",
        genes=[],
        motifs=[],
        protected_regions=[protected],
        codon_usage={},
    )

    result = genome.evaluate_mutation(Mutation(position=6, old="A", new="T"))

    assert result["valid"] is False
    assert result["reason"] == "Protected region"


def test_evaluate_mutation_rejects_non_synonymous_cds_edit():
    plus_seq = "ATG" + "AAA" + "TAA"  # Met-Lys-Stop
    gene = Gene("g", start=0, end=8, strand="+")
    genome = GenomeModel(
        sequence=plus_seq,
        topology="linear",
        genes=[gene],
        motifs=[],
        protected_regions=[],
        codon_usage={},
    )

    # A -> T at the first base of ATG (Met) gives TTG (Leu): not synonymous.
    result = genome.evaluate_mutation(Mutation(position=0, old="A", new="T"))

    assert result["valid"] is False
    assert result["reason"] == "Not synonymous"


def test_evaluate_mutation_accepts_synonymous_cds_edit():
    plus_seq = "ATG" + "AAA" + "TAA"  # Met-Lys-Stop
    gene = Gene("g", start=0, end=8, strand="+")
    genome = GenomeModel(
        sequence=plus_seq,
        topology="linear",
        genes=[gene],
        motifs=[],
        protected_regions=[],
        codon_usage={},
    )

    # AAA -> AAG: still Lys.
    result = genome.evaluate_mutation(Mutation(position=5, old="A", new="G"))

    assert result["valid"] is True
    assert result["edits"] == 1


def test_evaluate_mutation_accepts_neutral_edit_outside_any_region():
    genome = GenomeModel(
        sequence="A" * 20,
        topology="linear",
        genes=[],
        motifs=[],
        protected_regions=[],
        codon_usage={},
    )

    result = genome.evaluate_mutation(Mutation(position=10, old="A", new="T"))

    assert result["valid"] is True
