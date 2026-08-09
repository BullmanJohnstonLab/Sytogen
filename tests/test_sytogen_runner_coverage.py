"""
Direct coverage for four sytogen_runner.py functions that were previously
only exercised indirectly, through the full-pipeline integration tests in
test_sytogen_page.py:

  - strip_backbone(): removing a known vector backbone from a construct
  - _parse_genes(): reconstructing an origin-spanning CompoundLocation gene
  - _final_new_motif_check(): the whole-construct post-edit safety net
  - _parse_codon_usage(): inverting rank-based ('lower is better') columns

Each fixture below was verified interactively against the real functions
before being written as a test.
"""

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sytogen.scripts.genome_model import Motif
from sytogen.scripts.sytogen_runner import (
    _final_new_motif_check,
    _parse_codon_usage,
    _parse_genes,
    strip_backbone,
)


# =========================================================
# strip_backbone
# =========================================================

def test_strip_backbone_handles_circular_wraparound():
    # The backbone straddles the origin: the plasmid record is rotated so
    # the last 4 bases + first 4 bases of the backbone are split across
    # the end/start of the sequence, with the insert in between.
    backbone_seq = "AAAATTTT"
    insert_seq = "GGGGCCCC"
    construct_seq = backbone_seq[4:] + insert_seq + backbone_seq[:4]

    construct_record = SeqRecord(Seq(construct_seq), id="construct")
    backbone_record = SeqRecord(Seq(backbone_seq), id="backbone")

    result = strip_backbone(construct_record, backbone_record, topology="circular")

    assert str(result.seq) == insert_seq


def test_strip_backbone_handles_simple_linear_case():
    construct_record = SeqRecord(Seq("AAAATTTT" + "GGGGCCCC"), id="construct")
    backbone_record = SeqRecord(Seq("AAAATTTT"), id="backbone")

    result = strip_backbone(construct_record, backbone_record, topology="linear")

    assert str(result.seq) == "GGGGCCCC"


def test_strip_backbone_raises_when_backbone_not_found():
    construct_record = SeqRecord(Seq("GGGGCCCC"), id="construct")
    backbone_record = SeqRecord(Seq("AAAATTTT"), id="backbone")

    with pytest.raises(ValueError):
        strip_backbone(construct_record, backbone_record, topology="linear")


def test_strip_backbone_raises_when_backbone_found_more_than_once():
    construct_record = SeqRecord(Seq("AAAATTTT" + "GGGG" + "AAAATTTT"), id="construct")
    backbone_record = SeqRecord(Seq("AAAATTTT"), id="backbone")

    with pytest.raises(ValueError):
        strip_backbone(construct_record, backbone_record, topology="linear")


# =========================================================
# _parse_genes: origin-spanning CompoundLocation
# =========================================================

def test_parse_genes_reconstructs_origin_spanning_compound_location():
    # 20bp circular plasmid; a CDS wraps the origin as GenBank would write
    # join(15..20,1..4) - the segment right before the origin, then the
    # segment right after it, 10bp total.
    sequence = "A" * 20
    record = SeqRecord(Seq(sequence), id="test")

    part1 = FeatureLocation(14, 20, strand=1)  # 1-based 15-20
    part2 = FeatureLocation(0, 4, strand=1)    # 1-based 1-4
    compound_loc = CompoundLocation([part1, part2])

    record.features.append(
        SeqFeature(compound_loc, type="CDS", qualifiers={"gene": ["wrapGene"]})
    )

    genes = _parse_genes(record)

    assert len(genes) == 1
    gene = genes[0]
    assert gene.id == "wrapGene"
    # Raw, un-clamped coordinates signaling a wrap: start=14, end=23
    # (23 = sequence_length(20) + part2.end(4) - 1).
    assert gene.start == 14
    assert gene.end == 23
    assert gene.strand == "+"


def test_parse_genes_leaves_a_normal_single_part_gene_alone():
    sequence = "A" * 20
    record = SeqRecord(Seq(sequence), id="test")
    record.features.append(
        SeqFeature(
            FeatureLocation(2, 11, strand=1),
            type="CDS",
            qualifiers={"gene": ["normalGene"]},
        )
    )

    genes = _parse_genes(record)

    assert len(genes) == 1
    gene = genes[0]
    assert gene.start == 2
    assert gene.end == 10  # end is inclusive, .end - 1 of the exclusive FeatureLocation


# =========================================================
# _final_new_motif_check
# =========================================================

def test_final_new_motif_check_detects_a_newly_introduced_site():
    original = "AAAA" + "CCCCCCCC" + "TTTT"
    final    = "AAAA" + "GAATTCCC" + "TTTT"  # GAATTC introduced at position 4
    motifs = [Motif("GAATTC", start=0, end=5, strand="+")]

    new_hits = _final_new_motif_check(original, final, motifs, topology="linear")

    assert len(new_hits) == 1
    assert new_hits[0]["motif"] == "GAATTC"
    assert new_hits[0]["start"] == 4


def test_final_new_motif_check_ignores_a_preexisting_site():
    original = "AAAA" + "GAATTCCC" + "TTTT"
    final    = "AAAA" + "GAATTCCC" + "TTTT"  # unchanged
    motifs = [Motif("GAATTC", start=4, end=9, strand="+")]

    new_hits = _final_new_motif_check(original, final, motifs, topology="linear")

    assert new_hits == []


def test_final_new_motif_check_short_circuits_on_empty_motif_list():
    assert _final_new_motif_check("AAAA", "TTTT", [], topology="linear") == []


# =========================================================
# _parse_codon_usage: rank-based column inversion
# =========================================================

def test_parse_codon_usage_inverts_ranking_column():
    # 'ranking': rank 1 = most preferred. After inversion, a HIGHER
    # usage_score must still mean a MORE preferred codon, same convention
    # as 'fraction'/'frequency'.
    df = pd.DataFrame({"codon": ["TTT", "TTC"], "ranking": [1, 2]})

    usage = _parse_codon_usage(df)

    assert usage["TTT"] > usage["TTC"]


def test_parse_codon_usage_inverts_ranking_ratio_column():
    # Same convention for 'ranking_ratio': lower raw value = more preferred.
    df = pd.DataFrame({"codon": ["TTT", "TTC"], "ranking_ratio": [0.9, 0.1]})

    usage = _parse_codon_usage(df)

    # TTC has the lower raw ranking_ratio (0.1), so it's more preferred
    # and must end up with the HIGHER usage_score after inversion.
    assert usage["TTC"] > usage["TTT"]


def test_parse_codon_usage_does_not_invert_fraction_column():
    # 'fraction' is already higher-is-better - no inversion should happen.
    df = pd.DataFrame({"codon": ["TTT", "TTC"], "fraction": [0.9, 0.1]})

    usage = _parse_codon_usage(df)

    assert usage["TTT"] > usage["TTC"]
    assert usage["TTT"] == 0.9
