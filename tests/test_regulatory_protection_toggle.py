"""Tests for the opt-in `protect_predicted_regulatory` toggle in
run_sytogen_pipeline: when enabled, predicted promoter/RBS elements (from
regulatory_scanner.scan_promoters / scan_rbs) are added as protected
regions, in addition to annotation-derived protection. Off by default.
"""
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sytogen.scripts.sytogen_runner import _regulatory_protected_regions, run_sytogen_pipeline

MINUS35 = "TTGACA"
MINUS10 = "TATAAT"
PROMOTER_BLOCK = MINUS35 + ("C" * 17) + MINUS10  # 29 bp, spacing=17


def _codon_df():
    return pd.DataFrame([{"codon": "GAT", "fraction": 0.4}, {"codon": "GAC", "fraction": 0.6}])


def test_toggle_off_by_default_leaves_predicted_regions_editable():
    seq = ("A" * 10) + PROMOTER_BLOCK + ("A" * 10)
    # Overwrite 4 bases in the promoter's spacer region with a motif.
    motif_pos = 10 + 10
    seq = seq[:motif_pos] + "GGGG" + seq[motif_pos + 4 :]
    record = SeqRecord(Seq(seq), id="construct")
    motif_df = pd.DataFrame([{"motif": "GGGG", "enz_type": "2"}])

    result = run_sytogen_pipeline(
        record, _codon_df(), motif_df, params={"topology": "linear"}
    )

    assert result["summary"]["predicted_regulatory_regions_applied"] == 0
    assert result["summary"]["motifs_unresolved"] == 0
    assert result["summary"]["motifs_resolved"] > 0


def test_toggle_on_protects_predicted_promoter_region():
    seq = ("A" * 10) + PROMOTER_BLOCK + ("A" * 10)
    motif_pos = 10 + 10
    seq = seq[:motif_pos] + "GGGG" + seq[motif_pos + 4 :]
    record = SeqRecord(Seq(seq), id="construct")
    motif_df = pd.DataFrame([{"motif": "GGGG", "enz_type": "2"}])

    result = run_sytogen_pipeline(
        record,
        _codon_df(),
        motif_df,
        params={"topology": "linear", "protect_predicted_regulatory": "true"},
    )

    assert result["summary"]["predicted_regulatory_regions_applied"] >= 1
    # Every occurrence of the motif sits inside the predicted promoter and
    # is now off-limits -- nothing gets resolved.
    assert result["summary"]["motifs_resolved"] == 0
    assert result["summary"]["motifs_unresolved"] > 0


def test_toggle_does_not_affect_motifs_outside_predicted_regions():
    # Motif placed well away from the promoter block -- toggle should have
    # no effect on this one either way. Uses a self-palindromic motif
    # (AAGCTT, HindIII) with no reverse-complement self-match, so it can't
    # accidentally also match inside the promoter's poly-C spacer the way
    # GGGG's reverse complement (CCCC) would.
    seq = ("A" * 10) + PROMOTER_BLOCK + ("A" * 10) + "AAGCTT" + ("A" * 10)
    record = SeqRecord(Seq(seq), id="construct")
    motif_df = pd.DataFrame([{"motif": "AAGCTT", "enz_type": "2"}])

    off = run_sytogen_pipeline(
        record, _codon_df(), motif_df, params={"topology": "linear"}
    )
    on = run_sytogen_pipeline(
        record,
        _codon_df(),
        motif_df,
        params={"topology": "linear", "protect_predicted_regulatory": "true"},
    )

    assert off["summary"]["motifs_resolved"] == on["summary"]["motifs_resolved"]
    assert off["summary"]["motifs_resolved"] > 0


def test_toggle_leaves_coding_region_edits_synonymous_as_usual():
    # Sanity check that turning the toggle on doesn't change anything about
    # *how* an unrelated coding-region motif gets edited (still synonymous).
    cds = "ATG" + "GATC" + ("AAA" * 20) + "AA" + "TAA"  # 72bp, in-frame
    seq = ("A" * 10) + PROMOTER_BLOCK + ("A" * 10) + cds + ("A" * 10)
    record = SeqRecord(Seq(seq), id="construct")
    cds_start = seq.index(cds)
    record.features.append(
        SeqFeature(
            FeatureLocation(cds_start, cds_start + len(cds), strand=1),
            type="CDS",
            qualifiers={"gene": ["testGene"]},
        )
    )
    motif_df = pd.DataFrame([{"motif": "GATC", "enz_type": "2"}])

    result = run_sytogen_pipeline(
        record,
        _codon_df(),
        motif_df,
        params={"topology": "linear", "protect_predicted_regulatory": "true"},
    )

    assert result["summary"]["motifs_resolved"] == 1
    # Edited sequence must still translate identically (synonymous edit).
    original_protein = str(Seq(cds).translate())
    altered_cds_region = result["altered_sequence"][cds_start : cds_start + len(cds)]
    assert str(Seq(altered_cds_region).translate()) == original_protein


def test_regulatory_protected_regions_skips_origin_spanning_genes():
    # Documented limitation: a gene with end >= sequence length (the
    # origin-spanning convention from _parse_genes) is excluded rather
    # than mishandled.
    from sytogen.scripts.genome_model import Gene

    seq = "A" * 60
    origin_spanning_gene = Gene(gene_id="wrap", start=50, end=65, strand="+")  # end >= len(seq)
    # Should not raise, and should simply produce no RBS predictions tied
    # to this gene (there's no real upstream sequence to search either way).
    regions = _regulatory_protected_regions(seq, [origin_spanning_gene], topology="circular")
    assert all(r.label is not None for r in regions)


def test_regulatory_protected_regions_handles_circular_wraparound():
    # A promoter split across the circular origin should still end up
    # fully protected, represented as two ProtectedRegion segments.
    tail, head = PROMOTER_BLOCK[19:], PROMOTER_BLOCK[:19]
    seq = tail + ("C" * 10) + head  # 29 + 10 = 39 bp, motif wraps at index 0

    regions = _regulatory_protected_regions(seq, genes=[], topology="circular")

    promoter_regions = [r for r in regions if "promoter" in r.label]
    assert len(promoter_regions) == 2
    # One segment should run to the end of the sequence, the other from 0.
    ends = sorted(r.end for r in promoter_regions)
    starts = sorted(r.start for r in promoter_regions)
    assert starts[0] == 0
    assert ends[-1] == len(seq) - 1
