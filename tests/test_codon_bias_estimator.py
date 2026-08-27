"""Tests for sytogen.scripts.codon_bias_estimator.

This module previously had zero test coverage. Writing these tests
surfaced a real bug (fixed alongside this suite): parse_fasta_gff()
parsed the GFF strand column but never used it -- a '-'-strand CDS was
extracted as a raw forward-strand slice instead of being reverse-
complemented, so codon usage for minus-strand genes in FASTA+GFF mode
was computed from antisense sequence. GenBank-mode input (parse_genbank,
via Bio's feat.extract) was already strand-correct and is unaffected.
See test_parse_fasta_gff_reverse_complements_minus_strand_cds.
"""
import os
import tempfile

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sytogen.scripts.codon_bias_estimator import (
    compute_codon_usage,
    count_codons,
    load_fasta_records,
    parse_fasta_gff,
    parse_genbank,
    run_codon_bias,
    standard_codon_usage,
)


# ---------------------------------------------------------------------------
# count_codons
# ---------------------------------------------------------------------------

def test_count_codons_counts_non_overlapping_triplets():
    assert count_codons("ATGAAATAA") == {"ATG": 1, "AAA": 1, "TAA": 1}


def test_count_codons_ignores_trailing_partial_codon():
    # 10 bases -> 3 full codons, trailing "A" dropped.
    assert count_codons("ATGAAATAAA") == {"ATG": 1, "AAA": 1, "TAA": 1}


def test_count_codons_counts_repeats():
    assert count_codons("ATGATGATG") == {"ATG": 3}


# ---------------------------------------------------------------------------
# compute_codon_usage
# ---------------------------------------------------------------------------

def test_compute_codon_usage_proportions_and_ranking():
    # Met x2 (single-codon AA), Ile x1 of 3 possible codons, Lys x1 of 2,
    # Stop x1 of 3.
    seq = "ATGATGATCAAATAA"
    df = compute_codon_usage(seq, codon_table=11).set_index("Codon")

    # Single-codon amino acid: proportion/ranking are trivially 1.
    assert df.loc["ATG", "AA"] == "M"
    assert df.loc["ATG", "Count"] == 2
    assert df.loc["ATG", "Proportion"] == 1.0
    assert df.loc["ATG", "Ranking_ratio"] == 1.0

    # Ile: only ATC observed, out of {ATT, ATC, ATA}.
    assert df.loc["ATC", "Count"] == 1
    assert df.loc["ATC", "Proportion"] == 1.0
    assert df.loc["ATC", "Ranking"] == 1.0
    assert df.loc["ATT", "Count"] == 0
    assert df.loc["ATT", "Proportion"] == 0.0
    # Unused synonymous codons tie for rank 2 (dense ranking).
    assert df.loc["ATT", "Ranking"] == 2.0
    assert df.loc["ATA", "Ranking"] == 2.0

    # Lys: AAA observed once out of {AAA, AAG}.
    assert df.loc["AAA", "Proportion"] == 1.0
    assert df.loc["AAG", "Proportion"] == 0.0


def test_compute_codon_usage_unused_amino_acid_is_zero_not_nan():
    # A sequence using only ATG (Met) and TAA (Stop) -- every other amino
    # acid's codons must be present in the table with Proportion 0, not
    # NaN or a divide-by-zero error.
    df = compute_codon_usage("ATGTAA", codon_table=11)
    leu_rows = df[df["AA"] == "L"]
    assert len(leu_rows) > 0
    assert (leu_rows["Proportion"] == 0).all()
    assert not leu_rows["Proportion"].isna().any()


def test_compute_codon_usage_includes_every_codon_in_the_table_even_unobserved():
    df = compute_codon_usage("ATG", codon_table=11)
    # Table 11 has 64 codons total (61 sense + 3 stop).
    assert len(df) == 64
    assert df["Count"].sum() == 1


# ---------------------------------------------------------------------------
# standard_codon_usage
# ---------------------------------------------------------------------------

def test_standard_codon_usage_gives_equal_weight_within_each_amino_acid():
    df = standard_codon_usage(codon_table=11)
    # Every amino acid's fractions must sum to (approximately) 1.
    from Bio.Data import CodonTable

    codon_table_obj = CodonTable.unambiguous_dna_by_id[11]
    codon_to_aa = dict(codon_table_obj.forward_table)
    for stop in codon_table_obj.stop_codons:
        codon_to_aa[stop] = "Stop"

    df = df.assign(aa=df["codon"].map(codon_to_aa))
    totals = df.groupby("aa")["fraction"].sum()
    assert totals.apply(lambda x: abs(x - 1.0) < 1e-9).all()


def test_standard_codon_usage_single_codon_amino_acids_get_full_weight():
    df = standard_codon_usage(codon_table=11).set_index("codon")
    assert df.loc["ATG", "fraction"] == 1.0  # Met: only codon
    assert df.loc["TGG", "fraction"] == 1.0  # Trp: only codon


# ---------------------------------------------------------------------------
# parse_genbank
# ---------------------------------------------------------------------------

def _write_genbank(tmpdir, seq, features, record_id="chr1"):
    rec = SeqRecord(Seq(seq), id=record_id)
    rec.annotations["molecule_type"] = "DNA"
    for ftype, start, end, strand, qualifiers in features:
        rec.features.append(
            SeqFeature(FeatureLocation(start, end, strand=strand), type=ftype, qualifiers=qualifiers or {})
        )
    path = os.path.join(tmpdir, "g.gbk")
    SeqIO.write([rec], path, "genbank")
    return path


def test_parse_genbank_extracts_plus_strand_cds():
    with tempfile.TemporaryDirectory() as d:
        seq = "CCCCCC" + "ATGAAATAA" + "CCCCCC"
        path = _write_genbank(d, seq, [("CDS", 6, 15, 1, None)])
        assert parse_genbank(path) == "ATGAAATAA"


def test_parse_genbank_reverse_complements_minus_strand_cds():
    with tempfile.TemporaryDirectory() as d:
        orf = "ATGAAATAA"
        rc_orf = str(Seq(orf).reverse_complement())
        seq = "CCCCCC" + rc_orf + "CCCCCC"
        path = _write_genbank(d, seq, [("CDS", 6, 15, -1, None)])
        assert parse_genbank(path) == orf


def test_parse_genbank_applies_codon_start_offset():
    with tempfile.TemporaryDirectory() as d:
        # codon_start=2 means the first complete codon starts 1 base in.
        seq = "CCCCCC" + "X" + "ATGAAATAA" + "CCCCCC"
        seq = seq.replace("X", "A")  # leading junk base before the real ORF
        path = _write_genbank(
            d, seq, [("CDS", 6, 16, 1, {"codon_start": ["2"]})]
        )
        assert parse_genbank(path) == "ATGAAATAA"


def test_parse_genbank_accepts_orf_and_gene_feature_types():
    with tempfile.TemporaryDirectory() as d:
        seq = "CCCCCC" + "ATGAAATAA" + "CCCCCC"
        for ftype in ["ORF", "gene"]:
            path = _write_genbank(d, seq, [(ftype, 6, 15, 1, None)], record_id=ftype)
            assert parse_genbank(path) == "ATGAAATAA"


def test_parse_genbank_skips_cds_shorter_than_9bp():
    with tempfile.TemporaryDirectory() as d:
        seq = "CCCCCC" + "ATGTAA" + "CCCCCC"  # 6bp CDS, below the 9bp floor
        path = _write_genbank(d, seq, [("CDS", 6, 12, 1, None)])
        with pytest.raises(ValueError, match="No valid CDS/ORF"):
            parse_genbank(path)


def test_parse_genbank_skips_cds_with_ambiguous_bases():
    with tempfile.TemporaryDirectory() as d:
        seq = "CCCCCC" + "ATGNNNTAA" + "CCCCCC"
        path = _write_genbank(d, seq, [("CDS", 6, 15, 1, None)])
        with pytest.raises(ValueError, match="No valid CDS/ORF"):
            parse_genbank(path)


def test_parse_genbank_concatenates_multiple_cds():
    with tempfile.TemporaryDirectory() as d:
        cds1 = "ATGAAATAA"
        cds2 = "ATGCCCTAA"
        seq = cds1 + "CCCCCC" + cds2
        path = _write_genbank(
            d, seq,
            [("CDS", 0, 9, 1, None), ("CDS", 15, 24, 1, None)],
        )
        assert parse_genbank(path) == cds1 + cds2


def test_parse_genbank_raises_with_no_valid_features():
    with tempfile.TemporaryDirectory() as d:
        path = _write_genbank(d, "CCCCCCCCCCCCCCCC", [("misc_feature", 0, 10, 1, None)])
        with pytest.raises(ValueError, match="No valid CDS/ORF"):
            parse_genbank(path)


# ---------------------------------------------------------------------------
# parse_fasta_gff
# ---------------------------------------------------------------------------

def _write_fasta_gff(tmpdir, seq, gff_lines, record_id="chr1"):
    fasta_path = os.path.join(tmpdir, "g.fasta")
    gff_path = os.path.join(tmpdir, "g.gff3")
    with open(fasta_path, "w") as f:
        f.write(f">{record_id}\n{seq}\n")
    with open(gff_path, "w") as f:
        f.write("##gff-version 3\n")
        for line in gff_lines:
            f.write(line + "\n")
    return fasta_path, gff_path


def test_parse_fasta_gff_extracts_plus_strand_cds():
    with tempfile.TemporaryDirectory() as d:
        orf = "ATGAAATAA"
        seq = "CCCCCC" + orf + "CCCCCC"
        fasta_path, gff_path = _write_fasta_gff(
            d, seq, ["chr1\tsrc\tCDS\t7\t15\t.\t+\t0\tID=g1"]
        )
        assert parse_fasta_gff(fasta_path, gff_path) == orf


def test_parse_fasta_gff_reverse_complements_minus_strand_cds():
    # Regression test for a fixed bug: the GFF strand column was parsed
    # but silently ignored, so a '-'-strand CDS returned the raw forward
    # genomic slice (antisense to the real coding sequence) instead of
    # its reverse complement.
    with tempfile.TemporaryDirectory() as d:
        orf = "ATGAAATAA"
        rc_orf = str(Seq(orf).reverse_complement())
        seq = "CCCCCC" + rc_orf + "CCCCCC"
        fasta_path, gff_path = _write_fasta_gff(
            d, seq, ["chr1\tsrc\tCDS\t7\t15\t.\t-\t0\tID=g1"]
        )
        assert parse_fasta_gff(fasta_path, gff_path) == orf


def test_parse_fasta_gff_skips_short_and_ambiguous_features():
    with tempfile.TemporaryDirectory() as d:
        seq = "ATGTAA" + "NNNNNNNNN"  # 6bp (too short) + 9bp with N's
        fasta_path, gff_path = _write_fasta_gff(
            d, seq,
            [
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=short",
                "chr1\tsrc\tCDS\t7\t15\t.\t+\t0\tID=ambiguous",
            ],
        )
        with pytest.raises(ValueError, match="No valid CDS/ORF"):
            parse_fasta_gff(fasta_path, gff_path)


def test_parse_fasta_gff_ignores_features_on_unknown_seqid():
    with tempfile.TemporaryDirectory() as d:
        orf = "ATGAAATAA"
        seq = "CCCCCC" + orf + "CCCCCC"
        fasta_path, gff_path = _write_fasta_gff(
            d, seq, ["not_a_real_chrom\tsrc\tCDS\t7\t15\t.\t+\t0\tID=g1"]
        )
        with pytest.raises(ValueError, match="No valid CDS/ORF"):
            parse_fasta_gff(fasta_path, gff_path)


# ---------------------------------------------------------------------------
# load_fasta_records
# ---------------------------------------------------------------------------

def test_load_fasta_records_keys_by_short_id_and_full_header():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "g.fasta")
        with open(path, "w") as f:
            f.write(">chr1 some description here\nACGTACGT\n")
        records = load_fasta_records(path)
        assert "chr1" in records
        assert "chr1 some description here" in records
        assert str(records["chr1"].seq) == "ACGTACGT"


def test_load_fasta_records_raises_on_empty_file():
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "empty.fasta")
        open(path, "w").close()
        with pytest.raises(ValueError, match="No sequences found"):
            load_fasta_records(path)


# ---------------------------------------------------------------------------
# run_codon_bias (end-to-end)
# ---------------------------------------------------------------------------

def test_run_codon_bias_genbank_mode_produces_consistent_csv():
    with tempfile.TemporaryDirectory() as d:
        seq = "CCCCCC" + "ATGAAATAA" + "CCCCCC"
        genome_path = _write_genbank(d, seq, [("CDS", 6, 15, 1, None)])
        out_dir = os.path.join(d, "out")
        result = run_codon_bias(genome_path=genome_path, output_dir=out_dir)

        assert os.path.exists(result["csv"])
        df = pd.read_csv(result["csv"])
        expected = compute_codon_usage("ATGAAATAA", codon_table=11)
        pd.testing.assert_frame_equal(
            df.sort_values("Codon").reset_index(drop=True),
            expected.sort_values("Codon").reset_index(drop=True),
        )


def test_run_codon_bias_fasta_gff_mode_respects_minus_strand():
    with tempfile.TemporaryDirectory() as d:
        orf = "ATGAAATAA"
        rc_orf = str(Seq(orf).reverse_complement())
        seq = "CCCCCC" + rc_orf + "CCCCCC"
        fasta_path, gff_path = _write_fasta_gff(
            d, seq, ["chr1\tsrc\tCDS\t7\t15\t.\t-\t0\tID=g1"]
        )
        out_dir = os.path.join(d, "out")
        result = run_codon_bias(fasta_path=fasta_path, gff_path=gff_path, output_dir=out_dir)

        df = pd.read_csv(result["csv"])
        expected = compute_codon_usage(orf, codon_table=11)
        pd.testing.assert_frame_equal(
            df.sort_values("Codon").reset_index(drop=True),
            expected.sort_values("Codon").reset_index(drop=True),
        )


def test_run_codon_bias_requires_genome_or_fasta_gff_pair():
    with pytest.raises(ValueError, match="Provide either"):
        run_codon_bias(output_dir=tempfile.mkdtemp())
