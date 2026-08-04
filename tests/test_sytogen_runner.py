from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from types import SimpleNamespace

from sytogen.scripts.sytogen_runner import (
    _apply_protected_region_overrides,
    _candidate_priority,
    _make_matrix_row,
    _parse_protected_override_ranges,
    _parse_protected_regions,
    decision_matrix_to_tsv,
)


def test_parse_protected_override_ranges_parses_1_based_inclusive_windows():
    ranges = _parse_protected_override_ranges("10-20, 100-105", sequence_length=500)
    assert ranges == [(9, 19), (99, 104)]


def test_protected_override_ranges_remove_overlapping_annotation_regions():
    record = SeqRecord(Seq("A" * 600), id="example")
    record.features = [
        SeqFeature(
            FeatureLocation(9, 20),
            type="misc_feature",
            qualifiers={"note": ["promoter_A"]},
        ),
        SeqFeature(
            FeatureLocation(89, 110),
            type="misc_feature",
            qualifiers={"note": ["promoter_B"]},
        ),
        SeqFeature(
            FeatureLocation(199, 250),
            type="rep_origin",
            qualifiers={"note": ["origin_region"]},
        ),
    ]

    protected = _parse_protected_regions(record)
    overrides = _parse_protected_override_ranges("15-95", sequence_length=len(record.seq))
    filtered = _apply_protected_region_overrides(protected, overrides)

    assert [r.label for r in filtered] == ["origin_region"]


def test_candidate_priority_prefers_motif_destruction_over_usage_score():
    destroying_candidate = SimpleNamespace(
        result={"destroyed": 1, "edits": 1},
        usage_score=0.05,
    )
    non_destroying_candidate = SimpleNamespace(
        result={"destroyed": 0, "edits": 0},
        usage_score=0.99,
    )

    assert _candidate_priority(destroying_candidate, False, True) > _candidate_priority(non_destroying_candidate, True, True)


def test_candidate_priority_prefers_gc_preserving_when_enabled():
    gc_preserving_candidate = SimpleNamespace(
        result={"destroyed": 1, "edits": 1, "overlap_priority": 0},
        usage_score=0.5,
    )
    non_gc_preserving_candidate = SimpleNamespace(
        result={"destroyed": 1, "edits": 1, "overlap_priority": 0},
        usage_score=0.9,
    )

    assert _candidate_priority(gc_preserving_candidate, True, True) > _candidate_priority(non_gc_preserving_candidate, False, True)
    assert _candidate_priority(gc_preserving_candidate, True, False) < _candidate_priority(non_gc_preserving_candidate, False, False)


def test_decision_matrix_tsv_drops_columns_that_never_populate():
    matrix = [
        {
            "motif": "GATC",
            "motif_start": 100,
            "motif_end": 103,
            "motif_strand": "+",
            "edit_position": 101,
            "gene_id": "CDS_1",
            "gene_strand": "+",
            "before": "A",
            "after": "C",
            "original_codon": "GAT",
            "replacement_codon": "GCT",
            "AA_LetterCode": "D",
            "synonymous": True,
            "motifs_destroyed": 1,
            "reasoning": "Chosen candidate.",
            "motifs_created": 0,
            "usage_score": 0.2,
            "gc_preserving": False,
            "total_score": 1010.0,
            "chosen": True,
            "skip_reason": "",
            "attempted_count": "",
            "rejected_count": "",
            "top_rejection_reason": "",
            "top_rejection_count": "",
        }
    ]

    tsv = decision_matrix_to_tsv(matrix)
    header = tsv.splitlines()[0].split("\t")
    row = tsv.splitlines()[1].split("\t")
    index = {name: idx for idx, name in enumerate(header)}

    assert "before" in header
    assert "after" in header
    assert "attempted_count" not in header
    assert "rejected_count" not in header
    assert "top_rejection_reason" not in header
    assert "top_rejection_count" not in header
    assert row[index["motif_start"]] == "101"
    assert row[index["motif_end"]] == "104"
    assert row[index["edit_position"]] == "102"


def test_decision_matrix_tsv_keeps_columns_when_any_row_populates_them():
    matrix = [
        {
            "motif": "GATC",
            "motif_start": 100,
            "motif_end": 103,
            "motif_strand": "+",
            "edit_position": "",
            "gene_id": "",
            "gene_strand": "",
            "before": "",
            "after": "",
            "original_codon": "",
            "replacement_codon": "",
            "AA_LetterCode": "",
            "synonymous": "",
            "motifs_destroyed": 0,
            "reasoning": "No valid candidate could be constructed for this motif.",
            "motifs_created": 0,
            "usage_score": 0,
            "gc_preserving": "",
            "total_score": 0,
            "chosen": False,
            "skip_reason": "all_candidates_rejected",
            "attempted_count": 6,
            "rejected_count": 6,
            "top_rejection_reason": "does not destroy this motif",
            "top_rejection_count": 4,
        }
    ]

    tsv = decision_matrix_to_tsv(matrix)
    header = tsv.splitlines()[0].split("\t")

    assert "attempted_count" in header
    assert "rejected_count" in header
    assert "top_rejection_reason" in header
    assert "top_rejection_count" in header


def test_make_matrix_row_uses_genomic_bases_in_reasoning():
    genome = SimpleNamespace(
        sequence="AAACCCGGGTTT",
        find_gene=lambda position: SimpleNamespace(id="CDS_6", strand="-") if position == 5 else None,
    )
    motif = SimpleNamespace(motif="GGCC", start=4, end=7, strand="-")
    candidate = SimpleNamespace(
        mutation=SimpleNamespace(position=5, start=5, end=5, old="C", new="T"),
        result={"destroyed": 1, "created": 0, "edits": 1},
        codon="AGG",
        replacement="AGA",
        usage_score=0.6,
    )

    row = _make_matrix_row(motif, candidate, 101.0, True, genome)

    assert row["before"] == "C"
    assert row["after"] == "T"
    assert "C→T at position 6" in row["reasoning"]
    assert "AGG→AGA" not in row["reasoning"]
