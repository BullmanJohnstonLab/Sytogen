from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sytogen.scripts.sytogen_runner import (
    _apply_protected_region_overrides,
    _parse_protected_override_ranges,
    _parse_protected_regions,
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
