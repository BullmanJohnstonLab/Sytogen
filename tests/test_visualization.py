from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sytogen.scripts.sytogen_runner import _parse_protected_regions
from sytogen.scripts.visualization import _build_circular_figure, _extract_protected_regions, _hit_tracks


def test_long_annotated_origin_is_protected_and_uses_its_note_as_label():
    record = SeqRecord(Seq("A" * 1000), id="example")
    record.features = [
        SeqFeature(
            FeatureLocation(0, 600),
            type="misc_feature",
            qualifiers={"note": ["p15A origin of replication region"]},
        )
    ]

    protected = _parse_protected_regions(record)
    display_regions = _extract_protected_regions(record)

    assert [(region.start, region.end, region.label) for region in protected] == [
        (0, 599, "p15A origin of replication region")
    ]
    assert display_regions == [{
        "id": "p15A origin of replication region",
        "start": 0,
        "end": 599,
        "strand": "+",
    }]


def test_motiffinder_tracks_are_grouped_by_enzyme_type():
    tracks, _ = _hit_tracks([
        {"rec_seq": "GATC", "pos_0": 10, "strand": "+", "enz_type": "2"},
        {"rec_seq": "CCWGG", "pos_0": 20, "strand": "-", "enz_type": "2"},
        {"rec_seq": "AACNNNNNNGTGC", "pos_0": 30, "strand": "+", "enz_type": "1"},
    ])

    assert [track["label"] for track in tracks] == ["Type I", "Type II"]
    assert len(tracks[1]["points"]) == 2


def test_circular_figure_allocates_space_for_horizontal_legend():
    fig = _build_circular_figure([], [], [], [], 1000, "Test map")

    assert fig.layout.height == 700
    assert fig.layout.margin.b == 120
    assert fig.layout.legend.y == -0.18
    assert fig.layout.hovermode == "closest"
    assert fig.layout.hoverdistance == 30
