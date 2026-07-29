from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from sytogen.scripts.sytogen_runner import _parse_protected_regions
from sytogen.scripts.genome_model import Motif
from sytogen.scripts.visualization import (
    _build_circular_figure,
    _extract_protected_regions,
    _hit_tracks,
    build_plasmid_maps,
)


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


def test_misc_feature_protected_label_prefers_note_text():
    record = SeqRecord(Seq("A" * 500), id="example")
    record.features = [
        SeqFeature(
            FeatureLocation(10, 40),
            type="misc_feature",
            qualifiers={
                "gene": ["gene_name_should_not_win"],
                "note": ["Quoted misc_feature note"],
            },
        )
    ]

    display_regions = _extract_protected_regions(record)

    assert display_regions == [{
        "id": "Quoted misc_feature note",
        "start": 10,
        "end": 39,
        "strand": "+",
    }]


def test_circular_figure_allocates_space_for_horizontal_legend():
    fig = _build_circular_figure([], [], [], [], 1000, "Test map")

    assert fig.layout.height == 760
    assert fig.layout.autosize is True
    assert fig.layout.margin.b == 120
    assert fig.layout.legend.y == -0.18
    assert fig.layout.hovermode == "closest"
    assert fig.layout.hoverdistance == 30


def test_circular_track_contains_toggle_metadata_for_filters():
    tracks = [{
        "label": "GATC",
        "color": "#123456",
        "ring": "Type II",
        "points": [{
            "position": 100,
            "hover": "demo",
            "status": "unresolved",
            "motif": "GATC",
            "type": "Type II",
        }],
    }]

    fig = _build_circular_figure([], [], [], tracks, 1000, "Test map")
    motif_trace = next(t for t in fig.data if getattr(t, "name", "") == "GATC")

    assert motif_trace.meta["layer"] == "motif_markers"
    assert motif_trace.meta["ring"] == "Type II"
    assert motif_trace.customdata[0]["status"] == "unresolved"
    assert motif_trace.customdata[0]["motif"] == "GATC"


def test_unresolved_motif_hover_explains_why_it_stayed_unresolved():
    record = SeqRecord(Seq("A" * 100), id="example")
    motif = Motif("GATC", 10, 13, "+", enz_type="2")
    decision_matrix = [{
        "motif": motif.motif,
        "motif_start": motif.start,
        "motif_end": motif.end,
        "motif_strand": motif.strand,
        "skip_reason": "all_candidates_rejected",
        "reasoning": "Every non-coding substitution attempted at this position was rejected (see attempted/rejected columns for the tally).",
        "attempted_count": 6,
        "rejected_count": 6,
        "top_rejection_reason": "does not destroy this motif",
        "top_rejection_count": 4,
        "chosen": False,
    }]

    _, fig_after = build_plasmid_maps(
        record,
        [motif],
        [],
        decision_matrix,
        resolved_motif_keys=set(),
        sequence_length=100,
        topology="linear",
        mask_regions=None,
        title="Test map",
    )

    trace = next(t for t in fig_after.data if getattr(t, "name", "") == "GATC")
    hover_text = trace.hovertext[0]

    assert "Why unresolved:" in hover_text
    assert "does not destroy this motif" in hover_text
    assert "Attempted candidates: 6" in hover_text


def test_unresolved_motif_hover_explains_when_selected_edit_did_not_fully_silence_it():
    record = SeqRecord(Seq("A" * 1000), id="example")
    motif = Motif("SCNGS", 606, 610, "+", enz_type="1")
    decision_matrix = [{
        "motif": motif.motif,
        "motif_start": motif.start,
        "motif_end": motif.end,
        "motif_strand": motif.strand,
        "edit_position": 608,
        "gene_id": "CDS_6",
        "gene_strand": "-",
        "before": "C",
        "after": "T",
        "original_codon": "GCA",
        "replacement_codon": "GCT",
        "AA_LetterCode": "A",
        "synonymous": True,
        "motifs_destroyed": 0,
        "reasoning": "Chosen: C→T at position 609, codon-usage score 0.600 (highest-scoring valid option for this motif).",
        "motifs_created": 0,
        "usage_score": 0.6,
        "gc_preserving": True,
        "total_score": 50.0,
        "chosen": True,
        "skip_reason": "",
        "attempted_count": "",
        "rejected_count": "",
        "top_rejection_reason": "",
        "top_rejection_count": "",
    }]

    _, fig_after = build_plasmid_maps(
        record,
        [motif],
        [],
        decision_matrix,
        resolved_motif_keys=set(),
        sequence_length=1000,
        topology="linear",
        mask_regions=None,
        title="Test map",
    )

    trace = next(t for t in fig_after.data if getattr(t, "name", "") == "SCNGS")
    hover_text = trace.hovertext[0]

    assert "best available edit was applied" in hover_text
    assert "Selected edit: C→T" in hover_text
    assert "gene CDS_6" in hover_text
    assert "position 609" in hover_text


def test_build_plasmid_maps_shows_deprotected_override_regions():
    record = SeqRecord(Seq("A" * 500), id="example")
    motif = Motif("GATC", 100, 103, "+", enz_type="2")

    fig_before, _ = build_plasmid_maps(
        record,
        [motif],
        [],
        [],
        resolved_motif_keys=set(),
        sequence_length=500,
        topology="linear",
        mask_regions=None,
        protected_override_ranges=[(120, 150)],
        title="Test map",
    )

    trace_names = [getattr(trace, "name", "") for trace in fig_before.data]
    assert "Deprotected (override)" in trace_names
