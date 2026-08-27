"""Tests for sytogen.scripts.regulatory_scanner (scan_promoters / scan_rbs).

These are unit tests against hand-built sequences with known -35/-10 or RBS
placements, plus one integration test exercising the MotifFinder API endpoint
end-to-end with the `regulatory_scan` toggle enabled. Coordinates returned by
both functions are cross-checked against the sequence itself (reverse-
complementing for '-' strand hits) rather than hard-coded, so a test failure
means the returned coordinates don't actually point at the reported motif.

The three test_scan_rbs_minus_strand_* tests are regression tests for two
bugs (both in scan_rbs()'s handling of '-'-strand features) that were found
and fixed via this test suite: an incorrect coordinate transform for linear
topology, and a missing complement step (reverse without complement) for
circular topology. See the fix in regulatory_scanner.py's
_forward_slice/_oriented_sequence/_normalize_interval for details.
"""
import pytest
from Bio.Seq import Seq

from sytogen.scripts.regulatory_scanner import scan_promoters, scan_rbs


def _reconstruct(sequence, pred, circular=False):
    """Pull the sequence back out at pred's coordinates, undoing strand."""
    start, end = pred["start"], pred["end"]
    if circular and start >= end:
        sub = sequence[start:] + sequence[:end]
    else:
        sub = sequence[start:end]
    if pred["strand"] == "-":
        sub = str(Seq(sub).reverse_complement())
    return sub


# ---------------------------------------------------------------------------
# scan_promoters
# ---------------------------------------------------------------------------

MINUS35 = "TTGACA"
MINUS10 = "TATAAT"


def _promoter_seq(spacing, prefix="GCGCGCGCGC", suffix="GCGCGCGCGC", minus35=MINUS35):
    return prefix + minus35 + ("C" * spacing) + MINUS10 + suffix


def test_scan_promoters_finds_exact_plus_strand_hit():
    seq = _promoter_seq(spacing=17)
    preds = scan_promoters(seq, topology="linear")

    assert len(preds) == 1
    hit = preds[0]
    assert hit["type"] == "regulatory_promoter"
    assert hit["strand"] == "+"
    assert hit["minus35"] == MINUS35
    assert hit["minus10"] == MINUS10
    assert hit["spacing"] == 17
    assert hit["mismatches"] == 0
    assert hit["score"] == 1.0
    assert _reconstruct(seq, hit) == hit["sequence"]


def test_scan_promoters_finds_minus_strand_hit():
    block = MINUS35 + ("C" * 17) + MINUS10
    rc_block = str(Seq(block).reverse_complement())
    seq = "GCGCGCGCGC" + rc_block + "GCGCGCGCGC"

    preds = scan_promoters(seq, topology="linear")

    assert len(preds) == 1
    hit = preds[0]
    assert hit["strand"] == "-"
    assert hit["minus35"] == MINUS35
    assert hit["minus10"] == MINUS10
    assert hit["mismatches"] == 0
    # The coordinates must map back to the actual embedded block.
    assert _reconstruct(seq, hit) == hit["sequence"]


def test_scan_promoters_detects_wraparound_on_circular_topology():
    # Split the -35...-10 block across the origin: the tail goes at the
    # very start of the linear representation, the head at the very end.
    block = MINUS35 + ("C" * 17) + MINUS10  # 29 bp
    tail, head = block[19:], block[:19]
    seq = tail + "GCGCGCGCGC" + head  # motif spans the origin

    preds = scan_promoters(seq, topology="circular")

    assert len(preds) == 1
    hit = preds[0]
    assert hit["mismatches"] == 0
    # start > end signals a wraparound hit; reconstruction must still work.
    assert hit["start"] > hit["end"]
    assert _reconstruct(seq, hit, circular=True) == hit["sequence"]

    # The same sequence scanned as linear must NOT produce this hit, since
    # the block is genuinely split across what would be two disconnected
    # ends of a linear molecule.
    assert scan_promoters(seq, topology="linear") == []


@pytest.mark.parametrize("spacing,expect_hit", [(15, False), (16, True), (18, True), (19, False)])
def test_scan_promoters_respects_spacing_window(spacing, expect_hit):
    seq = _promoter_seq(spacing=spacing)
    preds = scan_promoters(seq, topology="linear")
    assert (len(preds) == 1) is expect_hit


def test_scan_promoters_allows_up_to_two_mismatches_in_minus35():
    # GTGACA vs TTGACA: 1 mismatch (position 0), still within tolerance.
    seq = _promoter_seq(spacing=17, minus35="GTGACA")
    preds = scan_promoters(seq, topology="linear")

    assert len(preds) == 1
    assert preds[0]["mismatches"] == 1
    assert preds[0]["minus35"] == "GTGACA"
    assert preds[0]["score"] < 1.0


def test_scan_promoters_rejects_three_or_more_mismatches():
    # GGGTCA vs TTGACA: 3 mismatches, exceeds max_mismatches=2.
    seq = _promoter_seq(spacing=17, minus35="GGGTCA")
    assert scan_promoters(seq, topology="linear") == []


def test_scan_promoters_ignores_sequence_with_no_hexamers():
    seq = "ACGT" * 30
    assert scan_promoters(seq, topology="linear") == []


# ---------------------------------------------------------------------------
# scan_rbs
# ---------------------------------------------------------------------------

RBS_MOTIF = "AGGAGG"  # also contains the shorter GGAG/GAGG sub-motifs


def test_scan_rbs_finds_upstream_motif_on_plus_strand_cds():
    # CDS starts (1-based) at 41; upstream window is [coding_start-20,
    # coding_start-4) = [20, 36) in 0-based coordinates. Place the RBS at
    # position 20.
    seq = ("C" * 20) + RBS_MOTIF + ("C" * (60 - 26))
    features = [{"type": "CDS", "strand": "+", "start": 41, "end": 100, "id": "geneA"}]

    preds = scan_rbs(seq, features, topology="linear")

    hit = next(p for p in preds if p["sequence"] == RBS_MOTIF)
    assert hit["type"] == "regulatory_rbs"
    assert hit["strand"] == "+"
    assert hit["associated_feature"] == "geneA"
    assert hit["start"], hit["end"] == (20, 26)
    assert _reconstruct(seq, hit) == hit["sequence"]
    # Also recovers the two shorter overlapping sub-motifs (GGAG, GAGG).
    assert {p["sequence"] for p in preds} == {"AGGAGG", "GGAG", "GAGG"}


def test_scan_rbs_reports_spacing_to_start_codon():
    seq = ("C" * 20) + RBS_MOTIF + ("C" * (60 - 26))
    features = [{"type": "CDS", "strand": "+", "start": 41, "end": 100, "id": "geneA"}]

    preds = scan_rbs(seq, features, topology="linear")
    hit = next(p for p in preds if p["sequence"] == RBS_MOTIF)

    # RBS ends at 0-based position 26; CDS starts at 0-based position 40.
    assert hit["spacing_to_start"] == 40 - 26


def test_scan_rbs_only_considers_gene_like_feature_types():
    seq = ("C" * 20) + RBS_MOTIF + ("C" * (60 - 26))
    non_gene_feature = [{"type": "misc_feature", "strand": "+", "start": 41, "end": 100, "id": "x"}]

    assert scan_rbs(seq, non_gene_feature, topology="linear") == []


def test_scan_rbs_returns_nothing_when_no_motif_present():
    seq = "C" * 60
    features = [{"type": "CDS", "strand": "+", "start": 41, "end": 100, "id": "geneA"}]

    assert scan_rbs(seq, features, topology="linear") == []


def test_scan_rbs_minus_strand_coordinates_on_linear_topology():
    # Regression test for a fixed bug: scan_rbs() used to return
    # coordinates that did not correspond to the matched sequence for
    # '-'-strand features on linear topology (a length-based coordinate
    # flip meant for a different calling convention was double-applied on
    # top of already-absolute coordinates). See git history on
    # regulatory_scanner.py's _map_oriented_interval/_normalize_interval
    # for the fix.
    #
    # CDS on '-' strand, 1-based 1..40 -> coding_start (5' end) = 40.
    # Upstream window (minus strand) = [coding_start+4, coding_start+20)
    #                                 = [44, 60).
    # Embed reverse_complement(RBS_MOTIF) at absolute position [50, 56),
    # i.e. offset 6 within the window, so that after scan_rbs's own
    # reverse-complement of the window the motif reads correctly.
    rc_motif = str(Seq(RBS_MOTIF).reverse_complement())
    window = ("C" * 6) + rc_motif + ("C" * 4)
    seq = ("C" * 44) + window + ("C" * 4)
    features = [{"type": "CDS", "strand": "-", "start": 1, "end": 40, "id": "geneB"}]

    preds = scan_rbs(seq, features, topology="linear")
    hit = next(p for p in preds if p["sequence"] == RBS_MOTIF)

    assert (hit["start"], hit["end"]) == (50, 56)
    assert _reconstruct(seq, hit) == hit["sequence"]


def test_scan_rbs_minus_strand_coordinates_on_circular_topology():
    # Same construction as the linear case above but on circular topology,
    # with no origin wraparound involved. This used to return zero hits: a
    # second, independent bug reverse-*ordered* the window without
    # complementing it, corrupting the window content itself for circular
    # '-'-strand scans (see _forward_slice / _oriented_sequence).
    rc_motif = str(Seq(RBS_MOTIF).reverse_complement())
    window = ("C" * 6) + rc_motif + ("C" * 4)
    seq = ("C" * 44) + window + ("C" * 4)
    features = [{"type": "CDS", "strand": "-", "start": 1, "end": 40, "id": "geneB"}]

    preds = scan_rbs(seq, features, topology="circular")
    hit = next(p for p in preds if p["sequence"] == RBS_MOTIF)

    assert (hit["start"], hit["end"]) == (50, 56)
    assert _reconstruct(seq, hit, circular=True) == hit["sequence"]


def test_scan_rbs_minus_strand_handles_origin_wraparound():
    # Place a '-'-strand gene such that its upstream window (which extends
    # to *higher* coordinates for a minus-strand gene) genuinely straddles
    # the origin on a small circular sequence, and embed the motif so it
    # spans the wrap point itself.
    seq_len = 30
    coding_start = 10  # feature end (1-based end == coding_start here)
    offset = 10
    rc_motif = str(Seq(RBS_MOTIF).reverse_complement())
    bases = list("C" * seq_len)
    window_start = coding_start + 4  # matches default upstream_window minimum
    for i, base in enumerate(rc_motif):
        bases[(window_start + offset + i) % seq_len] = base
    seq = "".join(bases)
    features = [{"type": "CDS", "strand": "-", "start": 1, "end": coding_start, "id": "geneC"}]

    preds = scan_rbs(seq, features, topology="circular")
    hit = next(p for p in preds if p["sequence"] == RBS_MOTIF)

    # Wrapped hit: end < start is the established wraparound convention
    # (also used by scan_promoters for origin-spanning promoter hits).
    assert hit["end"] < hit["start"]
    assert _reconstruct(seq, hit, circular=True) == hit["sequence"]


# ---------------------------------------------------------------------------
# Integration: MotifFinder API endpoint with regulatory_scan enabled
# ---------------------------------------------------------------------------


def test_motiffinder_api_runs_regulatory_scan_when_requested():
    from io import BytesIO

    from sytogen import create_app

    client = create_app().test_client()

    seq = _promoter_seq(spacing=17)
    fasta = f">construct\n{seq}\n".encode()
    motif_table = b"<enz_type>2<rec_seq>GATC<meth_base>A<>"

    response = client.post(
        "/api/motiffinder/run",
        data={
            "source_type": "fasta",
            "sequence_file": (BytesIO(fasta), "construct.fasta"),
            "motif_file": (BytesIO(motif_table), "motifs.txt"),
            "topology": "linear",
            "regulatory_scan": "true",
            "response_format": "json",
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    assert response.json["regulatory_predictions"] >= 1


def test_motiffinder_api_skips_regulatory_scan_by_default():
    from io import BytesIO

    from sytogen import create_app

    client = create_app().test_client()

    seq = _promoter_seq(spacing=17)
    fasta = f">construct\n{seq}\n".encode()
    motif_table = b"<enz_type>2<rec_seq>GATC<meth_base>A<>"

    response = client.post(
        "/api/motiffinder/run",
        data={
            "source_type": "fasta",
            "sequence_file": (BytesIO(fasta), "construct.fasta"),
            "motif_file": (BytesIO(motif_table), "motifs.txt"),
            "topology": "linear",
            "response_format": "json",
            # regulatory_scan omitted -> should default to off
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    assert response.json["regulatory_predictions"] == 0
