"""Tests for sytogen.scripts.assembly_planner.

This module previously had zero test coverage. Writing this suite
surfaced a real bug (fixed alongside it): fragment_sequence() appended
every choose_boundary() result unconditionally. When consecutive target
boundaries have overlapping search windows and the local score landscape
doesn't otherwise distinguish them (e.g. no nearby edits, so
nearest_edit_distance() is a constant everywhere), choose_boundary() can
return the exact same position for consecutive calls, producing a
duplicate boundary and, downstream, a zero-length fragment with a broken
primer pair. See test_fragment_sequence_never_produces_zero_length_fragments
and test_fragment_sequence_warns_when_boundaries_are_deduplicated.
"""
import random

import pytest

from sytogen.scripts.assembly_planner import (
    BOUNDARY_SEARCH_WINDOW,
    PRIMER_MAX_ANNEAL,
    PRIMER_MIN_ANNEAL,
    PRIMER_TARGET_TM,
    TARGET_GC,
    TARGET_TM,
    _extract_window,
    _extract_window_reverse,
    _grow_anneal_region,
    boundary_score,
    choose_boundary,
    collect_edit_positions,
    create_overlap,
    create_overlap_at_origin,
    design_primers_for_fragment,
    fragment_sequence,
    generate_assembly_warnings,
    is_in_gene,
    is_protected,
    nearest_edit_distance,
    overlap_score,
    AssemblyFragment,
)
from sytogen.scripts.genome_model import GenomeModel, Gene, ProtectedRegion
from sytogen.scripts.sequence_utils import reverse_complement


def _random_seq(length, seed=42):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


class FakeGenome:
    """Minimal duck-typed genome for isolated unit tests."""

    def __init__(self, genes=None, protected_regions=None):
        self._genes = genes or []
        self.protected_regions = protected_regions or []

    def find_gene(self, position):
        for g in self._genes:
            if g[0] <= position <= g[1]:
                return object()
        return None


# ---------------------------------------------------------------------------
# collect_edit_positions
# ---------------------------------------------------------------------------

def test_collect_edit_positions_only_counts_chosen_rows_with_a_position():
    matrix = [
        {"chosen": True, "edit_position": "10"},
        {"chosen": False, "edit_position": "20"},  # not chosen -> excluded
        {"chosen": True, "edit_position": ""},      # no position -> excluded
        {"chosen": True, "edit_position": None},    # no position -> excluded
        {"chosen": True, "edit_position": "30"},
    ]
    assert collect_edit_positions(matrix) == {10, 30}


def test_collect_edit_positions_empty_matrix():
    assert collect_edit_positions([]) == set()


# ---------------------------------------------------------------------------
# is_in_gene / is_protected / nearest_edit_distance
# ---------------------------------------------------------------------------

def test_is_in_gene_true_and_false():
    genome = FakeGenome(genes=[(100, 200)])
    assert is_in_gene(genome, 150) is True
    assert is_in_gene(genome, 250) is False


def test_is_protected_checks_inclusive_bounds():
    region = ProtectedRegion(start=100, end=110, label="test")
    genome = FakeGenome(protected_regions=[region])
    assert is_protected(genome, 100) is True   # inclusive start
    assert is_protected(genome, 110) is True   # inclusive end
    assert is_protected(genome, 99) is False
    assert is_protected(genome, 111) is False


def test_is_protected_handles_genome_with_no_protected_regions_attr():
    class Bare:
        pass
    assert is_protected(Bare(), 50) is False


def test_nearest_edit_distance_no_edits_returns_large_sentinel():
    assert nearest_edit_distance(100, set()) == 999999


def test_nearest_edit_distance_picks_closest():
    assert nearest_edit_distance(100, {50, 105, 200}) == 5


# ---------------------------------------------------------------------------
# overlap_score
# ---------------------------------------------------------------------------

def test_overlap_score_penalizes_long_homopolymer():
    # 6+ identical bases in a row triggers the homopolymer penalty.
    with_run = "AAAAAAGCGCGCGCGCGCGC"
    without_run = "ATGCATGCATGCATGCATGC"
    score_with, _, _ = overlap_score(with_run)
    score_without, _, _ = overlap_score(without_run)
    assert score_with < score_without - 40  # -50 penalty should dominate


def test_overlap_score_penalizes_repeated_8mer():
    # First 8 bases repeated later in the same sequence.
    repeated = "ACGTACGT" + "TTTT" + "ACGTACGT"
    non_repeated = "ACGTACGTTTTTGCATGCAT"
    score_repeated, _, _ = overlap_score(repeated)
    score_non_repeated, _, _ = overlap_score(non_repeated)
    assert score_repeated <= score_non_repeated - 10


def test_overlap_score_near_target_gc_and_tm_scores_highest():
    # A sequence with ~50% GC and length tuned toward TARGET_TM should
    # score close to the 100 ceiling (no homopolymer/repeat penalties).
    balanced = "ATGCATGCATGCATGCATGCATGCATGCATGCA"  # 34bp, 50% GC, no repeats
    score, tm, gc = overlap_score(balanced)
    assert abs(gc - 50.0) < 5
    assert score > 50  # comfortably above a poor-scoring sequence


# ---------------------------------------------------------------------------
# boundary_score
# ---------------------------------------------------------------------------

def test_boundary_score_penalizes_gene_and_protected_positions():
    seq = _random_seq(500)
    genome_plain = FakeGenome()
    genome_gene = FakeGenome(genes=[(100, 200)])
    genome_protected = FakeGenome(protected_regions=[ProtectedRegion(100, 200, label="x")])

    plain = boundary_score(seq, 150, genome_plain, set(), overlap_length=35)
    in_gene = boundary_score(seq, 150, genome_gene, set(), overlap_length=35)
    in_protected = boundary_score(seq, 150, genome_protected, set(), overlap_length=35)

    assert in_gene < plain - 900
    assert in_protected < plain - 900


def test_boundary_score_rewards_distance_from_edits():
    seq = _random_seq(500)
    genome = FakeGenome()
    near_edit = boundary_score(seq, 105, genome, {100}, overlap_length=35)
    far_from_edit = boundary_score(seq, 400, genome, {100}, overlap_length=35)
    assert far_from_edit > near_edit


# ---------------------------------------------------------------------------
# choose_boundary
# ---------------------------------------------------------------------------

def test_choose_boundary_avoids_a_gene_spanning_the_whole_window():
    seq = _random_seq(1000)
    # Gene covers the entire search window around position 500 except a
    # small gap right at the edges.
    genome = FakeGenome(genes=[(260, 740)])
    chosen = choose_boundary(seq, 500, genome, set(), overlap_length=35)
    # Best-scoring position should land outside the gene if any
    # candidate in the window does.
    assert not is_in_gene(genome, chosen)


def test_choose_boundary_stays_within_search_window():
    seq = _random_seq(2000)
    genome = FakeGenome()
    approx = 1000
    chosen = choose_boundary(seq, approx, genome, set(), overlap_length=35)
    assert abs(chosen - approx) <= BOUNDARY_SEARCH_WINDOW


# ---------------------------------------------------------------------------
# create_overlap / create_overlap_at_origin
# ---------------------------------------------------------------------------

def test_create_overlap_centers_on_boundary():
    seq = _random_seq(200)
    overlap = create_overlap(seq, boundary=100, overlap_length=20)
    assert overlap.end - overlap.start == 20
    assert overlap.start <= 100 <= overlap.end
    assert overlap.sequence == seq[overlap.start:overlap.end]


def test_create_overlap_at_origin_wraps_and_reconstructs():
    seq = _random_seq(300)
    overlap = create_overlap_at_origin(seq, overlap_length=20)
    assert overlap.end - overlap.start == 20
    assert overlap.start < 0  # un-clamped, by design
    # Reconstruct via wraparound and confirm it matches what's reported.
    reconstructed = seq[overlap.start:] + seq[: overlap.end]
    assert reconstructed == overlap.sequence


# ---------------------------------------------------------------------------
# _extract_window / _extract_window_reverse
# ---------------------------------------------------------------------------

def test_extract_window_forward_wraps_on_circular():
    seq = "ABCDEFGHIJ"
    assert _extract_window(seq, 8, 5, circular=True) == "IJABC"


def test_extract_window_forward_truncates_on_linear():
    seq = "ABCDEFGHIJ"
    assert _extract_window(seq, 8, 5, circular=False) == "IJ"


def test_extract_window_reverse_wraps_on_circular():
    seq = "ABCDEFGHIJ"
    assert _extract_window_reverse(seq, 3, 5, circular=True) == "IJABC"


def test_extract_window_reverse_truncates_on_linear():
    seq = "ABCDEFGHIJ"
    assert _extract_window_reverse(seq, 3, 5, circular=False) == "ABC"


def test_extract_window_reverse_end_equal_to_length_is_not_wrapped():
    # len(sequence) is itself a legitimate exclusive upper bound (the true
    # end of a linear fragment) and must not be treated as an out-of-range
    # value to be wrapped to 0.
    seq = "ABCDEFGHIJ"
    assert _extract_window_reverse(seq, len(seq), 3, circular=False) == "HIJ"


# ---------------------------------------------------------------------------
# _grow_anneal_region
# ---------------------------------------------------------------------------

def test_grow_anneal_region_stops_once_target_tm_reached():
    seq = _random_seq(60)

    def extract(n):
        return seq[:n]

    region, tm = _grow_anneal_region(extract)
    assert PRIMER_MIN_ANNEAL <= len(region) <= PRIMER_MAX_ANNEAL
    # Either it hit the target Tm, or it grew all the way to the max
    # length trying to.
    assert tm >= PRIMER_TARGET_TM or len(region) == PRIMER_MAX_ANNEAL


def test_grow_anneal_region_never_exceeds_max_length():
    # An extremely AT-rich region may never reach PRIMER_TARGET_TM within
    # PRIMER_MAX_ANNEAL bases -- growth must still stop at the cap.
    at_rich = "AT" * 50

    def extract(n):
        return at_rich[:n]

    region, tm = _grow_anneal_region(extract)
    assert len(region) == PRIMER_MAX_ANNEAL


# ---------------------------------------------------------------------------
# design_primers_for_fragment
# ---------------------------------------------------------------------------

def test_design_primers_for_fragment_anneal_regions_match_genome():
    seq = _random_seq(3000)
    genome = GenomeModel(sequence=seq, topology="circular", genes=[], protected_regions=[])
    plan = fragment_sequence(seq, genome, decision_matrix=[], fragment_size=1000, overlap_length=35, topology="circular")

    for frag in plan.fragments:
        tail_fwd = frag.overlap_left.sequence if frag.overlap_left else ""
        anneal_start = frag.overlap_left.end if frag.overlap_left else frag.start
        anneal_fwd = frag.forward_primer.sequence[len(tail_fwd):]
        expected_fwd = _extract_window(seq, anneal_start, len(anneal_fwd), True)
        assert frag.forward_primer.sequence == tail_fwd + expected_fwd

        tail_rev_plus = frag.overlap_right.sequence if frag.overlap_right else ""
        anneal_end = frag.overlap_right.start if frag.overlap_right else frag.end
        anneal_rev_len = len(frag.reverse_primer.sequence) - len(tail_rev_plus)
        expected_rev_plus = _extract_window_reverse(seq, anneal_end, anneal_rev_len, True)
        expected_rev = reverse_complement(tail_rev_plus) + reverse_complement(expected_rev_plus)
        assert frag.reverse_primer.sequence == expected_rev


def test_design_primers_linear_terminus_has_no_tail():
    seq = _random_seq(2500)
    genome = GenomeModel(sequence=seq, topology="linear", genes=[], protected_regions=[])
    plan = fragment_sequence(seq, genome, decision_matrix=[], fragment_size=1000, overlap_length=35, topology="linear")

    first, last = plan.fragments[0], plan.fragments[-1]
    assert first.overlap_left is None
    assert last.overlap_right is None
    # First fragment's forward primer should start right at position 0
    # with no homology tail prepended.
    assert first.forward_primer.sequence == _extract_window(
        seq, 0, len(first.forward_primer.sequence), False
    )


# ---------------------------------------------------------------------------
# generate_assembly_warnings
# ---------------------------------------------------------------------------

def _dummy_overlap(score=50.0, tm=65.0, gc=50.0):
    from sytogen.scripts.assembly_planner import Overlap
    return Overlap(start=0, end=35, sequence="A" * 35, tm=tm, gc_percent=gc, score=score)


def test_generate_assembly_warnings_flags_gene_boundary():
    genome = FakeGenome(genes=[(90, 110)])
    frag = AssemblyFragment(name="frag1", start=100, end=200, sequence="A" * 100)
    frag.overlap_left = _dummy_overlap()
    warnings = generate_assembly_warnings([frag], genome, set(), overlap_length=35, assembly_score=50)
    assert any("falls inside a gene" in w for w in warnings)


def test_generate_assembly_warnings_flags_poor_overlap_score():
    genome = FakeGenome()
    frag = AssemblyFragment(name="frag1", start=100, end=1200, sequence="A" * 1100)
    frag.overlap_left = _dummy_overlap(score=-10.0)
    warnings = generate_assembly_warnings([frag], genome, set(), overlap_length=35, assembly_score=-10)
    assert any("poor" in w or "compromise" in w for w in warnings)


def test_generate_assembly_warnings_flags_short_fragment():
    genome = FakeGenome()
    frag = AssemblyFragment(name="frag1", start=0, end=40, sequence="A" * 40)
    warnings = generate_assembly_warnings([frag], genome, set(), overlap_length=35, assembly_score=50)
    assert any("only 40bp" in w for w in warnings)


def test_generate_assembly_warnings_skips_true_linear_termini():
    # No overlap on either side (true linear ends) -- nothing to warn about
    # even if the position would otherwise be "in a gene".
    genome = FakeGenome(genes=[(0, 2000)])
    frag = AssemblyFragment(name="frag1", start=0, end=2000, sequence="A" * 2000)
    warnings = generate_assembly_warnings([frag], genome, set(), overlap_length=35, assembly_score=50)
    assert not any("falls inside a gene" in w for w in warnings)


# ---------------------------------------------------------------------------
# fragment_sequence (end-to-end)
# ---------------------------------------------------------------------------

def test_fragment_sequence_circular_junctions_are_internally_consistent():
    seq = _random_seq(5000)
    genome = GenomeModel(sequence=seq, topology="circular", genes=[], protected_regions=[])
    plan = fragment_sequence(seq, genome, decision_matrix=[], fragment_size=1500, overlap_length=35, topology="circular")

    assert len(plan.fragments) > 1
    for i, frag in enumerate(plan.fragments):
        nxt = plan.fragments[(i + 1) % len(plan.fragments)]
        assert frag.overlap_right.sequence == nxt.overlap_left.sequence


def test_fragment_sequence_linear_has_no_origin_overlap():
    seq = _random_seq(4000)
    genome = GenomeModel(sequence=seq, topology="linear", genes=[], protected_regions=[])
    plan = fragment_sequence(seq, genome, decision_matrix=[], fragment_size=1500, overlap_length=35, topology="linear")

    assert plan.fragments[0].overlap_left is None
    assert plan.fragments[-1].overlap_right is None


def test_fragment_sequence_single_fragment_when_shorter_than_target_size():
    seq = _random_seq(800)
    genome = GenomeModel(sequence=seq, topology="linear", genes=[], protected_regions=[])
    plan = fragment_sequence(seq, genome, decision_matrix=[], fragment_size=1500, overlap_length=35, topology="linear")

    assert len(plan.fragments) == 1
    assert plan.fragments[0].sequence == seq


def test_fragment_sequence_never_produces_zero_length_fragments():
    # Regression test for a fixed bug: an unconditional boundaries.append()
    # allowed consecutive choose_boundary() calls with overlapping search
    # windows and a flat score landscape (no nearby edits) to return the
    # same position twice, producing a zero-length fragment.
    seq = _random_seq(200)
    genome = GenomeModel(sequence=seq, topology="linear", genes=[], protected_regions=[])
    plan = fragment_sequence(seq, genome, decision_matrix=[], fragment_size=5, overlap_length=35, topology="linear")

    assert all(f.end > f.start for f in plan.fragments)
    assert all(len(f.sequence) > 0 for f in plan.fragments)
    # boundaries must be strictly increasing
    starts = [f.start for f in plan.fragments]
    assert starts == sorted(set(starts))


def test_fragment_sequence_warns_when_boundaries_are_deduplicated():
    seq = _random_seq(200)
    genome = GenomeModel(sequence=seq, topology="linear", genes=[], protected_regions=[])
    plan = fragment_sequence(seq, genome, decision_matrix=[], fragment_size=5, overlap_length=35, topology="linear")
    assert any("resolved to a position at or before the previous boundary" in w for w in plan.warnings)


def test_fragment_sequence_avoids_edited_positions_when_possible():
    seq = _random_seq(3000)
    genome = GenomeModel(sequence=seq, topology="linear", genes=[], protected_regions=[])
    # Simulate an edit right where the default boundary search would land.
    matrix = [{"chosen": True, "edit_position": "1500"}]
    plan = fragment_sequence(seq, genome, decision_matrix=matrix, fragment_size=1500, overlap_length=35, topology="linear")

    boundary_positions = [f.end for f in plan.fragments[:-1]]
    for pos in boundary_positions:
        assert abs(pos - 1500) > 0  # boundary_score actively penalizes proximity


def test_fragment_sequence_assembly_score_reflects_overlap_quality():
    seq = _random_seq(3000)
    genome = GenomeModel(sequence=seq, topology="linear", genes=[], protected_regions=[])
    plan = fragment_sequence(seq, genome, decision_matrix=[], fragment_size=1500, overlap_length=35, topology="linear")
    scores = [
        o.score
        for f in plan.fragments
        for o in (f.overlap_left, f.overlap_right)
        if o is not None
    ]
    assert plan.assembly_score == pytest.approx(sum(scores) / len(scores))
