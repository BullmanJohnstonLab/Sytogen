import pytest

from sytogen.scripts.sequence_utils import (
    is_gc_preserving_swap,
    reverse_complement,
    gc_percent,
    longest_homopolymer,
    wallace_tm,
    nearest_neighbor_tm,
)


def test_reverse_complement_handles_standard_dna_sequences():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("acgt") == "acgt"
    assert reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_rejects_unknown_bases():
    with pytest.raises(ValueError):
        reverse_complement("ACGX")


def test_is_gc_preserving_swap_handles_case_insensitive_dna_bases():
    assert is_gc_preserving_swap("g", "c") is True
    assert is_gc_preserving_swap("a", "t") is True
    assert is_gc_preserving_swap("g", "a") is False
    assert is_gc_preserving_swap("n", "n") is False


# ---------------------------------------------------------------------------
# gc_percent / longest_homopolymer
# ---------------------------------------------------------------------------

def test_gc_percent_basic():
    assert gc_percent("GCGC") == 100.0
    assert gc_percent("ATAT") == 0.0
    assert gc_percent("ATGC") == 50.0


def test_gc_percent_empty_string_is_zero_not_a_div_by_zero_error():
    assert gc_percent("") == 0.0


def test_gc_percent_is_case_insensitive():
    assert gc_percent("gcgc") == 100.0


def test_longest_homopolymer_finds_the_longest_run():
    assert longest_homopolymer("ATCGAAAATCG") == 4
    assert longest_homopolymer("ATCG") == 1
    assert longest_homopolymer("") == 0
    assert longest_homopolymer("AAAAAAAAAA") == 10


# ---------------------------------------------------------------------------
# wallace_tm / nearest_neighbor_tm
# ---------------------------------------------------------------------------

def test_wallace_tm_matches_the_classic_formula():
    # 2*(A+T) + 4*(G+C)
    assert wallace_tm("AATT") == 2 * 4
    assert wallace_tm("GGCC") == 4 * 4
    assert wallace_tm("ATGC") == 2 * 2 + 4 * 2


def test_nearest_neighbor_tm_gc_rich_beats_at_rich_of_same_length():
    at_rich = "ATATATATATATATATAT"
    gc_rich = "GCGCGCGCGCGCGCGCGCG"
    assert nearest_neighbor_tm(gc_rich) > nearest_neighbor_tm(at_rich)


def test_nearest_neighbor_tm_increases_with_length():
    short = "ATCGATCG"
    long_ = "ATCGATCGATCGATCGATCG"
    assert nearest_neighbor_tm(long_) > nearest_neighbor_tm(short)


def test_nearest_neighbor_tm_increases_with_monovalent_salt_concentration():
    seq = "ATCGATCGATCGATCGATCG"
    low_salt = nearest_neighbor_tm(seq, monovalent_conc_m=0.01)
    high_salt = nearest_neighbor_tm(seq, monovalent_conc_m=0.5)
    assert high_salt > low_salt


def test_nearest_neighbor_tm_falls_back_to_wallace_for_short_sequences():
    # Below 2bp there's no nearest-neighbor step to compute.
    assert nearest_neighbor_tm("A") == wallace_tm("A")


def test_nearest_neighbor_tm_falls_back_to_wallace_for_ambiguous_bases():
    assert nearest_neighbor_tm("ATCGN") == wallace_tm("ATCGN")


def test_nearest_neighbor_tm_is_case_insensitive():
    assert nearest_neighbor_tm("atcgatcgatcg") == nearest_neighbor_tm("ATCGATCGATCG")
