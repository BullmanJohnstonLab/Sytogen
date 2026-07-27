import pytest

from sytogen.scripts.sequence_utils import is_gc_preserving_swap, reverse_complement


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
