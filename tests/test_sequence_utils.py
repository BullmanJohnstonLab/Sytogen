import pytest

from sytogen.scripts.sequence_utils import reverse_complement


def test_reverse_complement_handles_standard_dna_sequences():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("acgt") == "acgt"
    assert reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_rejects_unknown_bases():
    with pytest.raises(ValueError):
        reverse_complement("ACGX")
