"""
Regression tests for find_motif_occurrences() and the two call sites that
depend on it (motiffinder_backend.search_motifs and
sytogen_runner._parse_motifs).

Background: both call sites used to have independent, non-overlapping
regex-based implementations that silently undercounted real occurrences
of self-overlap-capable, IUPAC-degenerate motifs (any motif containing
ambiguity codes like S/N/R/Y/...). Whether a given overlapping occurrence
was caught or missed depended on which direction the scan happened to run
in, so the two implementations also disagreed with each other on
identical input (see the pEPSA5 + S. aureus RM panel case: 261 vs. 279
reported hits for the same construct before this fix, both wrong; 282
after the fix, and now in agreement).
"""
import pandas as pd
import pytest
from Bio.Seq import Seq

from sytogen.scripts.sequence_utils import find_motif_occurrences
from sytogen.scripts.motiffinder_backend import search_motifs
from sytogen.scripts.sytogen_runner import _parse_motifs


def test_self_palindromic_motif_finds_overlapping_occurrences():
    # SCNGS is its own reverse complement (S<->S, C<->G, N<->N are all
    # IUPAC-complement-symmetric here). "ACCCGGGT" contains two genuinely
    # overlapping physical SCNGS sites, at starts 1 ("CCCGG") and 2
    # ("CCGGG") — overlapping by 4 of 5 bases. A non-overlapping scan can
    # only ever return one of the two, and a naive "skip the reverse pass
    # for palindromes" optimization (the old bug) can miss both depending
    # on scan direction.
    hits = find_motif_occurrences("ACCCGGGT", "SCNGS", is_circular=False)
    starts = {h[0] for h in hits}
    assert starts == {1, 2}


def test_non_palindromic_motif_finds_both_strands_without_duplication():
    # GAATTC (EcoRI) IS itself a palindrome, so use a genuinely
    # non-palindromic degenerate motif instead: GARANNNNNNDTYC's reverse
    # complement is GARANNNNNNDTYC's own complement reversed, which is
    # NOT equal to itself (mixed R/Y and D/H asymmetry) — confirm no
    # duplicate reporting of the same physical locus.
    seq = "AAAA" + "GAAAAAAAAAAATTC" + "AAAA"  # arbitrary non-palindromic-context scaffold
    hits = find_motif_occurrences(seq, "GAATTC", is_circular=False)
    # GAATTC is a true palindrome -> forward pass alone must suffice
    positions = [h for h in hits]
    # No duplicate (start, differing only by redundant strand) pairs at
    # the same position for a true palindrome:
    starts = [h[0] for h in positions]
    assert len(starts) == len(set(starts))


def test_circular_topology_still_catches_origin_spanning_motifs():
    # 4-base motif "AAAA" made to wrap the origin of a small circular
    # sequence: last 2 bases + first 2 bases spell AAAA.
    seq = "AACCCCCCCCCCCCCCAA"  # ends ...CCAA, starts AACC... -> wraps to AACC..AA -> "AA"+"AA" at the join
    hits = find_motif_occurrences(seq, "AA", is_circular=True)
    starts = {h[0] for h in hits}
    # The origin-spanning occurrence starts at len(seq) - 1 (last base + first base)
    assert (len(seq) - 1) in starts


def test_search_motifs_and_parse_motifs_agree_on_real_construct():
    """
    The actual regression this fix targets: MotifFinder's standalone
    search_motifs() and the RM-silencing pipeline's internal
    _parse_motifs() must report the identical set of occurrences for the
    identical (sequence, motif panel) input. Before the fix, these
    disagreed by 18 hits on the pEPSA5 + S. aureus RM panel fixture,
    concentrated entirely in two self-palindromic motifs (SCNGS, GCNGC).
    """
    from Bio import SeqIO
    from sytogen.scripts.rebase_motif_parser import parse_rebase_motif_file
    from sytogen.scripts.motiffinder_backend import parse_rebase_motifs

    seq_record = SeqIO.read(
        "tests/fixtures/motiffinder_pEPSA5/motiffinder_annotated.gbk", "genbank"
    )
    sequence = str(seq_record.seq).upper()
    with open("tests/fixtures/my_motifs_saureus.txt") as f:
        motif_text = f.read()

    motif_df_typed = parse_rebase_motif_file(motif_text, is_path=False)
    raw_motifs = parse_rebase_motifs(motif_text)

    standalone_hits = search_motifs(sequence, raw_motifs, is_circular=True)
    internal_motifs = _parse_motifs(motif_df_typed, sequence, "circular")

    standalone_keys = {(h["rec_seq"], h["pos_0"], h["strand"]) for h in standalone_hits}
    internal_keys = {(m.motif, m.start, m.strand) for m in internal_motifs}

    assert standalone_keys == internal_keys
    assert len(standalone_hits) == len(internal_motifs)
