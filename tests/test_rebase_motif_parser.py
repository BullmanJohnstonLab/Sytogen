"""Tests for sytogen.scripts.rebase_motif_parser.

Focused on the enz_name / methylase-only filtering fix: parse_rebase_
motif_file() previously extracted <rec_seq> from any REBASE-tagged
record regardless of <enz_name>, so a standalone methyltransferase
record (REBASE's "M." prefix convention, e.g. "M.Aod1ORFAP") was
included in the "motifs to silence" set identically to an actual
restriction enzyme -- even though a solitary M.-prefixed record with no
paired restriction enzyme poses no cutting risk on transformation.
"""
import pandas as pd

from sytogen.scripts.rebase_motif_parser import parse_rebase_motif_file


def _record(enz_name=None, enz_type="2", rec_seq="GATC", meth_base="2", meth_type="6"):
    parts = []
    if enz_name is not None:
        parts.append(f"<enz_name>{enz_name}")
    parts.append(f"<enz_type>{enz_type}")
    parts.append(f"<rec_seq>{rec_seq}")
    parts.append(f"<meth_base>{meth_base}")
    parts.append(f"<meth_type>{meth_type}")
    return "".join(parts) + "<>"


def test_methylase_only_record_is_excluded_by_default():
    text = _record(enz_name="M.OrphanMTase", rec_seq="GATC")
    df = parse_rebase_motif_file(text, is_path=False)
    assert df.empty


def test_methylase_prefix_check_is_case_insensitive():
    text = _record(enz_name="m.lowercaseprefix", rec_seq="GATC")
    df = parse_rebase_motif_file(text, is_path=False)
    assert df.empty


def test_genuine_restriction_enzyme_record_is_kept():
    text = _record(enz_name="EcoKI", rec_seq="AACNNNNNNGTGC")
    df = parse_rebase_motif_file(text, is_path=False)
    assert list(df["motif"]) == ["AACNNNNNNGTGC"]
    assert list(df["enz_name"]) == ["EcoKI"]


def test_paired_system_keeps_the_restriction_enzyme_when_methylase_is_dropped():
    # A realistic paired system: REBASE represents the methylase and the
    # restriction enzyme as two separate records sharing a rec_seq. Only
    # the M.-prefixed one should be dropped.
    text = (
        _record(enz_name="M.EcoKI", rec_seq="AACNNNNNNGTGC")
        + _record(enz_name="EcoKI", rec_seq="AACNNNNNNGTGC")
    )
    df = parse_rebase_motif_file(text, is_path=False)
    assert len(df) == 1
    assert df.iloc[0]["enz_name"] == "EcoKI"
    assert df.iloc[0]["motif"] == "AACNNNNNNGTGC"


def test_records_with_no_enz_name_are_not_affected():
    # enz_name is optional in REBASE exports; absence must not be
    # mistaken for a methylase-only record.
    text = _record(enz_name=None, rec_seq="GAATTC")
    df = parse_rebase_motif_file(text, is_path=False)
    assert list(df["motif"]) == ["GAATTC"]


def test_drop_methylase_only_can_be_disabled():
    text = _record(enz_name="M.OrphanMTase", rec_seq="GATC")
    df = parse_rebase_motif_file(text, is_path=False, drop_methylase_only=False)
    assert list(df["motif"]) == ["GATC"]
    assert list(df["enz_name"]) == ["M.OrphanMTase"]


def test_manually_entered_known_motifs_are_unaffected_by_methylase_filter():
    # The manual/known-enzyme-name entry path (used for e.g. domesticating
    # a Golden Gate site) is a completely separate code path with no
    # enz_name field at all -- it must never be caught by this filter.
    df = parse_rebase_motif_file("EcoRI", is_path=False)
    assert list(df["motif"]) == ["GAATTC"]


def test_output_includes_enz_name_column():
    text = _record(enz_name="EcoKI", rec_seq="AACNNNNNNGTGC")
    df = parse_rebase_motif_file(text, is_path=False)
    assert "enz_name" in df.columns
