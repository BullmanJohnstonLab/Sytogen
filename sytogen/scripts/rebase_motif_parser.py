"""
rebase_motif_parser.py
=======================
Converts REBASE-style tagged restriction-enzyme exports (e.g. my_motifs.txt)
into the plain DataFrame shape sytogen_runner._parse_motifs() expects
(a 'motif' column of IUPAC recognition sequences; no start/end needed —
_parse_motifs will regex-search the sequence on both strands itself).

Input format, one enzyme per record, records separated by "<>":
    <enz_type>2<rec_seq>ATGC<meth_base>C<meth_type>5mC<comp_meth_base>-<comp_meth_type>-<>

Only <rec_seq> is required for motif silencing; <enz_type> is kept as
metadata in case you want to filter (e.g. skip enz_type == -1, which
typically marks enzymes with no confirmed/classified specificity).

REBASE naming convention: an <enz_name> beginning with "M." denotes a
standalone methyltransferase (e.g. "M.Aod1ORFAP"), not a restriction
endonuclease. A solitary M.-prefixed record with no paired restriction
enzyme of the same system poses no cutting risk on transformation -- the
site gets methylated, but there's no cognate enzyme to cut it if it
isn't. Such records are excluded from the restriction-motif set by
default (see drop_methylase_only); this doesn't lose real restriction
risk when the system is a genuine R-M pair, since REBASE represents the
restriction enzyme itself as a separate record (typically unprefixed,
or "R."-prefixed) with its own <rec_seq> entry.
"""

import re
import pandas as pd

_FIELD_RE = {
    "enz_name": re.compile(r"<enz_name>([^<]*)"),
    "enz_type": re.compile(r"<enz_type>([^<]*)"),
    "rec_seq": re.compile(r"<rec_seq>([^<]*)"),
    "meth_base": re.compile(r"<meth_base>([^<]*)"),
    "meth_type": re.compile(r"<meth_type>([^<]*)"),
    "comp_meth_base": re.compile(r"<comp_meth_base>([^<]*)"),
    "comp_meth_type": re.compile(r"<comp_meth_type>([^<]*)"),
}

KNOWN_RESTRICTION_MOTIFS = {
    "ECORI": "GAATTC",
    "BAMHI": "GGATCC",
    "HINDIII": "AAGCTT",
    "NOTI": "GCGGCCGC",
    "XHOI": "CTCGAG",
    "SALI": "GTCGAC",
    "XBAI": "TCTAGA",
    "NCOI": "CCATGG",
    "KPNI": "GGTACC",
    "PSTI": "CTGCAG",
    "SPHI": "GCATGC",
}

KNOWN_RESTRICTION_MOTIFS = {
    "ECORI": "GAATTC",
    "BAMHI": "GGATCC",
    "HINDIII": "AAGCTT",
    "NOTI": "GCGGCCGC",
    "XHOI": "CTCGAG",
    "SALI": "GTCGAC",
    "XBAI": "TCTAGA",
    "NCOI": "CCATGG",
    "KPNI": "GGTACC",
    "PSTI": "CTGCAG",
    "SPHI": "GCATGC",
}


def _normalize_motif_sequence(seq):
    return re.sub(r"[^A-Za-z]", "", str(seq or "").upper())


def _parse_known_motif_line(line):
    cleaned = (line or "").split("#", 1)[0].strip()
    if not cleaned:
        return None

    if "<" in cleaned or ">" in cleaned:
        return None

    parts = [part for part in re.split(r"[\s:=,]+", cleaned) if part]
    if not parts:
        return None

    first = parts[0].strip().upper()
    if first in KNOWN_RESTRICTION_MOTIFS:
        motif = KNOWN_RESTRICTION_MOTIFS[first]
        suffix = " ".join(parts[1:])
        if suffix:
            seq = _normalize_motif_sequence(suffix)
            if re.fullmatch(r"[ACGTURYKMSWBDHVN]+", seq):
                motif = seq
        return {"motif": motif, "enz_type": ""}

    candidate_seq = _normalize_motif_sequence(first)
    if re.fullmatch(r"[ACGTURYKMSWBDHVN]+", candidate_seq):
        return {"motif": candidate_seq, "enz_type": ""}

    if len(parts) > 1:
        candidate_seq = _normalize_motif_sequence(parts[-1])
        if re.fullmatch(r"[ACGTURYKMSWBDHVN]+", candidate_seq):
            return {"motif": candidate_seq, "enz_type": ""}

    return None


def parse_known_motif_entries(text):
    rows = []
    for line in (text or "").splitlines():
        parsed = _parse_known_motif_line(line)
        if parsed:
            rows.append(parsed)
    return rows


def _normalize_motif_sequence(seq):
    return re.sub(r"[^A-Za-z]", "", str(seq or "").upper())


def _parse_known_motif_line(line):
    cleaned = (line or "").split("#", 1)[0].strip()
    if not cleaned:
        return None

    if "<" in cleaned or ">" in cleaned:
        return None

    parts = [part for part in re.split(r"[\s:=,]+", cleaned) if part]
    if not parts:
        return None

    first = parts[0].strip().upper()
    if first in KNOWN_RESTRICTION_MOTIFS:
        motif = KNOWN_RESTRICTION_MOTIFS[first]
        suffix = " ".join(parts[1:])
        if suffix:
            seq = _normalize_motif_sequence(suffix)
            if re.fullmatch(r"[ACGTURYKMSWBDHVN]+", seq):
                motif = seq
        return {"motif": motif, "enz_type": ""}

    candidate_seq = _normalize_motif_sequence(first)
    if re.fullmatch(r"[ACGTURYKMSWBDHVN]+", candidate_seq):
        return {"motif": candidate_seq, "enz_type": ""}

    if len(parts) > 1:
        candidate_seq = _normalize_motif_sequence(parts[-1])
        if re.fullmatch(r"[ACGTURYKMSWBDHVN]+", candidate_seq):
            return {"motif": candidate_seq, "enz_type": ""}

    return None


def parse_known_motif_entries(text):
    rows = []
    for line in (text or "").splitlines():
        parsed = _parse_known_motif_line(line)
        if parsed:
            rows.append(parsed)
    return rows


def parse_rebase_motif_file(path_or_text, is_path=True, drop_unclassified=False, drop_methylase_only=True):
    """
    Parameters
    ----------
    path_or_text : str
        File path, or raw text if is_path=False.
    is_path : bool
        Whether path_or_text is a filesystem path (True) or raw text (False).
    drop_unclassified : bool
        If True, drop records where enz_type == "-1" (no confirmed
        specificity in REBASE's convention). Default False — SyToGen still
        benefits from silencing sites even for enzymes of unknown type.
    drop_methylase_only : bool
        If True (default), drop records whose enz_name begins with "M."
        (case-insensitive) — REBASE's convention for a standalone
        methyltransferase with no paired restriction enzyme. See the
        module docstring for why this is safe: a genuine R-M pair's
        restriction enzyme is represented as its own separate record and
        is unaffected.

    Returns
    -------
    pd.DataFrame with columns: 'motif', 'enz_name', 'enz_type',
    'meth_base', 'meth_type', 'comp_meth_base', 'comp_meth_type'. Ready to
    pass straight into sytogen_runner._parse_motifs(df, sequence) (extra
    columns like 'enz_name' are simply ignored there).
    """
    if is_path:
        with open(path_or_text, "r") as f:
            content = f.read()
    else:
        content = path_or_text

    rows = []
    for record in content.split("<>"):
        if not record.strip():
            continue

        rec_seq_match = _FIELD_RE["rec_seq"].search(record)
        if not rec_seq_match:
            continue  # no recognition sequence in this record — nothing to silence

        motif = rec_seq_match.group(1).strip().upper()
        if not motif or motif == "-":
            continue

        enz_name_match = _FIELD_RE["enz_name"].search(record)
        enz_name = enz_name_match.group(1).strip() if enz_name_match else ""

        if drop_methylase_only and enz_name.upper().startswith("M."):
            continue  # standalone methyltransferase — no paired restriction enzyme, no cutting risk

        enz_type_match = _FIELD_RE["enz_type"].search(record)
        enz_type = enz_type_match.group(1).strip() if enz_type_match else ""

        meth_base_match = _FIELD_RE["meth_base"].search(record)
        meth_base = meth_base_match.group(1).strip() if meth_base_match else ""
        meth_type_match = _FIELD_RE["meth_type"].search(record)
        meth_type = meth_type_match.group(1).strip() if meth_type_match else ""
        comp_meth_base_match = _FIELD_RE["comp_meth_base"].search(record)
        comp_meth_base = comp_meth_base_match.group(1).strip() if comp_meth_base_match else ""
        comp_meth_type_match = _FIELD_RE["comp_meth_type"].search(record)
        comp_meth_type = comp_meth_type_match.group(1).strip() if comp_meth_type_match else ""

        if drop_unclassified and enz_type == "-1":
            continue

        rows.append({
            "motif": motif,
            "enz_name": enz_name,
            "enz_type": enz_type,
            "meth_base": meth_base,
            "meth_type": meth_type,
            "comp_meth_base": comp_meth_base,
            "comp_meth_type": comp_meth_type,
        })

    if not rows:
        rows = parse_known_motif_entries(content)

    return pd.DataFrame(
        rows,
        columns=["motif", "enz_name", "enz_type", "meth_base", "meth_type", "comp_meth_base", "comp_meth_type"],
    )


if __name__ == "__main__":
    import sys
    df = parse_rebase_motif_file(sys.argv[1])
    print(df.to_string(index=False))
    print(f"\n{len(df)} motifs parsed")
