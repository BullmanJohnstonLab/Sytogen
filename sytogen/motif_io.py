"""Framework-neutral motif profile parsing and normalization."""

from __future__ import annotations

import re

import pandas as pd

from sytogen.scripts.rebase_motif_parser import parse_rebase_motif_file

_MIJAMP_HEADER_RE = re.compile(
    r"^#\s*(?P<code>\d+)m(?P<base>[ACGT])\s+modified motifs\b",
    re.IGNORECASE,
)
_MIJAMP_MOTIF_RE = re.compile(
    r"^(?P<prefix>[ACGTURYKMSWBDHVN]+)\((?P<marker>[^()]*)\)(?P<suffix>[ACGTURYKMSWBDHVN]+)$",
    re.IGNORECASE,
)
_IUPAC_COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "U": "A",
    "R": "Y", "Y": "R", "K": "M", "M": "K", "S": "S", "W": "W",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


def _reverse_complement_iupac(seq):
    return "".join(_IUPAC_COMPLEMENT.get(base, base) for base in reversed(str(seq or "").upper()))


def _parse_mijamp_token(token, meth_type=None, meth_base_char=None):
    match = _MIJAMP_MOTIF_RE.match(str(token or "").strip())
    if not match:
        return None
    prefix = match.group("prefix").upper()
    suffix = match.group("suffix").upper()
    marker = match.group("marker").strip().upper()
    meth_base_char = meth_base_char or (re.search(r"([ACGT])$", marker) or [None, None])[1]
    if not meth_base_char:
        return None
    marker_type = re.search(r"(\d+)", marker)
    meth_type = meth_type or (f"m{marker_type.group(1)}{meth_base_char}" if marker_type else f"m{meth_base_char}")
    rec_seq = prefix + meth_base_char + suffix
    plus_position = len(prefix) + 1
    return {
        "rec_seq": rec_seq,
        "enz_type": "1" if "N" in rec_seq or any(c in rec_seq for c in "RYKMSWBDHV") else "2",
        "meth_base": str(plus_position),
        "meth_type": meth_type,
        "comp_meth_base": str(len(rec_seq) - plus_position + 1),
        "comp_meth_type": meth_type,
    }


def parse_mijamp_expected_output(text):
    rows = []
    meth_type = meth_base = None
    for line in (text or "").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            header = _MIJAMP_HEADER_RE.match(stripped)
            if header:
                meth_type = f"m{header.group('code')}{header.group('base').upper()}"
                meth_base = header.group("base").upper()
            continue
        columns = stripped.split("\t")
        parsed = _parse_mijamp_token(columns[0], meth_type, meth_base)
        if parsed:
            parsed["_total_counts"] = columns[2].strip() if len(columns) >= 3 else ""
            rows.append(parsed)

    merged = []
    consumed = set()
    for index, row in enumerate(rows):
        if index in consumed:
            continue
        current = dict(row)
        for candidate_index in range(index + 1, len(rows)):
            candidate = rows[candidate_index]
            if candidate_index in consumed or candidate["_total_counts"] != row["_total_counts"]:
                continue
            if candidate["rec_seq"] != _reverse_complement_iupac(row["rec_seq"]):
                continue
            if str(candidate["meth_base"]).isdigit():
                current["comp_meth_base"] = str(len(row["rec_seq"]) - int(candidate["meth_base"]) + 1)
            current["comp_meth_type"] = candidate["meth_type"]
            consumed.add(candidate_index)
            break
        current.pop("_total_counts", None)
        merged.append(current)
    return pd.DataFrame(merged, columns=["rec_seq", "enz_type", "meth_base", "meth_type", "comp_meth_base", "comp_meth_type"])


def _looks_like_motif(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return False
    cleaned = str(value).strip()
    if not cleaned or cleaned in {"-", "NA", "N/A", "nan", "None"}:
        return False
    letters = re.sub(r"[^ACGTURYKMSWBDHVN]", "", cleaned.upper())
    return 2 <= len(letters) <= 40 and len(letters) / max(1, len(cleaned)) >= 0.6


def parse_motif_text(text):
    mijamp = parse_mijamp_expected_output(text)
    if not mijamp.empty:
        return mijamp
    try:
        dataframe = pd.read_csv(pd.io.common.StringIO(text), sep=None, engine="python")
        if dataframe.empty:
            raise ValueError("empty table")
        columns = {str(column).strip().lower(): column for column in dataframe.columns}
        aliases = ("motif", "rec_seq", "recognition_motif", "recognition_sequence", "sequence", "seq")
        for alias in aliases:
            if alias in columns:
                source = columns[alias]
                return dataframe.rename(columns={source: "motif"}) if source != "motif" else dataframe
        for column in dataframe.columns:
            values = dataframe[column].dropna().astype(str)
            if any(token in str(column).lower() for token in ("motif", "recognition", "sequence", "site", "target")) and any(_looks_like_motif(value) for value in values):
                return dataframe.rename(columns={column: "motif"})
        for column in dataframe.columns:
            if any(_looks_like_motif(value) for value in dataframe[column].dropna().astype(str)):
                return dataframe.rename(columns={column: "motif"})
    except Exception:
        pass
    dataframe = parse_rebase_motif_file(text, is_path=False, drop_unclassified=False)
    if dataframe.empty:
        raise ValueError("Could not parse the restriction motif table.")
    for alias in ("motif", "rec_seq", "recognition_sequence", "sequence", "seq"):
        if alias in dataframe.columns:
            return dataframe.rename(columns={alias: "motif"})
    raise ValueError("Could not find a recognition-sequence field in the motif table.")


def motif_table_records(motif_df):
    aliases = {
        "rec_seq": ("rec_seq", "motif", "recognition_motif", "recognition_sequence", "sequence", "seq"),
        "enz_type": ("enz_type", "type"),
        "meth_base": ("meth_base", "methylated_base_plus", "methylated_base", "methylation_base", "methylation base"),
        "meth_type": ("meth_type", "methylated_base_plus_type", "methylation_type", "methylation", "methylated_type"),
        "comp_meth_base": ("comp_meth_base", "methylated_base_minus", "complementary_methylated_base", "complementary methylated base"),
        "comp_meth_type": ("comp_meth_type", "methylated_base_minus_type", "complementary_methylation_type", "complementary methylation type"),
    }
    columns = {str(column).strip().lower(): column for column in motif_df.columns}

    def value(row, field):
        source = next((columns[name] for name in aliases[field] if name in columns), None)
        return "" if source is None or pd.isna(row[source]) else str(row[source]).strip()

    def base_position(raw, sequence):
        raw = str(raw or "").strip()
        if not raw or raw in {"-", "NA", "N/A", "None", "nan"}:
            return "-"
        if raw.isdigit():
            return raw
        if raw.upper() in {"UNK", "UNKNOWN", "UNKN", "-99", "99"}:
            return "-"
        match = re.search(re.escape(raw.upper()), sequence.upper()) if raw.upper() in "ACGT" else None
        return str(match.start() + 1) if match else raw

    def meth_type(raw):
        text = str(raw or "").strip()
        normalized = text.upper()
        if not text or normalized in {"-", "NA", "N/A", "NONE", "NAN"}:
            return "-"
        if normalized in {"UNK", "UNKNOWN", "UNKN", "-99", "99"}:
            return "Unk"
        if normalized in {"M6A", "6MA", "6", "6A"}:
            return "m6A"
        if normalized in {"M5C", "5MC", "5", "5C"}:
            return "m5C"
        if normalized in {"M4C", "4MC", "4", "4C"}:
            return "m4C"
        return text

    records = []
    for _, row in motif_df.iterrows():
        rec_seq = value(row, "rec_seq").upper()
        if not rec_seq:
            continue
        plus = base_position(value(row, "meth_base"), rec_seq)
        minus = base_position(value(row, "comp_meth_base"), rec_seq)
        plus_type = meth_type(value(row, "meth_type"))
        minus_type = meth_type(value(row, "comp_meth_type"))
        if rec_seq == _reverse_complement_iupac(rec_seq):
            if minus == "-" and plus != "-":
                minus = str(len(rec_seq) - int(plus) + 1) if plus.isdigit() else plus
            if plus == "-" and minus != "-":
                plus = str(len(rec_seq) - int(minus) + 1) if minus.isdigit() else minus
            if minus_type == "-":
                minus_type = plus_type
            if plus_type == "-":
                plus_type = minus_type
        records.append({
            "rec_seq": rec_seq,
            "enz_type": value(row, "enz_type"),
            "meth_base": plus,
            "meth_type": plus_type,
            "comp_meth_base": minus,
            "comp_meth_type": minus_type,
        })
    return records
