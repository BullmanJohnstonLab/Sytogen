"""Conservative, report-only bacterial regulatory sequence scans."""

from __future__ import annotations

import re
from typing import Iterable

from Bio.Seq import Seq

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

DEFAULT_PROFILE = {
    "name": "generic_bacterial",
    "rbs": {
        "motifs": ["AGGAGG", "GGAG", "GAGG"],
        "upstream_window": (4, 20),
    },
    "promoter": {
        "minus35": "TTGACA",
        "minus10": "TATAAT",
        "max_mismatches": 2,
        "spacing": (16, 18),
    },
}


def _compile(pattern: str):
    return re.compile("".join(IUPAC.get(base, base) for base in pattern.upper()))


def _mismatches(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left.upper(), right.upper()))


def _oriented_sequence(sequence: str, start: int, end: int, strand: int, circular: bool) -> str:
    length = len(sequence)
    if circular:
        return "".join(sequence[index % length] for index in range(start, end))[::1] if strand == 1 else "".join(sequence[index % length] for index in range(end - 1, start - 1, -1))
    if start < 0 or end > length:
        return ""
    segment = sequence[start:end]
    return segment if strand == 1 else str(Seq(segment).reverse_complement())


def _map_oriented_interval(start: int, end: int, strand: int, length: int, circular: bool):
    if strand == 1:
        return start % length if circular else start, end % length if circular else end
    if circular:
        return (end - 1) % length, (start + 1) % length
    return length - end, length - start


def scan_rbs(sequence: str, features: Iterable[dict], topology: str = "linear", profile: dict | None = None) -> list[dict]:
    """Find RBS-like motifs upstream of annotated CDS features."""
    profile = profile or DEFAULT_PROFILE
    sequence = str(sequence).upper()
    circular = topology == "circular"
    length = len(sequence)
    predictions = []
    for feature in features:
        if str(feature.get("type", "")).upper() not in {"CDS", "ORF", "GENE", "MARKER"}:
            continue
        strand = 1 if feature.get("strand", "+") != "-" else -1
        feature_start = int(feature["start"]) - 1
        feature_end = int(feature["end"])
        coding_start = feature_start if strand == 1 else feature_end
        minimum, maximum = profile.get("rbs", {}).get("upstream_window", (4, 20))
        upstream_start = coding_start - maximum if strand == 1 else coding_start + minimum
        upstream_end = coding_start - minimum if strand == 1 else coding_start + maximum
        window = _oriented_sequence(sequence, upstream_start, upstream_end, strand, circular)
        if not window:
            continue
        for motif in profile.get("rbs", {}).get("motifs", []):
            for match in _compile(motif).finditer(window):
                oriented_start = upstream_start + match.start() if strand == 1 else upstream_end - match.end()
                oriented_end = oriented_start + len(motif)
                start, end = _map_oriented_interval(oriented_start, oriented_end, strand, length, circular)
                predictions.append({
                    "type": "regulatory_rbs",
                    "start": start,
                    "end": end,
                    "strand": "+" if strand == 1 else "-",
                    "sequence": match.group(0),
                    "score": 1.0,
                    "confidence": "candidate",
                    "profile": profile.get("name", "custom"),
                    "associated_feature": feature.get("id") or feature.get("name") or feature.get("attrs", ""),
                    "spacing_to_start": (len(window) - match.end()) + minimum,
                })
    return sorted(predictions, key=lambda row: (row["start"], row["strand"]))


def scan_promoters(sequence: str, topology: str = "linear", profile: dict | None = None) -> list[dict]:
    """Find paired -35/-10 promoter-like elements on both strands."""
    profile = profile or DEFAULT_PROFILE
    sequence = str(sequence).upper()
    circular = topology == "circular"
    length = len(sequence)
    promoter = profile.get("promoter", {})
    minus35 = promoter.get("minus35", "TTGACA").upper()
    minus10 = promoter.get("minus10", "TATAAT").upper()
    max_mismatches = int(promoter.get("max_mismatches", 2))
    min_spacing, max_spacing = promoter.get("spacing", (16, 18))
    scan_sequence = sequence + sequence[: len(minus35) + len(minus10) + max_spacing] if circular else sequence
    predictions = []
    for strand in (1, -1):
        oriented = scan_sequence if strand == 1 else str(Seq(scan_sequence).reverse_complement())
        for start35 in range(0, len(oriented) - len(minus35) - min_spacing - len(minus10) + 1):
            hit35 = oriented[start35:start35 + len(minus35)]
            mismatch35 = _mismatches(hit35, minus35)
            if mismatch35 > max_mismatches:
                continue
            for spacing in range(min_spacing, max_spacing + 1):
                start10 = start35 + len(minus35) + spacing
                hit10 = oriented[start10:start10 + len(minus10)]
                if len(hit10) != len(minus10):
                    continue
                mismatch10 = _mismatches(hit10, minus10)
                if mismatch10 > max_mismatches:
                    continue
                oriented_end = start10 + len(minus10)
                if strand == 1:
                    raw_start, raw_end = start35, oriented_end
                else:
                    raw_start, raw_end = length - oriented_end, length - start35
                if not circular and (raw_start < 0 or raw_end > length):
                    continue
                start = raw_start % length if circular else raw_start
                end = raw_end % length if circular else raw_end
                predictions.append({
                    "type": "regulatory_promoter",
                    "start": start,
                    "end": end,
                    "strand": "+" if strand == 1 else "-",
                    "sequence": oriented[start35:oriented_end],
                    "minus35": hit35,
                    "minus10": hit10,
                    "spacing": spacing,
                    "mismatches": mismatch35 + mismatch10,
                    "score": 1.0 - ((mismatch35 + mismatch10) / (len(minus35) + len(minus10))),
                    "confidence": "candidate",
                    "profile": profile.get("name", "custom"),
                })
    unique = {(row["start"], row["end"], row["strand"], row["minus35"], row["minus10"]): row for row in predictions}
    return sorted(unique.values(), key=lambda row: (row["start"], row["strand"]))


REGULATORY_TSV_HEADER = (
    "seqid\ttype\tstart\tend\tstrand\tsequence\tscore\tconfidence\t"
    "profile\tassociated_feature\tspacing_to_start\tminus35\tminus10\t"
    "spacing\tmismatches\n"
)


def predictions_to_tsv(predictions: list[dict]) -> str:
    fields = REGULATORY_TSV_HEADER.rstrip("\n").split("\t")
    rows = [REGULATORY_TSV_HEADER]
    for prediction in predictions:
        rows.append("\t".join(str(prediction.get(field, "")) for field in fields) + "\n")
    return "".join(rows)


def predictions_to_gff3(predictions: list[dict], sequence_lengths: dict[str, int]) -> str:
    lines = ["##gff-version 3\n"]
    for seqid, length in sequence_lengths.items():
        lines.append(f"##sequence-region {seqid} 1 {length}\n")
    for index, prediction in enumerate(predictions, start=1):
        attrs = [
            f"ID=regulatory_prediction_{index:04d}",
            f"Name={prediction['type']}",
            f"sequence={prediction.get('sequence', '')}",
            f"score={prediction.get('score', '')}",
            f"confidence={prediction.get('confidence', '')}",
            f"profile={prediction.get('profile', '')}",
        ]
        for key in ("associated_feature", "spacing_to_start", "minus35", "minus10", "spacing", "mismatches"):
            if prediction.get(key, "") != "":
                attrs.append(f"{key}={prediction[key]}")
        lines.append(
            f"{prediction['seqid']}\tMotifFinder\t{prediction['type']}\t"
            f"{int(prediction['start']) + 1}\t{int(prediction['end'])}\t"
            f"{prediction.get('score', '.')}\t{prediction['strand']}\t.\t"
            f"{';'.join(attrs)}\n"
        )
    return "".join(lines)
