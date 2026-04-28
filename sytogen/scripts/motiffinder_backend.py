"""
motiffinder_backend.py
Pure logic module — no Flask routes.

Exposes:
  parse_rebase_motifs(text)          -> list[dict]
  load_gff3_features(gff_text)       -> list[dict]
  features_overlapping(...)          -> list[str]
  search_motifs(seq_str, motifs, ...) -> list[dict]
  hits_to_gff3(hits, seqid, ...)     -> str
  hits_to_tsv(hits, seqid, ...)      -> str

Constants:
  GFF3_HEADER
  TSV_HEADER
"""

import re
from datetime import datetime

from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# IUPAC ambiguity code → regex character class
# ---------------------------------------------------------------------------
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def iupac_to_regex(seq: str) -> str:
    return "".join(IUPAC.get(b.upper(), b.upper()) for b in seq)


def compile_motif(motif_seq: str):
    return re.compile(iupac_to_regex(motif_seq))


# ---------------------------------------------------------------------------
# REBASE tag-format parser
# ---------------------------------------------------------------------------
def parse_rebase_motifs(text: str) -> list[dict]:
    """
    Parse the MyMotif REBASE tagged format.

    Each record is delimited by <> and contains tags like:
      <enz_type>2<rec_seq>GATCNAC<meth_base>6 ...

    Returns a list of dicts with keys: enz_type, rec_seq, meth_base,
    meth_type, comp_meth_base, comp_meth_type.
    """
    motifs = []
    records = re.split(r"<>", text.strip())
    tag_re = re.compile(r"<([^>]+)>([^<]*)")

    for record in records:
        record = record.strip()
        if not record:
            continue
        fields = {m.group(1): m.group(2).strip() for m in tag_re.finditer(record)}
        rec_seq = fields.get("rec_seq", "").strip().upper()
        if not rec_seq:
            continue
        motifs.append({
            "enz_type":       fields.get("enz_type", ""),
            "rec_seq":        rec_seq,
            "meth_base":      fields.get("meth_base", ""),
            "meth_type":      fields.get("meth_type", ""),
            "comp_meth_base": fields.get("comp_meth_base", ""),
            "comp_meth_type": fields.get("comp_meth_type", ""),
        })
    return motifs


# ---------------------------------------------------------------------------
# GFF3 annotation loader
# ---------------------------------------------------------------------------
def load_gff3_features(gff_text: str) -> list[dict]:
    """
    Minimal GFF3 parser — returns a list of feature dicts.
    Skips comment/pragma lines.
    """
    features = []
    for line in gff_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        features.append({
            "seqid":  parts[0],
            "source": parts[1],
            "type":   parts[2],
            "start":  int(parts[3]),  # 1-based
            "end":    int(parts[4]),  # 1-based inclusive
            "score":  parts[5],
            "strand": parts[6],
            "phase":  parts[7],
            "attrs":  parts[8],
        })
    return features


def features_overlapping(pos_1based: int, motif_len: int,
                          features: list[dict], seqid: str) -> list[str]:
    """Return names/IDs of features overlapping the motif position."""
    hit_end = pos_1based + motif_len - 1
    overlapping = []
    for f in features:
        if f["seqid"] != seqid:
            continue
        if f["start"] <= hit_end and f["end"] >= pos_1based:
            name = None
            for kv in f["attrs"].split(";"):
                kv = kv.strip()
                if kv.lower().startswith("name="):
                    name = kv.split("=", 1)[1]
                    break
                if kv.lower().startswith("id="):
                    name = kv.split("=", 1)[1]
            overlapping.append(name or f["type"])
    return overlapping


# ---------------------------------------------------------------------------
# Core motif search — both strands, circular-aware
# ---------------------------------------------------------------------------
def search_motifs(seq_str: str, motifs: list[dict],
                  is_circular: bool = False) -> list[dict]:
    """
    Search seq_str for all motif hits on both strands.

    For circular sequences, appends the first (motif_len - 1) bases to
    catch wrap-around hits.

    Returns list of hit dicts with keys:
      pos_0, strand, rec_seq, enz_type, meth_base, meth_type,
      comp_meth_base, comp_meth_type, hit_seq
    """
    seq_upper = seq_str.upper()
    seq_len   = len(seq_upper)
    hits      = []

    for motif in motifs:
        rec_seq    = motif["rec_seq"]
        motif_len  = len(rec_seq)
        pattern_fwd = compile_motif(rec_seq)

        rc_seq      = str(Seq(rec_seq).reverse_complement())
        pattern_rev = compile_motif(rc_seq)

        search_str = seq_upper
        if is_circular and motif_len > 1:
            search_str = seq_upper + seq_upper[: motif_len - 1]

        for m in pattern_fwd.finditer(search_str):
            pos = m.start() % seq_len
            hits.append({
                **motif,
                "pos_0":  pos,
                "strand": "+",
                "hit_seq": (
                    seq_upper[pos: pos + motif_len]
                    if pos + motif_len <= seq_len
                    else seq_upper[pos:] + seq_upper[: (pos + motif_len) % seq_len]
                ),
            })

        if rc_seq != rec_seq:
            for m in pattern_rev.finditer(search_str):
                pos = m.start() % seq_len
                hits.append({
                    **motif,
                    "pos_0":  pos,
                    "strand": "-",
                    "hit_seq": (
                        seq_upper[pos: pos + motif_len]
                        if pos + motif_len <= seq_len
                        else seq_upper[pos:] + seq_upper[: (pos + motif_len) % seq_len]
                    ),
                })

    hits.sort(key=lambda h: (h["pos_0"], h["strand"]))
    return hits


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------
GFF3_HEADER = (
    "##gff-version 3\n"
    f"# Generated by MotifFinder — {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}\n"
)

TSV_HEADER = (
    "seqid\tposition_1based\tstrand\tmotif\thit_seq\t"
    "enz_type\tmeth_base\tmeth_type\tcomp_meth_base\tcomp_meth_type\t"
    "overlapping_features\n"
)


def hits_to_gff3(hits: list[dict], seqid: str, seq_len: int,
                 features: list[dict]) -> str:
    lines = [GFF3_HEADER, f"##sequence-region {seqid} 1 {seq_len}\n"]

    for i, h in enumerate(hits, start=1):
        start_1 = h["pos_0"] + 1
        end_1   = h["pos_0"] + len(h["rec_seq"])

        overlaps    = features_overlapping(start_1, len(h["rec_seq"]), features, seqid)
        overlap_str = ",".join(overlaps) if overlaps else "."

        attrs = (
            f"ID=motif_hit_{i:04d};"
            f"Name={h['rec_seq']};"
            f"motif={h['rec_seq']};"
            f"hit_seq={h['hit_seq']};"
            f"enz_type={h['enz_type']};"
            f"meth_base={h['meth_base']};"
            f"meth_type={h['meth_type']};"
            f"comp_meth_base={h['comp_meth_base']};"
            f"comp_meth_type={h['comp_meth_type']};"
            f"overlapping_features={overlap_str}"
        )

        lines.append(
            f"{seqid}\tMotifFinder\tbiological_region\t"
            f"{start_1}\t{end_1}\t.\t{h['strand']}\t.\t{attrs}\n"
        )

    return "".join(lines)


def hits_to_tsv(hits: list[dict], seqid: str,
                features: list[dict]) -> str:
    rows = [TSV_HEADER]
    for h in hits:
        start_1  = h["pos_0"] + 1
        overlaps = features_overlapping(start_1, len(h["rec_seq"]), features, seqid)
        rows.append(
            f"{seqid}\t{start_1}\t{h['strand']}\t{h['rec_seq']}\t"
            f"{h['hit_seq']}\t{h['enz_type']}\t{h['meth_base']}\t"
            f"{h['meth_type']}\t{h['comp_meth_base']}\t{h['comp_meth_type']}\t"
            f"{','.join(overlaps) if overlaps else '.'}\n"
        )
    return "".join(rows)