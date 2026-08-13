"""Framework-neutral sequence and table input helpers."""

from __future__ import annotations

import io
import os

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

from sytogen.scripts.motiffinder_backend import load_gff3_features

MAX_CONSTRUCT_LENGTH = 3_000_000


def allowed_extension(filename, allowed):
    return os.path.splitext(filename)[1].lower() in allowed


def validate_construct_size(record):
    """Raise a clear validation error when a construct exceeds 3000 kb."""
    sequence_length = len(record.seq)
    if sequence_length > MAX_CONSTRUCT_LENGTH:
        name = record.id or record.name or "Uploaded construct"
        raise ValueError(
            f"{name} is {sequence_length:,} bp. The maximum supported construct "
            f"size is {MAX_CONSTRUCT_LENGTH:,} bp (3000 kb)."
        )
    return record


def _parse_gff3_attrs(attrs):
    qualifiers = {}
    for part in (attrs or "").split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        qualifiers.setdefault(key.strip(), []).append(value.strip())
    return qualifiers


def build_record_from_fasta_gff(fasta_text, gff_text, topology="circular"):
    """Build a single annotated SeqRecord from FASTA text and GFF3 text."""
    records = list(SeqIO.parse(io.StringIO(fasta_text), "fasta"))
    if len(records) != 1:
        raise ValueError(
            f"FASTA file must contain exactly one sequence, found {len(records)}."
        )
    record = records[0]
    seqid = record.id or record.name
    raw_features = load_gff3_features(gff_text)
    matching = [feature for feature in raw_features if feature.get("seqid") == seqid]
    if not matching and raw_features:
        distinct_seqids = {feature.get("seqid") for feature in raw_features}
        if len(distinct_seqids) == 1:
            matching = raw_features
        else:
            raise ValueError(
                f"No GFF3 features matched FASTA sequence id '{seqid}' "
                f"(GFF3 contains: {sorted(distinct_seqids)})."
            )

    record.features = []
    for feature in matching:
        strand_symbol = feature.get("strand", ".")
        strand = 1 if strand_symbol == "+" else -1 if strand_symbol == "-" else 1
        record.features.append(
            SeqFeature(
                FeatureLocation(
                    int(feature["start"]) - 1,
                    int(feature["end"]),
                    strand=strand,
                ),
                type=feature.get("type", "misc_feature"),
                qualifiers=_parse_gff3_attrs(feature.get("attrs", "")),
            )
        )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = topology
    return validate_construct_size(record)


def read_uploaded_table(file_storage):
    """Read an uploaded CSV/TSV file into a DataFrame."""
    text = file_storage.stream.read().decode("utf-8-sig")
    return pd.read_csv(io.StringIO(text), sep=None, engine="python")
