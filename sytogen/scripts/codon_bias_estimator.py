import os
from pathlib import Path
from collections import Counter

from numpy import record
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import SeqFeature, FeatureLocation


# =========================================================
# FILE TYPE HELPERS
# =========================================================

GENBANK_EXTENSIONS = {
    ".gb",
    ".gbk",
    ".genbank",
    ".gbff",
    ".gbf",
}

FASTA_EXTENSIONS = {
    ".fa",
    ".fasta",
    ".faa",
    ".ffn",
    ".fna",
}

GFF_EXTENSIONS = {
    ".gff",
    ".gff3",
}


def normalize_path(path_or_file):
    """
    Accept either:
      - filesystem path
      - Flask FileStorage object

    Returns filesystem path.
    """

    if hasattr(path_or_file, "filename"):

        filename = os.path.basename(path_or_file.filename)

        temp_dir = "/tmp"
        os.makedirs(temp_dir, exist_ok=True)

        out_path = os.path.join(temp_dir, filename)

        path_or_file.save(out_path)

        return out_path

    return str(path_or_file)


def detect_file_type(path):

    ext = Path(path).suffix.lower()

    if ext in GENBANK_EXTENSIONS:
        return "genbank"

    if ext in FASTA_EXTENSIONS:
        return "fasta"

    if ext in GFF_EXTENSIONS:
        return "gff"

    raise ValueError(f"Unsupported file extension: {ext}")


# =========================================================
# GENBANK PARSING
# =========================================================

def parse_genbank(path):

    path = normalize_path(path)

    sequences = list(SeqIO.parse(path, "genbank"))

    if not sequences:
        raise ValueError("No sequences found in GenBank file")

    cds_concat = []

    accepted_feature_types = {
        "CDS",
        "ORF",
        "orf",
        "open_reading_frame",
        "gene",
    }

    for record in sequences:

        for feat in record.features:

            if feat.type not in accepted_feature_types:
                continue

            try:

                seq = feat.extract(record.seq)

                # -----------------------------------------
                # codon_start handling
                # -----------------------------------------

                codon_start = int(
                    feat.qualifiers.get(
                        "codon_start",
                        ["1"]
                    )[0]
                ) - 1

                if codon_start > 0:
                    seq = seq[codon_start:]

                # -----------------------------------------
                # skip tiny ORFs
                # -----------------------------------------

                if len(seq) < 9:
                    continue

                # -----------------------------------------
                # trim to full codons
                # -----------------------------------------

                usable_len = (len(seq) // 3) * 3

                seq = seq[:usable_len]

                if usable_len == 0:
                    continue

                seq = str(seq).upper()

                # -----------------------------------------
                # remove ambiguous codons
                # -----------------------------------------

                if set(seq) - {"A", "T", "G", "C"}:
                    continue

                cds_concat.append(seq)

            except Exception:
                continue

    if not cds_concat:

        feature_summary = {}

        for record in sequences:
            for feat in record.features:
                feature_summary[feat.type] = (
                    feature_summary.get(feat.type, 0) + 1
                )

        raise ValueError(
            "No valid CDS/ORF features found.\n"
            f"Feature types present: {feature_summary}"
        )

    return "".join(cds_concat)

# =========================================================
# FASTA + GFF PARSING
# =========================================================
def parse_fasta_gff(fasta_path, gff_path):

    fasta_records = load_fasta_records(fasta_path)

    cds_concat = []

    accepted_feature_types = {
        "CDS",
        "ORF",
        "orf",
        "open_reading_frame",
        "gene",
    }

    try:

        with open(gff_path) as handle:

            for line in handle:

                if line.startswith("#"):
                    continue

                cols = line.rstrip("\n").split("\t")

                if len(cols) != 9:
                    continue

                (
                    seqid,
                    source,
                    feature_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attrs,
                ) = cols

                if feature_type not in accepted_feature_types:
                    continue

                if seqid not in fasta_records:
                    continue

                record = fasta_records[seqid]

                try:

                    seq = record.seq[int(start) - 1:int(end)]

                    # -----------------------------------------
                    # skip tiny ORFs
                    # -----------------------------------------

                    if len(seq) < 9:
                        continue

                    # -----------------------------------------
                    # trim to full codons
                    # -----------------------------------------

                    usable_len = (len(seq) // 3) * 3

                    seq = seq[:usable_len]

                    if usable_len == 0:
                        continue

                    seq = str(seq).upper()

                    # -----------------------------------------
                    # remove ambiguous codons
                    # -----------------------------------------

                    if set(seq) - {"A", "T", "G", "C"}:
                        continue

                    cds_concat.append(seq)

                except Exception:
                    continue

    except Exception as e:
        raise ValueError(
            f"Error parsing GFF file: {e}"
        )

    if not cds_concat:
        raise ValueError(
            "No valid CDS/ORF features found in GFF file"
        )

    return "".join(cds_concat)
def load_fasta_records(fasta_path):

    fasta_path = normalize_path(fasta_path)

    records = {}

    try:

        with open(
            fasta_path,
            "r",
            encoding="utf-8",
            errors="ignore",
        ) as handle:

            for rec in SeqIO.parse(handle, "fasta"):

                # full original header
                full_header = rec.description.strip()

                # primary token
                short_id = rec.id.strip()

                # store both forms
                records[short_id] = rec
                records[full_header] = rec

    except Exception as e:
        raise ValueError(
            f"Invalid FASTA file: {e}"
        )

    if not records:
        raise ValueError(
            "No sequences found in FASTA"
        )

    return records

# =========================================================
# CODON COUNTING
# =========================================================

def count_codons(sequence):

    return Counter(
        sequence[i:i + 3]
        for i in range(0, len(sequence) - 2, 3)
    )


# =========================================================
# CODON USAGE
# =========================================================

def compute_codon_usage(sequence, codon_table=11):

    codon_table_obj = CodonTable.unambiguous_dna_by_id[codon_table]

    codon_to_aa = codon_table_obj.forward_table.copy()

    for stop in codon_table_obj.stop_codons:
        codon_to_aa[stop] = "Stop"

    counts = count_codons(sequence)

    rows = []

    for codon, aa in codon_to_aa.items():

        rows.append({
            "AA": aa,
            "Codon": codon,
            "Count": counts.get(codon, 0),
        })

    df = pd.DataFrame(rows)

    df["Proportion"] = (
        df.groupby("AA")["Count"]
        .transform(lambda x: x / x.sum() if x.sum() else 0)
    )

    df["Ranking"] = (
        df.groupby("AA")["Count"]
        .rank(method="dense", ascending=False)
    )

    df["Ranking_ratio"] = (
        df.groupby("AA")["Ranking"]
        .transform(lambda x: x / len(x))
    )

    return df


# =========================================================
# FILE WRITERS
# =========================================================

def write_normalized_genbank(genome_path, output_dir):

    out_path = os.path.join(
        output_dir,
        "codonbias_input.gbk",
    )

    records = list(
        SeqIO.parse(genome_path, "genbank")
    )

    for record in records:
        record.annotations.setdefault(
            "molecule_type",
            "DNA",
        )

    SeqIO.write(records, out_path, "genbank")

    return out_path


def write_normalized_fasta_gff(
    fasta_path,
    gff_path,
    output_dir,
):

    fasta_out = os.path.join(
        output_dir,
        "codonbias_input.fasta",
    )

    gff_out = os.path.join(
        output_dir,
        "codonbias_input.gff3",
    )

    records = list(
        SeqIO.parse(fasta_path, "fasta")
    )

    SeqIO.write(records, fasta_out, "fasta")

    with open(gff_path) as src, open(gff_out, "w") as dst:

        has_header = False

        for line in src:

            if line.startswith("##gff-version"):
                has_header = True
                break

        src.seek(0)

        if not has_header:
            dst.write("##gff-version 3\n")

        for line in src:
            dst.write(line)

    return fasta_out, gff_out


def write_fasta_gff_from_genbank(
    genome_path,
    output_dir,
):

    fasta_path = os.path.join(
        output_dir,
        "codonbias_input.fasta",
    )

    gff_path = os.path.join(
        output_dir,
        "codonbias_input.gff3",
    )

    records = list(
        SeqIO.parse(genome_path, "genbank")
    )

    SeqIO.write(records, fasta_path, "fasta")

    with open(gff_path, "w") as gff:

        gff.write("##gff-version 3\n")

        for record in records:

            for idx, feat in enumerate(record.features):

                if feat.type == "source":
                    continue

                start = int(feat.location.start) + 1
                end = int(feat.location.end)

                strand = "."

                if feat.location.strand == 1:
                    strand = "+"

                elif feat.location.strand == -1:
                    strand = "-"

                attrs = f"ID={feat.type}_{idx}"

                gff.write(
                    "\t".join([
                        record.id,
                        "GenBank",
                        feat.type,
                        str(start),
                        str(end),
                        ".",
                        strand,
                        ".",
                        attrs,
                    ]) + "\n"
                )

    return fasta_path, gff_path


def write_genbank_from_fasta_gff(
    fasta_path,
    gff_path,
    output_dir,
):

    out_path = os.path.join(
        output_dir,
        "codonbias_input.gbk",
    )

    records = load_fasta_records(fasta_path)

    for record in records.values():
        record.annotations["molecule_type"] = "DNA"
        record.features = []

    with open(gff_path) as gff:

        for line in gff:

            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")

            if len(cols) != 9:
                continue

            (   seqid,
                source,
                feature_type,
                start,
                end,
                score,
                strand,
                phase,
                attrs,
            ) = cols

            record = records.get(seqid)

            if record is None: # fallback:
                # try first whitespace-delimited token
                record = records.get(seqid.split()[0])

            if record is None:
                continue

            strand_value = None

            if strand == "+":
                strand_value = 1

            elif strand == "-":
                strand_value = -1

            feature = SeqFeature(
                FeatureLocation(
                    int(start) - 1,
                    int(end),
                    strand=strand_value,
                ),
                type=feature_type,
            )

            records[seqid].features.append(feature)

    SeqIO.write(
        list(records.values()),
        out_path,
        "genbank",
    )

    return out_path


# =========================================================
# MAIN ENTRY
# =========================================================

def run_codon_bias(
    genome_path=None,
    fasta_path=None,
    gff_path=None,
    codon_table=11,
    output_dir=".",
):

    os.makedirs(output_dir, exist_ok=True)

    # -----------------------------------------------------
    # GENBANK INPUT
    # -----------------------------------------------------

    if genome_path:

        genome_path = normalize_path(genome_path)

        file_type = detect_file_type(genome_path)

        if file_type != "genbank":
            raise ValueError(
                "Single-file mode requires "
                "GenBank input (.gb/.gbk/.gbff)"
            )

        cds_seq = parse_genbank(genome_path)

        gbk_path = write_normalized_genbank(
            genome_path,
            output_dir,
        )

        fasta_path, gff_path = (
            write_fasta_gff_from_genbank(
                genome_path,
                output_dir,
            )
        )

    # -----------------------------------------------------
    # FASTA + GFF INPUT
    # -----------------------------------------------------

    elif fasta_path and gff_path:

        fasta_path = normalize_path(fasta_path)
        gff_path = normalize_path(gff_path)

        if detect_file_type(fasta_path) != "fasta":
            raise ValueError(
                "Invalid FASTA input file"
            )

        if detect_file_type(gff_path) != "gff":
            raise ValueError(
                "Invalid GFF input file"
            )

        cds_seq = parse_fasta_gff(
            fasta_path,
            gff_path,
        )

        fasta_path, gff_path = (
            write_normalized_fasta_gff(
                fasta_path,
                gff_path,
                output_dir,
            )
        )

        gbk_path = write_genbank_from_fasta_gff(
            fasta_path,
            gff_path,
            output_dir,
        )

    else:

        raise ValueError(
            "Provide either:\n"
            "  genome_path=<genbank>\n"
            "or\n"
            "  fasta_path=<fasta>, "
            "gff_path=<gff>"
        )

    df = compute_codon_usage(
        cds_seq,
        codon_table=codon_table,
    )

    out_path = os.path.join(
        output_dir,
        "codon_usage_table.csv",
    )

    df.to_csv(out_path, index=False)

    return {
        "csv": out_path,
        "genbank": gbk_path,
        "fasta": fasta_path,
        "gff": gff_path,
    }
