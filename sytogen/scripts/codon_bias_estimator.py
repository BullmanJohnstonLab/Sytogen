import os
from collections import Counter
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import pandas as pd


# -----------------------------
# PARSING
# -----------------------------
def parse_genbank(path):
    sequences = list(SeqIO.parse(path, "genbank"))
    if not sequences:
        raise ValueError("No sequences found in GenBank file")

    cds_concat = []

    for record in sequences:
        for feat in record.features:
            if feat.type != "CDS":
                continue

            seq = feat.extract(record.seq)
            codon_start = int(feat.qualifiers.get("codon_start", ["1"])[0]) - 1
            seq = feat.extract(record.seq)[codon_start:]

            usable_len = (len(seq) // 3) * 3
            seq = seq[:usable_len]

            if usable_len == 0:
                continue  # skip bad CDS

            cds_concat.append(str(seq).upper())

    if not cds_concat:
        raise ValueError("No valid CDS found")

    return "".join(cds_concat)


# -----------------------------
# CODON COUNTING (FAST)
# -----------------------------
def count_codons(sequence):
    return Counter(
        sequence[i:i+3]
        for i in range(0, len(sequence) - 2, 3)
    )


# -----------------------------
# CODON USAGE
# -----------------------------
def compute_codon_usage(sequence):
    # HARD-LOCK to NCBI standard bacterial genetic code
    codon_table = CodonTable.unambiguous_dna_by_id[11]

    codon_to_aa = codon_table.forward_table.copy()
    for stop in codon_table.stop_codons:
        codon_to_aa[stop] = "Stop"

    counts = count_codons(sequence)

    rows = []

    for codon, aa in codon_to_aa.items():
        rows.append({
            "AA": aa,
            "Codon": codon,
            "Count": counts.get(codon, 0)
        })

    df = pd.DataFrame(rows)

    # normalize within AA
    df["Proportion"] = df.groupby("AA")["Count"].transform(
        lambda x: x / x.sum() if x.sum() > 0 else 0
    )

    # ranking
    df["Ranking"] = df.groupby("AA")["Count"].rank(
        method="dense", ascending=False
    )

    df["Ranking_ratio"] = df.groupby("AA")["Ranking"].transform(
        lambda x: x / len(x)
    )

    return df


# -----------------------------
# MAIN ENTRY (Flask-safe)
# -----------------------------
def run_codon_bias(
    genome_path=None,
    fasta_path=None,
    gff_path=None,
    codon_table=11,
    output_dir="."
):
    if genome_path:
        cds_seq = parse_genbank(genome_path)
        gbk_path = write_normalized_genbank(genome_path, output_dir)
        fasta_path, gff_path = write_fasta_gff_from_genbank(genome_path, output_dir)
    elif fasta_path and gff_path:
        cds_seq = parse_fasta_gff(fasta_path, gff_path)
        fasta_path, gff_path = write_normalized_fasta_gff(fasta_path, gff_path, output_dir)
        gbk_path = write_genbank_from_fasta_gff(fasta_path, gff_path, output_dir)
    else:
        raise ValueError("Invalid input combination")

    df = compute_codon_usage(cds_seq)

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "codon_usage_table.csv")
    df.to_csv(out_path, index=False)

    return {
        "csv": out_path,
        "genbank": gbk_path,
        "fasta": fasta_path,
        "gff": gff_path,
    }


def write_normalized_genbank(genome_path, output_dir):
    out_path = os.path.join(output_dir, "codonbias_input.gbk")
    records = list(SeqIO.parse(genome_path, "genbank"))
    for record in records:
        record.annotations.setdefault("molecule_type", "DNA")
    SeqIO.write(records, out_path, "genbank")
    return out_path


def write_normalized_fasta_gff(fasta_path, gff_path, output_dir):
    fasta_out = os.path.join(output_dir, "codonbias_input.fasta")
    gff_out = os.path.join(output_dir, "codonbias_input.gff3")

    records = list(SeqIO.parse(fasta_path, "fasta"))
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


def _format_gff_attrs(qualifiers):
    attrs = []
    preferred_keys = ["ID", "Name", "gene", "locus_tag", "product", "protein_id"]

    for key in preferred_keys:
        if key in qualifiers:
            value = qualifiers[key]
            if isinstance(value, list):
                value = ",".join(str(v) for v in value)
            attrs.append(f"{key}={value}")

    return ";".join(attrs) if attrs else "."


def write_fasta_gff_from_genbank(genome_path, output_dir):
    fasta_path = os.path.join(output_dir, "codonbias_input.fasta")
    gff_path = os.path.join(output_dir, "codonbias_input.gff3")
    records = list(SeqIO.parse(genome_path, "genbank"))

    SeqIO.write(records, fasta_path, "fasta")

    with open(gff_path, "w") as gff:
        gff.write("##gff-version 3\n")

        for record in records:
            gff.write(f"##sequence-region {record.id} 1 {len(record.seq)}\n")

            for index, feat in enumerate(record.features, start=1):
                if feat.type == "source":
                    continue

                start = int(feat.location.start) + 1
                end = int(feat.location.end)
                strand = "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "."
                attrs = _format_gff_attrs(feat.qualifiers)

                if attrs == ".":
                    attrs = f"ID={feat.type}_{index}"

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


def _parse_gff_attrs(attrs):
    qualifiers = {}
    if attrs == ".":
        return qualifiers

    for item in attrs.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        qualifiers[key] = [value]

    return qualifiers


def write_genbank_from_fasta_gff(fasta_path, gff_path, output_dir):
    out_path = os.path.join(output_dir, "codonbias_input.gbk")
    records = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}

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

            seqid, source, feature_type, start, end, score, strand, phase, attrs = cols
            record = records.get(seqid)
            if not record:
                continue

            strand_value = 1 if strand == "+" else -1 if strand == "-" else None
            feature = SeqFeature(
                FeatureLocation(int(start) - 1, int(end), strand=strand_value),
                type=feature_type,
                qualifiers=_parse_gff_attrs(attrs),
            )

            if source and source != ".":
                feature.qualifiers.setdefault("source", [source])
            if phase and phase != ".":
                feature.qualifiers.setdefault("phase", [phase])

            record.features.append(feature)

    SeqIO.write(list(records.values()), out_path, "genbank")
    return out_path

# -----------------------------
# GFF PARSING (FASTA+GFF MODE)
# -----------------------------
def parse_fasta_gff(fasta_path, gff_path):
    records = {rec.id: rec for rec in SeqIO.parse(fasta_path, "fasta")}
    if not records:
        raise ValueError("No sequences in FASTA")

    cds_seqs = []

    with open(gff_path) as gff:
        for line in gff:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) != 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs = cols

            if feature_type != "CDS":
                continue

            if seqid not in records:
                continue

            start, end = int(start) - 1, int(end)
            seq = records[seqid].seq[start:end]

            if strand == "-":
                seq = seq.reverse_complement()

            if len(seq) % 3 != 0:
                continue

            cds_seqs.append(str(seq).upper())

    if not cds_seqs:
        raise ValueError("No valid CDS found in GFF")

    return "".join(cds_seqs)
