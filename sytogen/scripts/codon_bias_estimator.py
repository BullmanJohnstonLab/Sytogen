import os
from collections import Counter
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.SeqFeature import FeatureLocation
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
            else:
                print("Skipping CDS:", feat.qualifiers.get("locus_tag"))

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
    elif fasta_path and gff_path:
        cds_seq = parse_fasta_gff(fasta_path, gff_path)
    else:
        raise ValueError("Invalid input combination")

    df = compute_codon_usage(cds_seq)

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, "codon_usage_table.csv")
    df.to_csv(out_path, index=False)

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