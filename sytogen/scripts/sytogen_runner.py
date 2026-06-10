from Bio import SeqIO
import pandas as pd
import io

def extract_motifs(motif_df):
    """
    Extract motifs from MotifFinder output table.
    """
    # motif_df may be already a DataFrame or a file-like / bytes object
    if isinstance(motif_df, pd.DataFrame):
        df = motif_df.copy()
    else:
        # try to read from bytes / str or file-like
        if hasattr(motif_df, "read"):
            text = motif_df.read()
            if isinstance(text, bytes):
                text = text.decode("utf-8")
        elif isinstance(motif_df, (bytes, str)):
            text = motif_df.decode("utf-8") if isinstance(motif_df, bytes) else motif_df
        else:
            raise ValueError("Unsupported motif_df input")

        df = pd.read_csv(io.StringIO(text), sep="\t")

    df.columns = df.columns.str.strip().str.lower()

    if "motif" not in df.columns:
        raise ValueError("Motif table must contain 'motif' column")
    
    # Extract unique motifs, dropping any NaN values
    motifs = df["motif"].dropna().unique().tolist()
    return motifs

def extract_codon_table(codon_df):
    """
    Extract codon table from CodonBias output table.
    """
    codon_table = {}

    # detect column names flexibly
    codon_col = [c for c in codon_df.columns if "codon" in c.lower()][0]
    aa_col = [c for c in codon_df.columns if "aa" in c.lower()]

    aa_col = aa_col[0] if aa_col else None

    rank_col = None
    for c in codon_df.columns:
        if "rank" in c.lower() or "freq" in c.lower():
            rank_col = c
            break

    for _, row in codon_df.iterrows():
        codon = row[codon_col]

        codon_table[codon] = {
            "Translation": row.get(aa_col) if aa_col else None,
            "Inverse_ranking": row.get(rank_col, 0)
        }

    return codon_table


def run_sytogen_pipeline(
    seq_record,
    codon_df,
    motif_df,
    params=None,
    ):

    if params is None:
        params = {}

    # ============================
    # SEQUENCE
    # ============================

    sequence = str(seq_record.seq)

    annotations = []
    for feat in seq_record.features:
        if feat.type in ["CDS", "gene"]:
            annotations.append([
                int(feat.location.start),
                int(feat.location.end),
                "+" if feat.location.strand == 1 else "-"
            ])

    # ============================
    # MOTIFS AND CODON TABLE
    # ============================

    motifs = extract_motifs(motif_df)
    print("CSV MOTIFS:")
    print(motifs)
    motif_hits = []
    
    for feat in seq_record.features:
        if feat.type == "misc_feature" and "motif" in feat.qualifiers:
            motif_hits.append({
                "motif": feat.qualifiers["motif"][0],
                "start": int(feat.location.start),
                "end": int(feat.location.end)})
        print("GBK MOTIF HITS:")
        print(motif_hits[:5])
    
    codon_table = extract_codon_table(codon_df)

    # ============================
    # RUN SYTOGEN
    # ============================

    from sytogen.scripts.optimizer import run_step1
    from sytogen.scripts.heuristic import run_step2

    inputs = {"sequence": sequence,
              "annotations": annotations,
              "codon_table": codon_table,
              "motifs": motifs}

    step1 = run_step1(inputs, params)
    result_seq = run_step2(step1, params)

    return result_seq