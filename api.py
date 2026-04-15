from flask import Blueprint, request, send_file, abort
import io
from Bio import SeqIO

from sytogen.scripts.motiffinder_backend import (
    parse_rebase_motifs,
    search_motifs,
    hits_to_gff3,
    hits_to_tsv,
    load_gff3_features,
    GFF3_HEADER,
    TSV_HEADER,
)

# ------------------------------------------------------------------
# Blueprint
# ------------------------------------------------------------------

api = Blueprint("api", __name__)

# ------------------------------------------------------------------
# MotifFinder endpoint
# ------------------------------------------------------------------

@api.route("/motiffinder/run", methods=["POST"])
def run_motiffinder_sync():
    """
    Run MotifFinder synchronously and return a ZIP of results.

    Expects multipart/form-data with:
      sequence_file   — GenBank (.gb/.gbk) or FASTA (.fa/.fasta/.fna)
      motif_file      — MyMotif REBASE tagged format (.txt)
      source_type     — "genbank" | "fasta"
      annotation_file — optional GFF3 (only when source_type == "fasta")

    Returns:
      ZIP containing motiffinder_results.gff3 + motiffinder_summary.tsv
    """

    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    missing = []
    if "sequence_file" not in request.files:
        missing.append("sequence_file")
    if "motif_file" not in request.files:
        missing.append("motif_file")
    if missing:
        abort(400, f"Missing required files: {', '.join(missing)}")

    source_type = request.form.get("source_type", "").lower()
    if source_type not in ("genbank", "fasta"):
        abort(400, "source_type must be 'genbank' or 'fasta'")

    seq_file   = request.files["sequence_file"]
    motif_file = request.files["motif_file"]
    ann_file   = request.files.get("annotation_file")

    # ------------------------------------------------------------------
    # Parse motifs
    # ------------------------------------------------------------------
    motif_text = motif_file.read().decode("utf-8", errors="replace")
    motifs = parse_rebase_motifs(motif_text)
    if not motifs:
        abort(400, "No valid motifs found in motif_file")

    # ------------------------------------------------------------------
    # Parse sequence + features
    # ------------------------------------------------------------------
    seq_text = seq_file.read().decode("utf-8", errors="replace")
    features = []
    records = []

    if source_type == "genbank":
        try:
            records = list(SeqIO.parse(io.StringIO(seq_text), "genbank"))
        except Exception as e:
            abort(400, f"Failed to parse GenBank file: {e}")

        for rec in records:
            for feat in rec.features:
                loc   = feat.location
                seqid = rec.id or rec.name
                attrs = (
                    f"ID={feat.type}_{int(loc.start)+1}_{int(loc.end)};"
                    f"Name={feat.qualifiers.get('gene', [feat.type])[0]}"
                )
                features.append({
                    "seqid":  seqid,
                    "source": "GenBank",
                    "type":   feat.type,
                    "start":  int(loc.start) + 1,
                    "end":    int(loc.end),
                    "score":  ".",
                    "strand": "+" if loc.strand == 1 else "-" if loc.strand == -1 else ".",
                    "phase":  ".",
                    "attrs":  attrs,
                })

    else:  # fasta
        try:
            records = list(SeqIO.parse(io.StringIO(seq_text), "fasta"))
        except Exception as e:
            abort(400, f"Failed to parse FASTA file: {e}")

        if ann_file:
            gff_text = ann_file.read().decode("utf-8", errors="replace")
            features = load_gff3_features(gff_text)

    if not records:
        abort(400, "No sequences found in sequence_file")

    # ------------------------------------------------------------------
    # Run search across all records, assemble outputs
    # ------------------------------------------------------------------
    all_gff3_parts = [GFF3_HEADER]
    all_tsv_parts  = [TSV_HEADER]

    for rec in records:
        seqid   = rec.id or rec.name or "unknown"
        seq_str = str(rec.seq)
        seq_len = len(seq_str)

        is_circular = False
        if source_type == "genbank":
            is_circular = rec.annotations.get("topology", "").lower() == "circular"

        hits = search_motifs(seq_str, motifs, is_circular=is_circular)

        rec_features = [f for f in features if f["seqid"] == seqid]

        # GFF3 — strip repeated header lines from per-record call
        gff3_body  = hits_to_gff3(hits, seqid, seq_len, rec_features)
        body_lines = [l for l in gff3_body.splitlines(keepends=True) if not l.startswith("#")]
        all_gff3_parts.append(f"##sequence-region {seqid} 1 {seq_len}\n")
        all_gff3_parts.extend(body_lines)

        # TSV — skip repeated header row
        tsv_body = hits_to_tsv(hits, seqid, rec_features)
        all_tsv_parts.extend(tsv_body.splitlines(keepends=True)[1:])

    final_gff3 = "".join(all_gff3_parts)
    final_tsv  = "".join(all_tsv_parts)

    # ------------------------------------------------------------------
    # Pack into ZIP and return
    # ------------------------------------------------------------------
    import zipfile
    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("motiffinder_results.gff3", final_gff3)
        zf.writestr("motiffinder_summary.tsv",  final_tsv)
    zip_buf.seek(0)

    return send_file(
        zip_buf,
        mimetype="application/zip",
        as_attachment=True,
        download_name="motiffinder_output.zip",
    )

from flask import Blueprint, request, jsonify, send_file
import os, tempfile
from .codon_bias_estimator import read_input_seq, codon_usage  # adjust import path

web = Blueprint('web', __name__)

@web.route('/codonbias/run', methods=['POST'])
def run_codonbias():
    f = request.files.get('genome_file')
    codon_table = request.form.get('codon_table', '11')
    if not f:
        return jsonify({'success': False, 'error': 'No file provided'}), 400

    with tempfile.NamedTemporaryFile(suffix='.gbk', delete=False) as tmp:
        f.save(tmp.name)
        tmp_path = tmp.name

    try:
        genome = read_input_seq(tmp_path)
        if not genome['correct']:
            return jsonify({'success': False, 'error': '; '.join(genome['messages'])})

        result = codon_usage(genome['total_cds_concats'], codon_table)
        cds_count = sum(len(v['cds']) for v in genome['annotations'].values())

        out_dir = tempfile.mkdtemp()
        csv_path = os.path.join(out_dir, 'codon_usage_table.csv')
        result['Table'].to_csv(csv_path, index=False)

        # store path in session or return download URL
        return jsonify({'success': True, 'cds_count': cds_count, 'csv_path': csv_path})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})
    finally:
        os.unlink(tmp_path)