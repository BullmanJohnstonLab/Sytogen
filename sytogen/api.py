from flask import Blueprint, request, send_file, abort, jsonify
from werkzeug.exceptions import HTTPException
from Bio.SeqFeature import FeatureLocation
import io
import os
import uuid
import tempfile
import zipfile
from werkzeug.utils import secure_filename
from Bio import SeqIO
from threading import Thread

from sytogen.scripts.motiffinder_backend import (
    parse_rebase_motifs,
    search_motifs,
    hits_to_gff3,
    hits_to_tsv,
    load_gff3_features,
    GFF3_HEADER,
    TSV_HEADER,
)

from sytogen.scripts.optimizer import run_step1
from sytogen.scripts.heuristic import run_step2

from sytogen.scripts.codon_bias_estimator import run_codon_bias
# ------------------------------------------------------------------
# Blueprint
# ------------------------------------------------------------------

api = Blueprint("api", __name__)

JOBS = {}

# ------------------------------------------------------------------
# Global error handler (JSON errors everywhere)
# ------------------------------------------------------------------

@api.app_errorhandler(HTTPException)
def handle_http_exception(e):
    return jsonify(error=e.description), e.code


# ------------------------------------------------------------------
# MotifFinder endpoint
# ------------------------------------------------------------------

@api.route("/motiffinder/run", methods=["POST"])
def run_motiffinder_sync():
    """
    Run MotifFinder synchronously and return a ZIP of results.
    """

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

    motif_text = motif_file.read().decode("utf-8", errors="replace")
    motifs = parse_rebase_motifs(motif_text)
    if not motifs:
        abort(400, "No valid motifs found in motif_file")

    seq_text = seq_file.read().decode("utf-8", errors="replace")

    try:
        if source_type == "genbank":
            records = list(SeqIO.parse(io.StringIO(seq_text), "genbank"))
        else:
            records = list(SeqIO.parse(io.StringIO(seq_text), "fasta"))
    except Exception as e:
        abort(400, f"Failed to parse sequence file: {e}")

    if not records:
        abort(400, "No sequences found in sequence_file")

    features = []
    if source_type == "genbank":
        for rec in records:
            for feat in rec.features:
                loc = feat.location
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
    else:
        if ann_file:
            gff_text = ann_file.read().decode("utf-8", errors="replace")
            features = load_gff3_features(gff_text)

    all_gff3_parts = [GFF3_HEADER]
    all_tsv_parts  = [TSV_HEADER]

    for rec in records:
        seqid = rec.id or rec.name or "unknown"
        seq_str = str(rec.seq)
        seq_len = len(seq_str)

        is_circular = (
            source_type == "genbank"
            and rec.annotations.get("topology", "").lower() == "circular"
        )

        hits = search_motifs(seq_str, motifs, is_circular=is_circular)
        rec_features = [f for f in features if f["seqid"] == seqid]

        gff3_body = hits_to_gff3(hits, seqid, seq_len, rec_features)
        body_lines = [
            l for l in gff3_body.splitlines(keepends=True)
            if not l.startswith("#")
        ]

        all_gff3_parts.append(f"##sequence-region {seqid} 1 {seq_len}\n")
        all_gff3_parts.extend(body_lines)

        tsv_body = hits_to_tsv(hits, seqid, rec_features)
        all_tsv_parts.extend(tsv_body.splitlines(keepends=True)[1:])

    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("motiffinder_results.gff3", "".join(all_gff3_parts))
        zf.writestr("motiffinder_summary.tsv", "".join(all_tsv_parts))

    zip_buf.seek(0)

    return send_file(
        zip_buf,
        mimetype="application/zip",
        as_attachment=True,
        download_name="motiffinder_output.zip",
    )

# ------------------------------------------------------------------
# CodonBias endpoint
# ------------------------------------------------------------------

@api.route("/codonbias/run", methods=["POST"])
def run_codonbias():

    # ----------------------------
    # Parse inputs
    # ----------------------------
    try:
        codon_table = int(request.form.get("codon_table", 11))
    except ValueError:
        return jsonify(error="Invalid codon table"), 400

    source_type = request.form.get("source_type")

    if source_type not in {"genbank", "fasta"}:
        return jsonify(error="Invalid source_type"), 400

    # ----------------------------
    # Temp workspace
    # ----------------------------
    with tempfile.TemporaryDirectory() as tmpdir:

        try:
            # ============================
            # FASTA + GFF MODE
            # ============================
            if source_type == "fasta":

                fasta = request.files.get("fasta_file")
                gff   = request.files.get("gff_file")

                if not fasta or not gff:
                    return jsonify(error="FASTA + GFF required"), 400

                if not fasta.filename or not gff.filename:
                    return jsonify(error="Missing filenames"), 400

                fasta_path = os.path.join(tmpdir, secure_filename(fasta.filename))
                gff_path   = os.path.join(tmpdir, secure_filename(gff.filename))

                fasta.save(fasta_path)
                gff.save(gff_path)

                csv_path = run_codon_bias(
                    fasta_path=fasta_path,
                    gff_path=gff_path,
                    codon_table=codon_table,
                    output_dir=tmpdir,
                )

            # ============================
            # GENBANK MODE
            # ============================
            else:
                gbk = request.files.get("genome_file")

                if not gbk or not gbk.filename:
                    return jsonify(error="GenBank file required"), 400

                gbk_path = os.path.join(tmpdir, secure_filename(gbk.filename))
                gbk.save(gbk_path)

                csv_path = run_codon_bias(
                    genome_path=gbk_path,
                    codon_table=codon_table,
                    output_dir=tmpdir,
                )

            # ----------------------------
            # Create ZIP
            # ----------------------------
            zip_path = os.path.join(tmpdir, "codonbias_output.zip")

            with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
                z.write(csv_path, arcname="codon_usage_table.csv")

            # ----------------------------
            # Return file
            # ----------------------------
            return send_file(
                zip_path,
                mimetype="application/zip",
                as_attachment=True,
                download_name="codonbias_output.zip",
            )

        except Exception as e:
            # critical: return JSON, not HTML
            return jsonify(error=str(e)), 500
        
def worker(job_id, paths, params, tmpdir):

    try:
        JOBS[job_id]["status"] = "running"
        params = {**params, "output_dir": tmpdir}

        step1 = run_step1(paths, params)
        result = run_step2(step1, params)

        result_path = os.path.join(tmpdir, "result.fasta")
        with open(result_path, "w") as f:
            f.write(result)

        JOBS[job_id]["status"] = "done"
        JOBS[job_id]["result"] = result_path

    except Exception as e:
        JOBS[job_id]["status"] = "error"
        JOBS[job_id]["error"] = str(e)
    
@api.route("/status/<job_id>", methods=["GET"])
def status(job_id):

    job = JOBS.get(job_id)

    if not job:
        return jsonify({"status": "unknown"}), 404

    return jsonify({
        "status": job["status"],
        "error": job.get("error")
    })
    
@api.route("/sytogen/run", methods=["POST"])
def run_sytogen():

    job_id = str(uuid.uuid4())
    tmpdir = tempfile.mkdtemp()

    paths = {}

    for key in [
        "genbank",
        "fasta",
        "gff",
        "codon_usage",
        "motif_table",
        "genome",
        "strain_genome",
    ]:
        f = request.files.get(key)
        if f:
            path = os.path.join(tmpdir, secure_filename(f.filename))
            f.save(path)
            paths[key] = path

    params = {
        "avoid_regions": request.form.get("avoid_regions"),
        "flex_regions": request.form.get("flex_regions"),
        "forced_edits": request.form.get("forced_edits"),
        "excluded_edits": request.form.get("excluded_edits"),
        "preserve_gc": request.form.get("preserve_gc") == "true",
        "avoid_new_motifs": request.form.get("avoid_new_motifs") == "true",
        "strict_synonymous": request.form.get("strict_synonymous") == "true",
        "codon_table": request.form.get("codon_table", 11),
    }

    JOBS[job_id] = {
    "status": "queued",   # queued | running | done | error
    "result": None,
    "error": None,
    "tmpdir": tmpdir
}

    Thread(target=worker, args=(job_id, paths, params, tmpdir)).start()

    return jsonify({"job_id": job_id})

@api.route("/sytogen/result/<job_id>", methods=["GET"])
def result(job_id):

    job = JOBS.get(job_id)

    if not job:
        return jsonify({"error": "invalid job"}), 404

    if job["status"] != "done":
        return jsonify({"error": "not ready"}), 202

    if not job.get("result") or not os.path.exists(job["result"]):
        return jsonify({"error": "result missing"}), 500

    return send_file(
        job["result"],
        mimetype="text/plain",
        as_attachment=True,
        download_name="sytogen_result.fasta"
    )
