import io
import json
import os
import uuid
import shutil
import zipfile
import tempfile
import traceback
import pandas as pd
import copy

from threading import Thread

from flask import (
    Blueprint,
    request,
    send_file,
    abort,
    jsonify,
)

from werkzeug.datastructures import FileStorage
from werkzeug.exceptions import HTTPException
from werkzeug.utils import secure_filename

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (
    FeatureLocation,
    SeqFeature,
)

from sytogen.scripts.motiffinder_backend import (
    parse_rebase_motifs,
    search_motifs,
    hits_to_gff3,
    hits_to_tsv,
    load_gff3_features,
    GFF3_HEADER,
    TSV_HEADER,
)

from sytogen.scripts.codon_bias_estimator import (
    run_codon_bias,
)
from sytogen.scripts.sytogen_runner import (
    run_sytogen_pipeline,
    decision_matrix_to_tsv,
)

# =========================================================
# Blueprint
# =========================================================

api = Blueprint("api", __name__)

JOBS = {}

# =========================================================
# Allowed extensions
# =========================================================

GENBANK_EXTENSIONS = {
    ".gb",
    ".gbk",
    ".gbff",
    ".genbank",
    ".gbf",
}

FASTA_EXTENSIONS = {
    ".fa",
    ".fasta",
    ".fna",
    ".ffn",
    ".faa",
}

GFF_EXTENSIONS = {
    ".gff",
    ".gff3",
}


# =========================================================
# Helpers
# =========================================================

def allowed_extension(filename, allowed):

    ext = os.path.splitext(filename)[1].lower()

    return ext in allowed


def read_uploaded_table(file_storage):
    text = file_storage.stream.read().decode("utf-8-sig")
    return pd.read_csv(
        io.StringIO(text),
        sep=None,
        engine="python",
    )


# =========================================================
# Global JSON error handler
# =========================================================

@api.app_errorhandler(HTTPException)
def handle_http_exception(e):

    return jsonify(
        error=e.description
    ), e.code


# =========================================================
# MotifFinder endpoint
# =========================================================

@api.route("/motiffinder/run", methods=["POST"])
def run_motiffinder_sync():

    missing = []

    if "sequence_file" not in request.files:
        missing.append("sequence_file")

    if "motif_file" not in request.files:
        missing.append("motif_file")

    if missing:
        abort(
            400,
            f"Missing required files: "
            f"{', '.join(missing)}"
        )

    source_type = request.form.get(
        "source_type",
        "",
    ).lower()

    if source_type not in {
        "genbank",
        "fasta",
    }:
        abort(
            400,
            "source_type must be "
            "'genbank' or 'fasta'"
        )

    seq_file = request.files["sequence_file"]
    motif_file = request.files["motif_file"]
    ann_file = request.files.get(
        "annotation_file"
    )

    motif_text = motif_file.read().decode(
        "utf-8",
        errors="replace",
    )

    motifs = parse_rebase_motifs(
        motif_text
    )

    if not motifs:
        abort(
            400,
            "No valid motifs found"
        )

    seq_text = seq_file.read().decode(
        "utf-8",
        errors="replace",
    )

    try:

        if source_type == "genbank":

            records = list(
                SeqIO.parse(
                    io.StringIO(seq_text),
                    "genbank",
                )
            )

        else:

            records = list(
                SeqIO.parse(
                    io.StringIO(seq_text),
                    "fasta",
                )
            )

    except Exception as e:

        abort(
            400,
            f"Failed to parse sequence "
            f"file: {e}"
        )

    if not records:
        abort(
            400,
            "No sequences found"
        )

    features = []

    if source_type == "genbank":

        for rec in records:

            for feat in rec.features:

                loc = feat.location

                seqid = rec.id or rec.name

                attrs = (
                    f"ID={feat.type}_"
                    f"{int(loc.start)+1}_"
                    f"{int(loc.end)};"
                    f"Name="
                    f"{feat.qualifiers.get('gene', [feat.type])[0]}"
                )

                features.append({
                    "seqid": seqid,
                    "source": "GenBank",
                    "type": feat.type,
                    "start": int(loc.start) + 1,
                    "end": int(loc.end),
                    "score": ".",
                    "strand":
                        "+"
                        if loc.strand == 1
                        else "-"
                        if loc.strand == -1
                        else ".",
                    "phase": ".",
                    "attrs": attrs,
                })

    else:

        if ann_file:

            gff_text = ann_file.read().decode(
                "utf-8",
                errors="replace",
            )

            features = load_gff3_features(
                gff_text
            )

    all_gff3_parts = [GFF3_HEADER]
    all_tsv_parts = [TSV_HEADER]

    for rec in records:

        rec.annotations.setdefault(
            "molecule_type",
            "DNA",
        )

        seqid = (
            rec.id
            or rec.name
            or "unknown"
        )

        seq_str = str(rec.seq)

        seq_len = len(seq_str)

        is_circular = (
            source_type == "genbank"
            and rec.annotations.get(
                "topology",
                "",
            ).lower() == "circular"
        )

        hits = search_motifs(
            seq_str,
            motifs,
            is_circular=is_circular,
        )

        rec_features = [
            f for f in features
            if f["seqid"] == seqid
        ]

        for i, hit in enumerate(
            hits,
            start=1,
        ):

            start = hit["pos_0"]

            end = (
                start
                + len(hit["rec_seq"])
            )

            if end > seq_len:
                continue

            rec.features.append(
                SeqFeature(
                    FeatureLocation(
                        start,
                        end,
                        strand=(
                            1
                            if hit["strand"] == "+"
                            else -1
                        ),
                    ),
                    type="misc_feature",
                    qualifiers={
                        "ID": [
                            f"motif_hit_{i:04d}"
                        ],
                        "note": [
                            f"MotifFinder hit "
                            f"{hit['rec_seq']}"
                        ],
                    },
                )
            )

        gff3_body = hits_to_gff3(
            hits,
            seqid,
            seq_len,
            rec_features,
        )

        body_lines = [
            l
            for l in gff3_body.splitlines(
                keepends=True
            )
            if not l.startswith("#")
        ]

        all_gff3_parts.append(
            f"##sequence-region "
            f"{seqid} 1 {seq_len}\n"
        )

        all_gff3_parts.extend(
            body_lines
        )

        tsv_body = hits_to_tsv(
            hits,
            seqid,
            seq_len,
            rec_features,
        )

        all_tsv_parts.extend(
            tsv_body.splitlines(
                keepends=True
            )[1:]
        )

    zip_buf = io.BytesIO()

    gbk_buf = io.StringIO()

    SeqIO.write(
        records,
        gbk_buf,
        "genbank",
    )

    with zipfile.ZipFile(
        zip_buf,
        "w",
        zipfile.ZIP_DEFLATED,
    ) as zf:

        zf.writestr(
            "motiffinder_annotated.gbk",
            gbk_buf.getvalue(),
        )

        zf.writestr(
            "motiffinder_results.gff3",
            "".join(all_gff3_parts),
        )

        zf.writestr(
            "motiffinder_summary.tsv",
            "".join(all_tsv_parts),
        )

    zip_buf.seek(0)

    return send_file(
        zip_buf,
        mimetype="application/zip",
        as_attachment=True,
        download_name="motiffinder_output.zip",
    )


# =========================================================
# CodonBias endpoint
# =========================================================

@api.route("/codonbias/run", methods=["POST"])
def run_codonbias():

    try:

        codon_table = int(
            request.form.get(
                "codon_table",
                11,
            )
        )

    except ValueError:

        return jsonify(
            error="Invalid codon table"
        ), 400

    source_type = request.form.get(
        "source_type"
    )

    if source_type not in {
        "genbank",
        "fasta",
    }:
        return jsonify(
            error="Invalid source_type"
        ), 400

    tmpdir = tempfile.mkdtemp(
        prefix="codonbias_"
    )

    try:

        # -------------------------------------------------
        # FASTA + GFF MODE
        # -------------------------------------------------

        if source_type == "fasta":

            fasta = request.files.get(
                "fasta_file"
            )

            gff = request.files.get(
                "gff_file"
            )

            if not fasta or not gff:
                return jsonify(
                    error="FASTA + GFF required"
                ), 400

            fasta_name = secure_filename(
                fasta.filename
            )

            gff_name = secure_filename(
                gff.filename
            )

            if not allowed_extension(
                fasta_name,
                FASTA_EXTENSIONS,
            ):
                return jsonify(
                    error="Invalid FASTA extension"
                ), 400

            if not allowed_extension(
                gff_name,
                GFF_EXTENSIONS,
            ):
                return jsonify(
                    error="Invalid GFF extension"
                ), 400

            fasta_path = os.path.join(
                tmpdir,
                fasta_name,
            )

            gff_path = os.path.join(
                tmpdir,
                gff_name,
            )

            fasta.save(fasta_path)
            gff.save(gff_path)

            output_paths = run_codon_bias(
                fasta_path=fasta_path,
                gff_path=gff_path,
                codon_table=codon_table,
                output_dir=tmpdir,
            )

        # -------------------------------------------------
        # GENBANK MODE
        # -------------------------------------------------

        else:

            gbk = request.files.get(
                "genome_file"
            )

            if not gbk:
                return jsonify(
                    error="GenBank file required"
                ), 400

            gbk_name = secure_filename(
                gbk.filename
            )

            if not allowed_extension(
                gbk_name,
                GENBANK_EXTENSIONS,
            ):
                return jsonify(
                    error="Invalid GenBank extension"
                ), 400

            gbk_path = os.path.join(
                tmpdir,
                gbk_name,
            )

            gbk.save(gbk_path)

            output_paths = run_codon_bias(
                genome_path=gbk_path,
                codon_table=codon_table,
                output_dir=tmpdir,
            )

        zip_path = os.path.join(
            tmpdir,
            "codonbias_output.zip",
        )

        with zipfile.ZipFile(
            zip_path,
            "w",
            compression=zipfile.ZIP_DEFLATED,
        ) as z:

            z.write(
                output_paths["csv"],
                arcname="codon_usage_table.csv",
            )

            z.write(
                output_paths["genbank"],
                arcname="codonbias_input.gbk",
            )

            z.write(
                output_paths["fasta"],
                arcname="codonbias_input.fasta",
            )

            z.write(
                output_paths["gff"],
                arcname="codonbias_input.gff3",
            )

        return send_file(
            zip_path,
            mimetype="application/zip",
            as_attachment=True,
            download_name="codonbias_output.zip",
        )

    except ValueError as e:

        return jsonify(
            error=str(e)
        ), 400

    except Exception as e:

        traceback.print_exc()

        return jsonify(
            error=str(e),
            traceback=traceback.format_exc(),
        ), 500


# =========================================================
# Worker
# =========================================================

def worker(job_id, paths, params, tmpdir):

    try:
        from sytogen.scripts.optimizer import run_step1
        from sytogen.scripts.heuristic import run_step2

        JOBS[job_id]["status"] = "running"

        params = {
            **params,
            "output_dir": tmpdir,
        }

        step1 = run_step1(
            paths,
            params,
        )

        result = run_step2(
            step1,
            params,
        )

        result_path = os.path.join(
            tmpdir,
            "result.fasta",
        )

        with open(result_path, "w") as f:
            f.write(result)

        JOBS[job_id]["status"] = "done"
        JOBS[job_id]["result"] = result_path

    except Exception as e:

        traceback.print_exc()

        JOBS[job_id]["status"] = "error"
        JOBS[job_id]["error"] = str(e)
        JOBS[job_id]["traceback"] = (
            traceback.format_exc()
        )


# =========================================================
# Status endpoint
# =========================================================

@api.route("/status/<job_id>", methods=["GET"])
def status(job_id):

    job = JOBS.get(job_id)

    if not job:

        return jsonify({
            "status": "unknown"
        }), 404

    return jsonify({
        "status": job["status"],
        "error": job.get("error"),
    })


# =========================================================
# Run Sytogen
# =========================================================

@api.route("/sytogen/run", methods=["POST"])
def run_sytogen():
    gbk_file = request.files.get("genbank")
    codon_file = request.files.get("codon_usage")
    motif_file = request.files.get("motif_table")

    if not gbk_file or not codon_file or not motif_file:
        return jsonify(error="Missing uploaded files"), 400

    topology = request.form.get("topology", "circular").lower()
    if topology not in {"circular", "linear"}:
        return jsonify(error="topology must be 'circular' or 'linear'"), 400

    try:
        # =================================================
        # PARSE OBJECTS
        # =================================================

        # Convert GenBank → SeqRecord
        seq_record = SeqIO.read(io.TextIOWrapper(gbk_file.stream), "genbank")

        # Convert uploaded tables to DataFrames, accepting CSV or TSV output.
        codon_df = read_uploaded_table(codon_file)
        motif_df = read_uploaded_table(motif_file)

        # =================================================
        # RUN PIPELINE
        # =================================================

        result = run_sytogen_pipeline(
            seq_record,
            codon_df,
            motif_df,
            params={
                "topology":    topology,
                "preserve_gc": request.form.get("preserve_gc") == "true",
            }
        )

        # =================================================
        # RETURN RESULT
        # =================================================


        zip_buffer = io.BytesIO()
        output_record = copy.deepcopy(seq_record)
        output_record.seq = Seq(result["altered_sequence"])
        output_record.id = f"{seq_record.id}_sytogen"
        output_record.name = f"{seq_record.name}_sytogen"
        output_record.description = f"{seq_record.description} | SyToGen result"

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("sytogen_result.fasta",    result["altered_fasta"])
            zf.writestr("sytogen_result.gbk",      output_record.format("genbank"))
            zf.writestr("original_sequence.fasta", result["original_fasta"])
            zf.writestr("input_sequence.gbk",      seq_record.format("genbank"))
            zf.writestr(
                "decision_matrix.tsv",
                decision_matrix_to_tsv(result["decision_matrix"]),
            )
            zf.writestr(
                "summary.json",
                json.dumps(result["summary"], indent=2),
            )

        zip_buffer.seek(0)

        return send_file(
            zip_buffer,
            mimetype="application/zip",
            as_attachment=True,
            download_name="sytogen_output.zip",
        )
    except Exception as e:
        traceback.print_exc()
        return jsonify(error=str(e)), 500


# =========================================================
# Sytogen async submit endpoint (was entirely missing)
# =========================================================

@api.route("/sytogen/submit", methods=["POST"])
def submit_sytogen():
    gbk_file = request.files.get("genbank")
    codon_file = request.files.get("codon_usage")
    motif_file = request.files.get("motif_table")

    if not gbk_file or not codon_file or not motif_file:
        return jsonify(error="Missing uploaded files"), 400

    topology = request.form.get("topology", "circular").lower()
    if topology not in {"circular", "linear"}:
        return jsonify(error="topology must be 'circular' or 'linear'"), 400

    tmpdir = tempfile.mkdtemp(prefix="sytogen_")

    # Save uploads to disk so the worker thread can read them
    gbk_path   = os.path.join(tmpdir, secure_filename(gbk_file.filename))
    codon_path = os.path.join(tmpdir, secure_filename(codon_file.filename))
    motif_path = os.path.join(tmpdir, secure_filename(motif_file.filename))
    gbk_file.save(gbk_path)
    codon_file.save(codon_path)
    motif_file.save(motif_path)

    job_id = str(uuid.uuid4())
    JOBS[job_id] = {"status": "queued"}

    paths = {
        "genbank":     gbk_path,
        "codon_usage": codon_path,
        "motif_table": motif_path,
    }
    params = {
        "topology":    topology,
        "preserve_gc": request.form.get("preserve_gc") == "true",
    }

    Thread(
        target=worker,
        args=(job_id, paths, params, tmpdir),
        daemon=True,
    ).start()

    return jsonify(job_id=job_id), 202


# =========================================================
# Result endpoint
# =========================================================

@api.route(
    "/sytogen/result/<job_id>",
    methods=["GET"],
)
def result(job_id):

    job = JOBS.get(job_id)

    if not job:

        return jsonify({
            "error": "invalid job"
        }), 404

    if job["status"] != "done":

        return jsonify({
            "error": "not ready"
        }), 202

    result_path = job.get("result")

    if (
        not result_path
        or not os.path.exists(result_path)):

        return jsonify({
            "error": "result missing"
        }), 500

    return send_file(
        result_path,
        mimetype="text/plain",
        as_attachment=True,
        download_name="sytogen_result.fasta",
    )
