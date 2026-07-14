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
    assembly_plan_to_tsv,
    assembly_plan_fragments_fasta,
    assembly_plan_summary,
    assembly_primers_to_tsv,
    strip_backbone,
)
from sytogen.scripts.visualization import render_plasmid_maps

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


def read_backbone_record(file_storage):
    """
    Parse an optional vector-backbone FASTA upload. Must contain exactly
    one record — same validation legacy_sytogen.sequence_preprocess()
    applied — so it's unambiguous which sequence to locate and strip out
    of the construct.
    """
    text = file_storage.stream.read().decode("utf-8-sig")
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    if len(records) != 1:
        raise ValueError(
            f"Backbone file must contain exactly one FASTA record, found {len(records)}."
        )
    return records[0]


def _parse_gff3_attrs(attrs):
    """
    'ID=CDS_1_300;Name=geneA;locus_tag=geneA_1' -> {'ID': ['CDS_1_300'],
    'Name': ['geneA'], 'locus_tag': ['geneA_1']}. GFF3's attribute column
    is a ';'-separated list of key=value pairs; this is the same format
    load_gff3_features()'s own callers already assume (see the manual
    'ID=...;Name=...' construction in run_motiffinder_sync above).
    """
    qualifiers = {}
    for part in (attrs or "").split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        qualifiers.setdefault(key.strip(), []).append(value.strip())
    return qualifiers


def build_record_from_fasta_gff(fasta_text, gff_text, topology="circular"):
    """
    Build a Bio.SeqRecord (with real Bio.SeqFeature features, not GFF3
    dicts) from a FASTA sequence + GFF3 annotation pair, so it's a drop-in
    replacement for SeqIO.read(..., 'genbank') everywhere downstream —
    run_sytogen_pipeline / _parse_genes / _parse_protected_regions never
    need to know which input format the person actually uploaded.

    Reuses load_gff3_features() (already used by the MotifFinder/CodonBias
    FASTA+GFF3 path) for the actual GFF3 parsing rather than writing a
    second parser.
    """
    records = list(SeqIO.parse(io.StringIO(fasta_text), "fasta"))
    if len(records) != 1:
        raise ValueError(
            f"FASTA file must contain exactly one sequence, found {len(records)}."
        )
    record = records[0]
    seqid = record.id or record.name

    raw_features = load_gff3_features(gff_text)
    # A GFF3 can carry annotations for multiple contigs/sequences; only
    # the ones matching this FASTA record's id are relevant. If nothing
    # matches by id but the GFF only names one seqid anyway (a common
    # mismatch — e.g. FASTA header has extra description text GFF3
    # export truncated), fall back to using all of it rather than
    # silently producing zero genes.
    matching = [f for f in raw_features if f.get("seqid") == seqid]
    if not matching and raw_features:
        distinct_seqids = {f.get("seqid") for f in raw_features}
        if len(distinct_seqids) == 1:
            matching = raw_features
        else:
            raise ValueError(
                f"No GFF3 features matched FASTA sequence id '{seqid}' "
                f"(GFF3 contains: {sorted(distinct_seqids)})."
            )

    features = []
    for f in matching:
        strand_symbol = f.get("strand", ".")
        strand = 1 if strand_symbol == "+" else -1 if strand_symbol == "-" else 1
        start = int(f["start"]) - 1  # GFF3 is 1-based inclusive -> BioPython 0-based
        end = int(f["end"])
        features.append(SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type=f.get("type", "misc_feature"),
            qualifiers=_parse_gff3_attrs(f.get("attrs", "")),
        ))

    record.features = features
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = topology
    return record


def read_uploaded_table(file_storage):
    text = file_storage.stream.read().decode("utf-8-sig")
    return pd.read_csv(
        io.StringIO(text),
        sep=None,
        engine="python",
    )


def parse_motif_text(text):
    """
    Core motif-table parsing logic, operating on raw text so it can be
    shared by both the synchronous upload path (read_motif_table, below)
    and the async worker path (which reads the same file back from disk).

    Returns a DataFrame with a 'motif' column, as sytogen_runner._parse_motifs()
    expects. Accepts either a plain delimited table with a 'motif' column,
    or a REBASE-style tagged export (e.g. "<enz_type>2<rec_seq>ATGC...<>"),
    falling back to the existing parse_rebase_motifs() parser for the latter.
    """
    # --- Attempt 1: plain delimited table with a 'motif' column ---
    try:
        df = pd.read_csv(io.StringIO(text), sep=None, engine="python")
        normalized_cols = {str(c).strip().lower() for c in df.columns}
        if "motif" in normalized_cols:
            return df
    except Exception:
        pass  # not a plain delimited table — fall through to REBASE parsing

    # --- Attempt 2: REBASE-style tagged export ---
    motifs = parse_rebase_motifs(text)
    if not motifs:
        raise ValueError(
            "Could not parse the restriction motif table. Expected either "
            "a delimited file with a 'motif' column, or a REBASE-style "
            "tagged export (e.g. containing '<rec_seq>...' entries)."
        )

    motif_df = pd.DataFrame(motifs)

    # parse_rebase_motifs() returns REBASE field names (e.g. 'rec_seq');
    # normalize whichever recognition-sequence field it used to 'motif'.
    if "motif" not in motif_df.columns:
        for candidate in ("rec_seq", "recognition_sequence", "sequence", "seq"):
            if candidate in motif_df.columns:
                motif_df = motif_df.rename(columns={candidate: "motif"})
                break

    if "motif" not in motif_df.columns:
        raise ValueError(
            "Could not find a recognition-sequence field in the motif table."
        )

    return motif_df


def read_motif_table(file_storage):
    """Parse an uploaded motif-table file (see parse_motif_text)."""
    text = file_storage.stream.read().decode("utf-8-sig")
    return parse_motif_text(text)


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
        JOBS[job_id]["status"] = "running"

        # Parse inputs the same way the synchronous /sytogen/run endpoint
        # does, reusing the same helpers (including the REBASE-format
        # motif-table fallback) so both code paths behave identically.
        source_type = params.get("source_type", "genbank")
        if source_type == "genbank":
            seq_record = SeqIO.read(paths["genbank"], "genbank")
        else:
            with open(paths["fasta_file"], "r", encoding="utf-8-sig") as f:
                fasta_text = f.read()
            with open(paths["gff_file"], "r", encoding="utf-8-sig") as f:
                gff_text = f.read()
            seq_record = build_record_from_fasta_gff(
                fasta_text, gff_text, params.get("topology", "circular")
            )

        with open(paths["codon_usage"], "r", encoding="utf-8-sig") as f:
            codon_df = pd.read_csv(io.StringIO(f.read()), sep=None, engine="python")

        with open(paths["motif_table"], "r", encoding="utf-8-sig") as f:
            motif_df = parse_motif_text(f.read())

        # Optional vector-backbone removal — same as the sync endpoint.
        if paths.get("backbone"):
            with open(paths["backbone"], "r", encoding="utf-8-sig") as f:
                backbone_records = list(SeqIO.parse(io.StringIO(f.read()), "fasta"))
            if len(backbone_records) != 1:
                raise ValueError(
                    f"Backbone file must contain exactly one FASTA record, "
                    f"found {len(backbone_records)}."
                )
            seq_record = strip_backbone(
                seq_record, backbone_records[0], params.get("topology", "circular")
            )

        result = run_sytogen_pipeline(
            seq_record,
            codon_df,
            motif_df,
            params={
                "topology":              params.get("topology", "circular"),
                "preserve_gc":           params.get("preserve_gc", False),
                "include_assembly_plan": params.get("include_assembly_plan", False),
            },
        )

        # Build the same output bundle the sync endpoint returns, so async
        # jobs also get the full decision matrix, summary, and mutated
        # GenBank/FASTA — not just a bare sequence.
        output_record = copy.deepcopy(seq_record)
        output_record.seq = Seq(result["altered_sequence"])
        output_record.id = f"{seq_record.id}_sytogen"
        output_record.name = f"{seq_record.name}_sytogen"
        output_record.description = f"{seq_record.description} | SyToGen result"
        for mutation in result["applied_mutations"]:
            output_record.features.append(
                SeqFeature(
                    FeatureLocation(
                        mutation.position,
                        mutation.position + len(mutation.new),
                    ),
                    type="SyT",
                    qualifiers={
                        "label": [f"{mutation.old} --> {mutation.new}"],
                    },
                )
            )
        motifs_used = motif_df.to_csv(sep="\t", index=False)

        plasmid_maps = render_plasmid_maps(
            output_record,
            result["motifs"],
            title=seq_record.id,
        )

        zip_path = os.path.join(tmpdir, "sytogen_output.zip")
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("sytogen_result.fasta",    result["altered_fasta"])
            zf.writestr("sytogen_result.gbk",      output_record.format("genbank"))
            zf.writestr("original_sequence.fasta", result["original_fasta"])
            zf.writestr("input_sequence.gbk",      seq_record.format("genbank"))
            zf.writestr("motifs_used.tsv",         motifs_used)
            zf.writestr(
                "decision_matrix.tsv",
                decision_matrix_to_tsv(result["decision_matrix"]),
            )
            for filename, contents in plasmid_maps.items():
                zf.writestr(filename, contents)
            zf.writestr(
                "summary.json",
                json.dumps(result["summary"], indent=2),
            )
            zf.writestr(
                "new_motifs_check.json",
                json.dumps({
                    "new_motifs_introduced": result["summary"]["new_motifs_introduced"],
                    "new_motifs": result["new_motifs"],
                }, indent=2),
            )
            if result.get("assembly_plan"):
                zf.writestr(
                    "assembly_plan.tsv",
                    assembly_plan_to_tsv(result["assembly_plan"]),
                )
                zf.writestr(
                    "assembly_fragments.fasta",
                    assembly_plan_fragments_fasta(result["assembly_plan"]),
                )
                zf.writestr(
                    "assembly_primers.tsv",
                    assembly_primers_to_tsv(result["assembly_plan"]),
                )
                zf.writestr(
                    "assembly_plan_summary.json",
                    json.dumps(assembly_plan_summary(result["assembly_plan"]), indent=2),
                )

        JOBS[job_id]["status"] = "done"
        JOBS[job_id]["result"] = zip_path

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
    source_type = request.form.get("source_type", "genbank").lower()
    if source_type not in {"genbank", "fasta"}:
        return jsonify(error="source_type must be 'genbank' or 'fasta'"), 400

    gbk_file    = request.files.get("genbank")
    fasta_file  = request.files.get("fasta_file")
    gff_file    = request.files.get("gff_file")
    codon_file  = request.files.get("codon_usage")
    motif_file  = request.files.get("motif_table")

    if source_type == "genbank":
        if not gbk_file:
            return jsonify(error="Missing uploaded GenBank file"), 400
    else:
        if not fasta_file or not gff_file:
            return jsonify(error="FASTA + GFF3 mode requires both files"), 400

    if not codon_file or not motif_file:
        return jsonify(error="Missing uploaded files"), 400

    topology = request.form.get("topology", "circular").lower()
    if topology not in {"circular", "linear"}:
        return jsonify(error="topology must be 'circular' or 'linear'"), 400

    try:
        # =================================================
        # PARSE OBJECTS
        # =================================================

        if source_type == "genbank":
            seq_record = SeqIO.read(io.TextIOWrapper(gbk_file.stream), "genbank")
        else:
            fasta_text = fasta_file.stream.read().decode("utf-8-sig")
            gff_text   = gff_file.stream.read().decode("utf-8-sig")
            seq_record = build_record_from_fasta_gff(fasta_text, gff_text, topology)

        # Convert uploaded tables to DataFrames, accepting CSV or TSV output,
        # and REBASE-tagged exports for the motif table.
        codon_df = read_uploaded_table(codon_file)
        motif_df = read_motif_table(motif_file)

        # Optional vector-backbone removal — if provided, locate and strip
        # the backbone before anything else runs, so RM-silencing and
        # assembly planning only ever see the insert.
        backbone_file = request.files.get("backbone")
        if backbone_file and backbone_file.filename:
            backbone_record = read_backbone_record(backbone_file)
            seq_record = strip_backbone(seq_record, backbone_record, topology)

        # =================================================
        # RUN PIPELINE
        # =================================================

        result = run_sytogen_pipeline(
            seq_record,
            codon_df,
            motif_df,
            params={
                "topology":              topology,
                "preserve_gc":           request.form.get("preserve_gc") == "true",
                "include_assembly_plan": request.form.get("include_assembly_plan") == "true",
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
        for mutation in result["applied_mutations"]:
            output_record.features.append(
                SeqFeature(
                    FeatureLocation(
                        mutation.position,
                        mutation.position + len(mutation.new),
                    ),
                    type="SyT",
                    qualifiers={
                        "label": [f"{mutation.old} --> {mutation.new}"],
                    },
                )
            )
        motifs_used = motif_df.to_csv(sep="\t", index=False)

        plasmid_maps = render_plasmid_maps(
            output_record,
            result["motifs"],
            title=seq_record.id,
        )

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("sytogen_result.fasta",    result["altered_fasta"])
            zf.writestr("sytogen_result.gbk",      output_record.format("genbank"))
            zf.writestr("original_sequence.fasta", result["original_fasta"])
            zf.writestr("input_sequence.gbk",      seq_record.format("genbank"))
            zf.writestr("motifs_used.tsv",         motifs_used)
            zf.writestr(
                "decision_matrix.tsv",
                decision_matrix_to_tsv(result["decision_matrix"]),
            )
            for filename, contents in plasmid_maps.items():
                zf.writestr(filename, contents)
            zf.writestr(
                "summary.json",
                json.dumps(result["summary"], indent=2),
            )
            zf.writestr(
                "new_motifs_check.json",
                json.dumps({
                    "new_motifs_introduced": result["summary"]["new_motifs_introduced"],
                    "new_motifs": result["new_motifs"],
                }, indent=2),
            )
            if result.get("assembly_plan"):
                zf.writestr(
                    "assembly_plan.tsv",
                    assembly_plan_to_tsv(result["assembly_plan"]),
                )
                zf.writestr(
                    "assembly_fragments.fasta",
                    assembly_plan_fragments_fasta(result["assembly_plan"]),
                )
                zf.writestr(
                    "assembly_primers.tsv",
                    assembly_primers_to_tsv(result["assembly_plan"]),
                )
                zf.writestr(
                    "assembly_plan_summary.json",
                    json.dumps(assembly_plan_summary(result["assembly_plan"]), indent=2),
                )

        zip_buffer.seek(0)

        response = send_file(
            zip_buffer,
            mimetype="application/zip",
            as_attachment=True,
            download_name="sytogen_output.zip",
        )
        # So the page can warn immediately without unzipping the download
        # to find new_motifs_check.json — see the new-motifs banner logic
        # in sytogen.html.
        response.headers["X-New-Motifs-Introduced"] = str(result["summary"]["new_motifs_introduced"])
        return response
    except ValueError as e:
        return jsonify(error=str(e)), 400
    except Exception as e:
        traceback.print_exc()
        return jsonify(error=str(e)), 500


# =========================================================
# Sytogen async submit endpoint (was entirely missing)
# =========================================================

@api.route("/sytogen/submit", methods=["POST"])
def submit_sytogen():
    source_type = request.form.get("source_type", "genbank").lower()
    if source_type not in {"genbank", "fasta"}:
        return jsonify(error="source_type must be 'genbank' or 'fasta'"), 400

    gbk_file    = request.files.get("genbank")
    fasta_file  = request.files.get("fasta_file")
    gff_file    = request.files.get("gff_file")
    codon_file  = request.files.get("codon_usage")
    motif_file  = request.files.get("motif_table")

    if source_type == "genbank":
        if not gbk_file:
            return jsonify(error="Missing uploaded GenBank file"), 400
    else:
        if not fasta_file or not gff_file:
            return jsonify(error="FASTA + GFF3 mode requires both files"), 400

    if not codon_file or not motif_file:
        return jsonify(error="Missing uploaded files"), 400

    topology = request.form.get("topology", "circular").lower()
    if topology not in {"circular", "linear"}:
        return jsonify(error="topology must be 'circular' or 'linear'"), 400

    tmpdir = tempfile.mkdtemp(prefix="sytogen_")

    # Save uploads to disk so the worker thread can read them
    codon_path = os.path.join(tmpdir, secure_filename(codon_file.filename))
    motif_path = os.path.join(tmpdir, secure_filename(motif_file.filename))
    codon_file.save(codon_path)
    motif_file.save(motif_path)

    gbk_path = fasta_path = gff_path = None
    if source_type == "genbank":
        gbk_path = os.path.join(tmpdir, secure_filename(gbk_file.filename))
        gbk_file.save(gbk_path)
    else:
        fasta_path = os.path.join(tmpdir, secure_filename(fasta_file.filename))
        gff_path   = os.path.join(tmpdir, secure_filename(gff_file.filename))
        fasta_file.save(fasta_path)
        gff_file.save(gff_path)

    backbone_file = request.files.get("backbone")
    backbone_path = None
    if backbone_file and backbone_file.filename:
        backbone_path = os.path.join(tmpdir, secure_filename(backbone_file.filename))
        backbone_file.save(backbone_path)

    job_id = str(uuid.uuid4())
    JOBS[job_id] = {"status": "queued"}

    paths = {
        "genbank":     gbk_path,     # None in fasta mode
        "fasta_file":  fasta_path,   # None in genbank mode
        "gff_file":    gff_path,     # None in genbank mode
        "codon_usage": codon_path,
        "motif_table": motif_path,
        "backbone":    backbone_path,   # None if not provided
    }
    params = {
        "source_type":           source_type,
        "topology":              topology,
        "preserve_gc":           request.form.get("preserve_gc") == "true",
        "include_assembly_plan": request.form.get("include_assembly_plan") == "true",
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
        mimetype="application/zip",
        as_attachment=True,
        download_name="sytogen_output.zip",
    )
