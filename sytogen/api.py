import io
import re
import json
import os
import uuid
import shutil
import time
import zipfile
import base64
import tempfile
import traceback
import pandas as pd
import copy

from threading import Thread, Lock, Semaphore
from concurrent.futures import ThreadPoolExecutor

from flask import (
    Blueprint,
    request,
    send_file,
    abort,
    jsonify,
    after_this_request,
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
from sytogen.scripts.rebase_motif_parser import parse_rebase_motif_file
from sytogen.scripts.sytogen_runner import (
    run_sytogen_pipeline,
    decision_matrix_to_tsv,
    motif_summary_to_tsv,
    assembly_plan_to_tsv,
    assembly_plan_fragments_fasta,
    assembly_plan_summary,
    assembly_primers_to_tsv,
)
from sytogen.scripts.visualization import build_plasmid_maps, build_motiffinder_map
from sytogen import job_store

# =========================================================
# Blueprint
# =========================================================

api = Blueprint("api", __name__)

# Job state (status/result/tmpdir/etc for the async /sytogen pipeline) is
# stored in sytogen/job_store.py (SQLite on shared disk) rather than an
# in-process dict, so it works correctly across multiple worker processes.

# How long a job (and its temp directory) is kept after it finishes if the
# client never calls /sytogen/result/<job_id> to download it. Downloaded
# jobs are already cleaned up immediately (see the result() endpoint) -
# this TTL is the safety net for jobs that are submitted and then
# abandoned. Both values are overridable per-deployment via env vars.
JOB_TTL_SECONDS = int(os.environ.get("SYTOGEN_JOB_TTL_SECONDS", 60 * 60))
JOB_SWEEP_INTERVAL_SECONDS = int(os.environ.get("SYTOGEN_JOB_SWEEP_INTERVAL_SECONDS", 5 * 60))

# The sytogen pipeline is CPU-intensive (regex motif scanning + candidate
# generation across the whole construct). Previously each /sytogen/submit
# spawned a brand new, unbounded Thread that started that work immediately -
# a burst of submissions could spawn unboundedly many concurrently-running
# CPU-bound threads and degrade or take down the whole process. A fixed-size
# pool caps how many jobs actually run at once per worker process; a
# semaphore separately caps how many can be queued+running combined, so a
# flood of submissions gets a clear "try again later" instead of piling up
# an unbounded backlog. Both are per-process, so with N gunicorn workers the
# real ceiling is N * SYTOGEN_MAX_CONCURRENT_JOBS - size accordingly.
MAX_CONCURRENT_JOBS = int(os.environ.get("SYTOGEN_MAX_CONCURRENT_JOBS", 4))
MAX_QUEUED_JOBS = int(os.environ.get("SYTOGEN_MAX_QUEUED_JOBS", 20))

JOB_EXECUTOR = ThreadPoolExecutor(
    max_workers=MAX_CONCURRENT_JOBS,
    thread_name_prefix="sytogen-worker",
)
JOB_ADMISSION = Semaphore(MAX_QUEUED_JOBS)

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

# The interactive motif and redesign workflows are intended for plasmids and
# similarly sized constructs. Keeping this cap explicit prevents accidental
# whole-genome uploads from exhausting the synchronous analysis endpoints.
MAX_CONSTRUCT_LENGTH = 20_000


# =========================================================
# Helpers
# =========================================================

def allowed_extension(filename, allowed):

    ext = os.path.splitext(filename)[1].lower()

    return ext in allowed


def validate_construct_size(record):
    """Raise a clear validation error when a construct exceeds 20 kb."""
    sequence_length = len(record.seq)
    if sequence_length > MAX_CONSTRUCT_LENGTH:
        name = record.id or record.name or "Uploaded construct"
        raise ValueError(
            f"{name} is {sequence_length:,} bp. The maximum supported construct "
            f"size is {MAX_CONSTRUCT_LENGTH:,} bp (20 kb)."
        )
    return record


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
    return validate_construct_size(record)


def read_uploaded_table(file_storage):
    text = file_storage.stream.read().decode("utf-8-sig")
    return pd.read_csv(
        io.StringIO(text),
        sep=None,
        engine="python",
    )


_MIJAMP_MOTIF_RE = re.compile(
    r"^(?P<prefix>[ACGTURYKMSWBDHVN]+)\((?P<code>[456]m[AGC])\)(?P<suffix>[ACGTURYKMSWBDHVN]+)$",
    re.IGNORECASE,
)


def _parse_mijamp_motif_token(token):
    match = _MIJAMP_MOTIF_RE.match(str(token or "").strip())
    if not match:
        return None

    code = match.group("code").strip().upper()
    meth_type = f"m{code[0]}{code[-1]}"
    meth_base_char = code[-1]

    prefix = match.group("prefix").upper()
    suffix = match.group("suffix").upper()
    rec_seq = prefix + meth_base_char + suffix

    plus_position = len(prefix) + 1
    minus_position = len(rec_seq) - plus_position + 1

    if "N" in rec_seq or any(char in rec_seq for char in "RYKMSWBDHV"):
        enz_type = "1"
    else:
        enz_type = "2"

    return {
        "rec_seq": rec_seq,
        "enz_type": enz_type,
        "meth_base": str(plus_position),
        "meth_type": meth_type or "-",
        "comp_meth_base": str(minus_position),
        "comp_meth_type": meth_type or "-",
    }


def parse_mijamp_expected_output(text):
    rows = []
    for line in (text or "").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue

        motif_token = stripped.split("\t", 1)[0].strip()
        parsed = _parse_mijamp_motif_token(motif_token)
        if parsed:
            rows.append(parsed)

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(
        rows,
        columns=["rec_seq", "enz_type", "meth_base", "meth_type", "comp_meth_base", "comp_meth_type"],
    )


def parse_motif_text(text):
    """
    Core motif-table parsing logic, operating on raw text so it can be
    shared by both the synchronous upload path (read_motif_table, below)
    and the async worker path (which reads the same file back from disk).

    Returns a DataFrame with a 'motif' column, as sytogen_runner._parse_motifs()
    expects. Accepts either a plain delimited table with a motif-like column,
    REBASE-style tagged exports (e.g. "<enz_type>2<rec_seq>ATGC...<>"), or
    simple known-enzyme names / MIJAMP-style tables where motif sequences are
    stored in a column such as 'Motif', 'Recognition sequence', or 'Sequence',
    including MIJAMP expected-output files with parenthesized methylation marks.
    """

    def looks_like_motif(value):
        if value is None:
            return False
        if isinstance(value, float) and pd.isna(value):
            return False
        cleaned = str(value).strip()
        if not cleaned or cleaned in {"-", "NA", "N/A", "nan", "None"}:
            return False
        letters = re.sub(r"[^ACGTURYKMSWBDHVN]", "", cleaned.upper())
        return 2 <= len(letters) <= 40 and len(letters) / max(1, len(cleaned)) >= 0.6

    # --- Attempt 1: MIJAMP expected-output format ---
    mijamp_df = parse_mijamp_expected_output(text)
    if not mijamp_df.empty:
        return mijamp_df

    # --- Attempt 2: plain delimited motif table ---
    try:
        df = pd.read_csv(io.StringIO(text), sep=None, engine="python")
        if df.empty:
            raise ValueError("empty table")

        normalized_cols = {str(c).strip().lower(): c for c in df.columns}
        for candidate in ("motif", "rec_seq", "recognition_motif", "recognition_sequence", "sequence", "seq"):
            if candidate not in normalized_cols:
                continue
            source = normalized_cols[candidate]
            if source != "motif":
                df = df.rename(columns={source: "motif"})
            return df

        for column in df.columns:
            normalized = str(column).strip().lower()
            if any(token in normalized for token in ("motif", "recognition", "sequence", "site", "target")):
                values = df[column].dropna().astype(str)
                motif_like = [v for v in values if looks_like_motif(v)]
                if motif_like and len(motif_like) >= max(1, min(3, len(values))):
                    df = df.rename(columns={column: "motif"})
                    return df

        for column in df.columns:
            values = df[column].dropna().astype(str)
            motif_like = [v for v in values if looks_like_motif(v)]
            if motif_like and len(motif_like) >= max(1, min(3, len(values))):
                df = df.rename(columns={column: "motif"})
                return df
    except Exception:
        pass  # not a plain delimited table — fall through to REBASE parsing

    # --- Attempt 3: REBASE-style tagged export or simple known-enzyme names ---
    motif_df = parse_rebase_motif_file(text, is_path=False, drop_unclassified=False)
    if motif_df.empty:
        raise ValueError(
            "Could not parse the restriction motif table. Expected either "
            "a delimited file with a 'motif' column, or a REBASE-style "
            "tagged export (e.g. containing '<rec_seq>...' entries), or "
            "plain enzyme names such as 'EcoRI' or 'BamHI'."
        )

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


def summarize_motif_hits(hits):
    """Return a compact, stable summary for the MotifFinder result panel."""
    grouped = {}
    for hit in hits:
        key = (hit["rec_seq"], str(hit.get("enz_type") or "-1"))
        row = grouped.setdefault(key, {
            "motif": key[0],
            "enzyme_type": key[1],
            "hits": 0,
            "forward_hits": 0,
            "reverse_hits": 0,
        })
        row["hits"] += 1
        if hit.get("strand") == "+":
            row["forward_hits"] += 1
        else:
            row["reverse_hits"] += 1
    return sorted(grouped.values(), key=lambda row: (-row["hits"], row["motif"], row["enzyme_type"]))


def motif_table_records(motif_df):
    """Convert a parsed motif table to the fields used by MyMotif's editor."""
    aliases = {
        "rec_seq": ("rec_seq", "motif", "recognition_motif", "recognition_sequence", "sequence", "seq"),
        "enz_type": ("enz_type", "type"),
        "meth_base": (
            "meth_base",
            "methylated_base_plus",
            "methylated_base",
            "methylated base",
            "methylation_base",
            "methylation base",
        ),
        "meth_type": (
            "meth_type",
            "methylated_base_plus_type",
            "methylation_type",
            "methylation type",
            "methylation",
            "methylated_type",
            "methylated type",
        ),
        "comp_meth_base": (
            "comp_meth_base",
            "methylated_base_minus",
            "complementary_methylated_base",
            "complementary methylated base",
        ),
        "comp_meth_type": (
            "comp_meth_type",
            "methylated_base_minus_type",
            "complementary_methylation_type",
            "complementary methylation type",
        ),
    }
    columns = {str(column).strip().lower(): column for column in motif_df.columns}

    def value(row, field):
        source = next((columns[name] for name in aliases[field] if name in columns), None)
        if source is None or pd.isna(row[source]):
            return ""
        return str(row[source]).strip()

    records = []
    for _, row in motif_df.iterrows():
        rec_seq = value(row, "rec_seq").upper()
        if rec_seq:
            records.append({
                "rec_seq": rec_seq,
                "enz_type": value(row, "enz_type") or "-1",
                "meth_base": value(row, "meth_base") or "-",
                "meth_type": value(row, "meth_type") or "-",
                "comp_meth_base": value(row, "comp_meth_base") or "-",
                "comp_meth_type": value(row, "comp_meth_type") or "-",
            })
    return records


@api.route("/mymotif/parse", methods=["POST"])
def parse_mymotif_files():
    """Parse CSV/TSV tables and REBASE tagged exports for the MyMotif editor."""
    files = request.files.getlist("motif_files")
    if not files:
        return jsonify(error="Choose at least one CSV, TSV, or REBASE motif file."), 400

    motifs = []
    errors = []
    for uploaded_file in files:
        if not uploaded_file or not uploaded_file.filename:
            continue
        try:
            text = uploaded_file.read().decode("utf-8-sig")
            parsed = motif_table_records(parse_motif_text(text))
            if not parsed:
                raise ValueError("No recognition motifs were found.")
            motifs.extend(parsed)
        except (UnicodeDecodeError, ValueError, pd.errors.ParserError) as exc:
            errors.append({"file": uploaded_file.filename, "error": str(exc)})

    if not motifs:
        return jsonify(
            error="No motifs could be imported.",
            file_errors=errors,
        ), 400

    return jsonify(motifs=motifs, file_errors=errors)


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
    response_format = request.form.get("response_format", "").lower()

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

    try:
        for record in records:
            validate_construct_size(record)
    except ValueError as exc:
        abort(400, str(exc))

    features = []

    if source_type == "genbank":

        for rec in records:

            for feat in rec.features:

                loc = feat.location

                seqid = rec.id or rec.name

                # Prefer the GenBank /note text for non-gene annotations:
                # it is the human-readable feature description (e.g. an
                # origin name), unlike a generated "misc_feature_16" ID.
                feature_label = (
                    feat.qualifiers.get("note", [None])[0]
                    or feat.qualifiers.get("gene", [None])[0]
                    or feat.qualifiers.get("label", [None])[0]
                    or feat.type
                )
                feature_label = str(feature_label).replace(";", ",")
                attrs = (
                    f"ID={feat.type}_"
                    f"{int(loc.start)+1}_"
                    f"{int(loc.end)};"
                    f"Name={feature_label}"
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
    record_plots = {}   # seqid -> plotly Figure, one per record
    record_motif_summaries = {}  # seqid -> compact per-motif hit counts

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

        topology_token = str(
            rec.annotations.get(
                "topology",
                "",
            )
        ).strip().lower()
        if not topology_token:
            # Some GenBank exporters put "circular" in LOCUS where BioPython
            # maps it to data_file_division instead of topology.
            topology_token = str(
                rec.annotations.get(
                    "data_file_division",
                    "",
                )
            ).strip().lower()

        is_circular = (
            source_type == "genbank"
            and topology_token == "circular"
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

        record_plots[seqid] = build_motiffinder_map(
            rec_features,
            hits,
            seq_len,
            "circular" if is_circular else "linear",
            title=seqid,
        )
        record_motif_summaries[seqid] = summarize_motif_hits(hits)

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

        # One self-contained interactive HTML map per record — no
        # image-export dependency needed (fig.to_html embeds Plotly.js
        # via CDN reference), same approach as SyToGen's plasmid maps.
        for plot_seqid, fig in record_plots.items():
            safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", plot_seqid)
            zf.writestr(
                f"motif_map_{safe_name}.html",
                fig.to_html(full_html=True, include_plotlyjs="cdn"),
            )

    zip_buf.seek(0)

    # First record's map for the live page (the common case is a single
    # sequence per run); every record still gets its own HTML file above.
    first_seqid = records[0].id or records[0].name or "unknown"
    first_plot = record_plots.get(first_seqid)

    if response_format == "json":
        return jsonify({
            "zip_base64": base64.b64encode(zip_buf.getvalue()).decode("ascii"),
            "plot": json.loads(first_plot.to_json()) if first_plot else None,
            "motif_summary": record_motif_summaries.get(first_seqid, []),
            "annotated_gbk": gbk_buf.getvalue(),
        })

    return jsonify({
        "zip_base64": base64.b64encode(zip_buf.getvalue()).decode("ascii"),
        "plot": json.loads(first_plot.to_json()) if first_plot else None,
        "motif_summary": record_motif_summaries.get(first_seqid, []),
    })


# =========================================================
# CodonBias endpoint
# =========================================================

@api.route("/codonbias/run", methods=["POST"])
def run_codonbias():

    try:

        response_format = request.form.get("response_format", "").lower()

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

    with tempfile.TemporaryDirectory(prefix="codonbias_") as tmpdir:
        return _run_codonbias_in(tmpdir, source_type, codon_table, response_format)


def _run_codonbias_in(tmpdir, source_type, codon_table, response_format):
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

        if response_format == "json":
            with open(zip_path, "rb") as zip_handle:
                zip_bytes = zip_handle.read()

            with open(output_paths["csv"], "r", encoding="utf-8-sig") as csv_handle:
                codon_usage_csv = csv_handle.read()

            return jsonify({
                "zip_base64": base64.b64encode(zip_bytes).decode("ascii"),
                "codon_usage_csv": codon_usage_csv,
            })

        with open(zip_path, "rb") as zip_handle:
            zip_bytes = zip_handle.read()

        return send_file(
            io.BytesIO(zip_bytes),
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
        ), 500


# =========================================================
# Worker
# =========================================================

def worker(job_id, paths, params, tmpdir):

    job_store.update_job(job_id, tmpdir=tmpdir)

    try:
        job_store.update_job(job_id, status="running")

        # Parse inputs the same way the synchronous /sytogen/run endpoint
        # does, reusing the same helpers (including the REBASE-format
        # motif-table fallback) so both code paths behave identically.
        source_type = params.get("source_type", "genbank")
        if source_type == "genbank":
            seq_record = validate_construct_size(SeqIO.read(paths["genbank"], "genbank"))
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

        # Build params dict, including optional assembly options (use None for unspecified)
        pipeline_params = {
            "topology":              params.get("topology", "circular"),
            "preserve_gc":           params.get("preserve_gc", False),
            "include_assembly_plan": params.get("include_assembly_plan", False),
            "mask_ranges":           params.get("mask_ranges", ""),
            "protected_override_ranges": params.get("protected_override_ranges", ""),
        }
        
        # Add optional assembly parameters if provided
        if "fragment_size" in params:
            try:
                pipeline_params["fragment_size"] = int(params["fragment_size"])
            except (ValueError, TypeError):
                pass
        if "overlap_length" in params:
            try:
                pipeline_params["overlap_length"] = int(params["overlap_length"])
            except (ValueError, TypeError):
                pass
        
        # Note: target_gc, target_tm, primer_* parameters are reserved for future use
        # Currently they require modifying assembly_planner constants
        
        result = run_sytogen_pipeline(
            seq_record,
            codon_df,
            motif_df,
            params=pipeline_params,
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

        fig_before, fig_after = build_plasmid_maps(
            output_record,
            result["motifs"],
            result["new_motifs"],
            result["decision_matrix"],
            result["resolved_motif_keys"],
            len(result["altered_sequence"]),
            params.get("topology", "circular"),
            mask_regions=result["mask_regions"],
            protected_override_ranges=result.get("protected_override_ranges", []),
            title=seq_record.id,
        )

        zip_path = os.path.join(tmpdir, "sytogen_output.zip")
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("sytogen_result.fasta",    result["altered_fasta"])
            zf.writestr("sytogen_result.gbk",      output_record.format("genbank"))
            zf.writestr("original_sequence.fasta", result["original_fasta"])
            zf.writestr("input_sequence.gbk",      seq_record.format("genbank"))
            zf.writestr("motifs_used.tsv",         motifs_used)
            zf.writestr("motif_summary.tsv",       motif_summary_to_tsv(result["motif_summary"]))
            zf.writestr(
                "decision_matrix.tsv",
                decision_matrix_to_tsv(result["decision_matrix"]),
            )
            zf.writestr(
                "plasmid_map_before.html",
                fig_before.to_html(full_html=True, include_plotlyjs="cdn"),
            )
            zf.writestr(
                "plasmid_map_after.html",
                fig_after.to_html(full_html=True, include_plotlyjs="cdn"),
            )
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

        job_store.update_job(
            job_id,
            status="done",
            result=zip_path,
            finished_at=time.time(),
        )

    except Exception as e:

        traceback.print_exc()

        job_store.update_job(
            job_id,
            status="error",
            error=str(e),
            traceback=traceback.format_exc(),
            finished_at=time.time(),
        )

        shutil.rmtree(tmpdir, ignore_errors=True)
        job_store.update_job(job_id, tmpdir=None)


def _run_worker_and_release_slot(job_id, paths, params, tmpdir):
    """Run a job on the bounded JOB_EXECUTOR pool, then release its
    JOB_ADMISSION slot. The slot is held for the job's whole lifetime
    (queued in the pool + actually running), not just released once it's
    handed off - otherwise JOB_ADMISSION would only bound how fast jobs
    get *submitted*, not how many are actually queued or in flight."""
    try:
        worker(job_id, paths, params, tmpdir)
    finally:
        JOB_ADMISSION.release()


# =========================================================
# Job TTL sweep
# =========================================================
# Jobs that ARE downloaded clean up their own temp directory immediately
# (see the result() endpoint below). This sweep is the safety net for
# jobs that finish but are never collected - without it, an abandoned
# job's temp directory (and its row in the shared job store) would sit
# around forever. Every worker process runs this loop independently;
# since they all share the same SQLite-backed job store, whichever one's
# timer fires first does the cleanup for all of them.

def _sweep_loop():
    while True:
        time.sleep(JOB_SWEEP_INTERVAL_SECONDS)
        try:
            job_store.sweep_expired_jobs(JOB_TTL_SECONDS)
        except Exception:
            traceback.print_exc()


_sweeper_started = False
_sweeper_lock = Lock()


def start_job_sweeper():
    """Start the background TTL sweep thread. Safe to call more than once
    per process (e.g. if create_app() runs multiple times in tests, or
    once per gunicorn worker) - only the first call in a given process
    actually starts a thread."""
    global _sweeper_started
    with _sweeper_lock:
        if _sweeper_started:
            return
        _sweeper_started = True
    Thread(target=_sweep_loop, daemon=True).start()


# =========================================================
# Status endpoint
# =========================================================

@api.route("/status/<job_id>", methods=["GET"])
def status(job_id):

    job = job_store.get_job(job_id)

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
            seq_record = validate_construct_size(
                SeqIO.read(io.TextIOWrapper(gbk_file.stream), "genbank")
            )
        else:
            fasta_text = fasta_file.stream.read().decode("utf-8-sig")
            gff_text   = gff_file.stream.read().decode("utf-8-sig")
            seq_record = build_record_from_fasta_gff(fasta_text, gff_text, topology)

        # Convert uploaded tables to DataFrames, accepting CSV or TSV output,
        # and REBASE-tagged exports for the motif table.
        codon_df = read_uploaded_table(codon_file)
        motif_df = read_motif_table(motif_file)

        # =================================================
        # RUN PIPELINE
        # =================================================

        # Build params dict, including optional assembly options
        pipeline_params = {
            "topology":              topology,
            "preserve_gc":           request.form.get("preserve_gc") == "true",
            "include_assembly_plan": request.form.get("include_assembly_plan") == "true",
            "mask_ranges":           request.form.get("mask_ranges", ""),
            "protected_override_ranges": request.form.get("protected_override_ranges", ""),
        }
        
        # Add optional assembly parameters if provided
        if "fragment_size" in request.form:
            try:
                pipeline_params["fragment_size"] = int(request.form["fragment_size"])
            except (ValueError, TypeError):
                pass
        if "overlap_length" in request.form:
            try:
                pipeline_params["overlap_length"] = int(request.form["overlap_length"])
            except (ValueError, TypeError):
                pass
        
        # Note: target_gc, target_tm, primer_* parameters are reserved for future use

        result = run_sytogen_pipeline(
            seq_record,
            codon_df,
            motif_df,
            params=pipeline_params
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

        fig_before, fig_after = build_plasmid_maps(
            output_record,
            result["motifs"],
            result["new_motifs"],
            result["decision_matrix"],
            result["resolved_motif_keys"],
            len(result["altered_sequence"]),
            topology,
            mask_regions=result["mask_regions"],
            protected_override_ranges=result.get("protected_override_ranges", []),
            title=seq_record.id,
        )

        with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("sytogen_result.fasta",    result["altered_fasta"])
            zf.writestr("sytogen_result.gbk",      output_record.format("genbank"))
            zf.writestr("original_sequence.fasta", result["original_fasta"])
            zf.writestr("input_sequence.gbk",      seq_record.format("genbank"))
            zf.writestr("motifs_used.tsv",         motifs_used)
            zf.writestr("motif_summary.tsv",       motif_summary_to_tsv(result["motif_summary"]))
            zf.writestr(
                "decision_matrix.tsv",
                decision_matrix_to_tsv(result["decision_matrix"]),
            )
            zf.writestr(
                "plasmid_map_before.html",
                fig_before.to_html(full_html=True, include_plotlyjs="cdn"),
            )
            zf.writestr(
                "plasmid_map_after.html",
                fig_after.to_html(full_html=True, include_plotlyjs="cdn"),
            )
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

        # JSON response rather than a raw zip blob: the page renders
        # plot_after live via Plotly.js, and decodes zip_base64 into a
        # Blob itself for the download button. fig_before isn't sent here
        # — the "before" map now lives on MotifFinder's own page (that's
        # the tool that actually finds the motifs SyToGen starts from) —
        # but it's still included as a self-contained interactive HTML
        # file in the zip (plasmid_map_before.html) as a reference.
        return jsonify({
            "zip_base64": base64.b64encode(zip_buffer.getvalue()).decode("ascii"),
            "plot_after": json.loads(fig_after.to_json()),
            "summary": result["summary"],
            "motif_summary": result["motif_summary"],
            "new_motifs_introduced": result["summary"]["new_motifs_introduced"],
        })
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

    if not JOB_ADMISSION.acquire(blocking=False):
        return jsonify(
            error="Server is busy processing other jobs. Please try again shortly."
        ), 503

    try:
        # Validate before persisting an async job so an oversized upload gets
        # the same immediate, actionable response as the synchronous route.
        if source_type == "genbank":
            genbank_text = gbk_file.stream.read().decode("utf-8-sig")
            validate_construct_size(SeqIO.read(io.StringIO(genbank_text), "genbank"))
            gbk_file.stream.seek(0)
        else:
            fasta_text = fasta_file.stream.read().decode("utf-8-sig")
            gff_text = gff_file.stream.read().decode("utf-8-sig")
            build_record_from_fasta_gff(fasta_text, gff_text, topology)
            fasta_file.stream.seek(0)
            gff_file.stream.seek(0)
    except (UnicodeDecodeError, ValueError) as exc:
        JOB_ADMISSION.release()
        return jsonify(error=str(exc)), 400

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

    job_id = str(uuid.uuid4())
    job_store.create_job(job_id)

    paths = {
        "genbank":     gbk_path,     # None in fasta mode
        "fasta_file":  fasta_path,   # None in genbank mode
        "gff_file":    gff_path,     # None in genbank mode
        "codon_usage": codon_path,
        "motif_table": motif_path,
    }
    params = {
        "source_type":           source_type,
        "topology":              topology,
        "preserve_gc":           request.form.get("preserve_gc") == "true",
        "include_assembly_plan": request.form.get("include_assembly_plan") == "true",
        "mask_ranges":           request.form.get("mask_ranges", ""),
        "protected_override_ranges": request.form.get("protected_override_ranges", ""),
    }

    JOB_EXECUTOR.submit(_run_worker_and_release_slot, job_id, paths, params, tmpdir)

    return jsonify(job_id=job_id), 202


# =========================================================
# Result endpoint
# =========================================================

@api.route(
    "/sytogen/result/<job_id>",
    methods=["GET"],
)
def result(job_id):

    job = job_store.get_job(job_id)

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

    tmpdir = job.get("tmpdir")

    @after_this_request
    def _cleanup(response):
        if tmpdir:
            shutil.rmtree(tmpdir, ignore_errors=True)
        job_store.delete_job(job_id)
        return response

    return send_file(
        result_path,
        mimetype="application/zip",
        as_attachment=True,
        download_name="sytogen_output.zip",
    )
