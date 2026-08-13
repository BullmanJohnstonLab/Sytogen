"""Command-line interface for the modern SyToGen workflows."""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import zipfile
from copy import deepcopy
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

from sytogen.io import build_record_from_fasta_gff, validate_construct_size
from sytogen.motif_io import motif_table_records, parse_motif_text
from sytogen.scripts.codon_bias_estimator import run_codon_bias
from sytogen.scripts.motiffinder_backend import (
    hits_to_gff3,
    hits_to_tsv,
    load_gff3_features,
    search_motifs,
)


def _package_version() -> str:
    try:
        return version("sytogen")
    except PackageNotFoundError:
        return "0.1.0"
from sytogen.scripts.sytogen_runner import (
    assembly_plan_fragments_fasta,
    assembly_plan_summary,
    assembly_plan_to_tsv,
    assembly_primers_to_tsv,
    decision_matrix_to_tsv,
    motif_summary_to_tsv,
    run_sytogen_pipeline,
)


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8-sig")


def _write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def _parse_sequence(args: argparse.Namespace):
    if args.source_type == "genbank":
        with args.sequence.open(encoding="utf-8") as handle:
            record = SeqIO.read(handle, "genbank")
        return validate_construct_size(record)

    if args.gff is None:
        raise ValueError("FASTA mode requires --gff.")
    return build_record_from_fasta_gff(
        _read_text(args.sequence),
        _read_text(args.gff),
        args.topology,
    )


def _add_mutation_features(record, mutations):
    output_record = deepcopy(record)
    output_record.seq = Seq(str(output_record.seq))
    output_record.id = record.id
    output_record.name = record.name
    output_record.description = f"{record.description} | SyToGen result"
    for mutation in mutations:
        output_record.features.append(
            SeqFeature(
                FeatureLocation(mutation.position, mutation.position + len(mutation.new)),
                type="SyT",
                qualifiers={"label": [f"{mutation.old} --> {mutation.new}"],},
            )
        )
    return output_record


def _write_sytogen_outputs(output_dir: Path, result, input_record, motif_df, zip_path: Path | None):
    output_dir.mkdir(parents=True, exist_ok=True)
    output_record = _add_mutation_features(input_record, result["applied_mutations"])
    files = {
        "sytogen_result.fasta": result["altered_fasta"],
        "sytogen_result.gbk": output_record.format("genbank"),
        "original_sequence.fasta": result["original_fasta"],
        "input_sequence.gbk": input_record.format("genbank"),
        "motifs_used.tsv": motif_df.to_csv(sep="\t", index=False),
        "motif_summary.tsv": motif_summary_to_tsv(result["motif_summary"]),
        "decision_matrix.tsv": decision_matrix_to_tsv(result["decision_matrix"]),
        "summary.json": json.dumps(result["summary"], indent=2) + "\n",
        "new_motifs_check.json": json.dumps(
            {
                "new_motifs_introduced": result["summary"]["new_motifs_introduced"],
                "new_motifs": result["new_motifs"],
            },
            indent=2,
        )
        + "\n",
    }
    if result.get("assembly_plan"):
        plan = result["assembly_plan"]
        files.update(
            {
                "assembly_plan.tsv": assembly_plan_to_tsv(plan),
                "assembly_fragments.fasta": assembly_plan_fragments_fasta(plan),
                "assembly_primers.tsv": assembly_primers_to_tsv(plan),
                "assembly_plan_summary.json": json.dumps(assembly_plan_summary(plan), indent=2) + "\n",
            }
        )

    for name, content in files.items():
        _write_text(output_dir / name, content)

    if zip_path:
        zip_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
            for name in files:
                archive.write(output_dir / name, name)


def run_command(args: argparse.Namespace) -> int:
    sequence_record = _parse_sequence(args)
    codon_df = pd.read_csv(args.codon_usage, sep=None, engine="python")
    motif_df = parse_motif_text(_read_text(args.motifs))
    params = {
        "topology": args.topology,
        "preserve_gc": args.preserve_gc,
        "include_assembly_plan": args.assembly_plan,
        "mask_ranges": args.mask_ranges,
        "protected_override_ranges": args.protected_override_ranges,
    }
    if args.fragment_size is not None:
        params["fragment_size"] = args.fragment_size
    if args.overlap_length is not None:
        params["overlap_length"] = args.overlap_length

    result = run_sytogen_pipeline(sequence_record, codon_df, motif_df, params=params)
    _write_sytogen_outputs(args.output_dir, result, sequence_record, motif_df, args.zip)
    summary = {
        **result["summary"],
        "output_dir": str(args.output_dir),
        "zip": str(args.zip) if args.zip else None,
    }
    print(json.dumps(summary if args.json else result["summary"], indent=2))
    return 0


def parse_mymotifs_command(args: argparse.Namespace) -> int:
    records = []
    errors = []
    for path in args.inputs:
        try:
            parsed = motif_table_records(parse_motif_text(_read_text(path)))
            if not parsed:
                raise ValueError("No recognition motifs were found.")
            records.extend(parsed)
        except (UnicodeDecodeError, ValueError, pd.errors.ParserError) as exc:
            errors.append(f"{path}: {exc}")

    if not records:
        raise ValueError("No motifs could be imported." + (f" {'; '.join(errors)}" if errors else ""))

    if args.format == "json":
        content = json.dumps({"motifs": records, "file_errors": errors}, indent=2) + "\n"
    else:
        fields = ["rec_seq", "enz_type", "meth_base", "meth_type", "comp_meth_base", "comp_meth_type"]
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
        content = buffer.getvalue()

    summary = {
        "input_files": [str(path) for path in args.inputs],
        "motifs": len(records),
        "file_errors": errors,
        "output": None if args.output == Path("-") else str(args.output),
    }
    if args.json and args.output == Path("-"):
        print(json.dumps(summary, indent=2))
        return 0

    if args.output == Path("-"):
        sys.stdout.write(content)
    else:
        _write_text(args.output, content)
    if errors:
        print("Warning: " + "; ".join(errors), file=sys.stderr)
    if args.json:
        print(json.dumps(summary, indent=2))
    return 0


def codon_bias_command(args: argparse.Namespace) -> int:
    if bool(args.genome) == bool(args.fasta):
        raise ValueError("Provide exactly one input mode: --genome or --fasta with --gff.")
    if args.fasta and args.gff is None:
        raise ValueError("FASTA mode requires --gff.")

    result = run_codon_bias(
        genome_path=str(args.genome) if args.genome else None,
        fasta_path=str(args.fasta) if args.fasta else None,
        gff_path=str(args.gff) if args.gff else None,
        codon_table=args.codon_table,
        output_dir=str(args.output_dir),
    )
    summary = {
        "input_mode": "genbank" if args.genome else "fasta_gff3",
        "codon_table": args.codon_table,
        "output_dir": str(args.output_dir),
        "outputs": result,
    }
    _write_text(args.output_dir / "summary.json", json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary if args.json else result, indent=2))
    return 0


def _genbank_features(record) -> list[dict]:
    seqid = record.id or record.name or "sequence"
    features = []
    for feature in record.features:
        location = feature.location
        label = (
            feature.qualifiers.get("note", [None])[0]
            or feature.qualifiers.get("gene", [None])[0]
            or feature.qualifiers.get("label", [None])[0]
            or feature.type
        )
        features.append(
            {
                "seqid": seqid,
                "source": "GenBank",
                "type": feature.type,
                "start": int(location.start) + 1,
                "end": int(location.end),
                "score": ".",
                "strand": "+" if location.strand == 1 else "-" if location.strand == -1 else ".",
                "phase": ".",
                "attrs": f"ID={feature.type}_{int(location.start) + 1}_{int(location.end)};Name={str(label).replace(';', ',')}",
            }
        )
    return features


def motif_finder_command(args: argparse.Namespace) -> int:
    if args.source_type == "genbank":
        with args.sequence.open(encoding="utf-8") as handle:
            records = list(SeqIO.parse(handle, "genbank"))
        features = [feature for record in records for feature in _genbank_features(record)]
    else:
        records = list(SeqIO.parse(str(args.sequence), "fasta"))
        features = load_gff3_features(_read_text(args.gff)) if args.gff else []

    if not records:
        raise ValueError("No sequences found.")
    records = [validate_construct_size(record) for record in records]
    motif_df = parse_motif_text(_read_text(args.motifs))
    motifs = motif_table_records(motif_df)
    if not motifs:
        raise ValueError("No valid motifs found.")

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    tsv_parts = []
    gff_parts = []
    record_summaries = []
    for index, record in enumerate(records):
        seqid = record.id or record.name or f"sequence_{index + 1}"
        hits = search_motifs(str(record.seq), motifs, is_circular=args.topology == "circular")
        record_features = [feature for feature in features if feature["seqid"] == seqid]
        tsv = hits_to_tsv(hits, seqid, record_features)
        gff = hits_to_gff3(hits, seqid, len(record.seq), record_features)
        tsv_parts.append(tsv if index == 0 else "\n".join(tsv.splitlines()[1:]) + ("\n" if hits else ""))
        gff_parts.append(
            gff
            if index == 0
            else "".join(
                line
                for line in gff.splitlines(keepends=True)
                if not line.startswith("##gff-version") and not line.startswith("# Generated by")
            )
        )
        record_summaries.append(
            {
                "sequence_id": seqid,
                "sequence_length": len(record.seq),
                "hits": len(hits),
            }
        )

    total_hits = sum(summary["hits"] for summary in record_summaries)
    _write_text(output_dir / "motif_hits.tsv", "".join(tsv_parts))
    _write_text(output_dir / "motif_hits.gff3", "".join(gff_parts))
    summary = {
        "topology": args.topology,
        "motifs_input": len(motifs),
        "records": record_summaries,
        "hits": total_hits,
        **(
            {
                "sequence_id": record_summaries[0]["sequence_id"],
                "sequence_length": record_summaries[0]["sequence_length"],
            }
            if len(record_summaries) == 1
            else {}
        ),
    }
    _write_text(output_dir / "summary.json", json.dumps(summary, indent=2) + "\n")
    short_summary = {"records": len(records), "hits": total_hits, "output_dir": str(output_dir)}
    print(json.dumps(summary if args.json else short_summary, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sytogen", description="Run SyToGen workflows from the command line.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {_package_version()}")
    commands = parser.add_subparsers(dest="command", required=True)

    motifs = commands.add_parser("mymotifs", help="Parse motif profile files.")
    motifs_commands = motifs.add_subparsers(dest="mymotifs_command", required=True)
    parse = motifs_commands.add_parser("parse", help="Normalize CSV, TSV, REBASE, or MIJAMP motif files.")
    parse.add_argument("inputs", type=Path, nargs="+", help="Input motif files.")
    parse.add_argument("-o", "--output", type=Path, default=Path("-"), help="Output CSV/JSON path, or - for stdout.")
    parse.add_argument("--format", choices=("csv", "json"), default="csv")
    parse.add_argument("--json", action="store_true", help="Print a JSON processing summary.")
    parse.set_defaults(func=parse_mymotifs_command)

    run = commands.add_parser("run", help="Run the SyToGen optimization pipeline.")
    run.add_argument("--sequence", type=Path, required=True, help="Input GenBank or FASTA sequence.")
    run.add_argument("--source-type", choices=("genbank", "fasta"), default="genbank")
    run.add_argument("--gff", type=Path, help="GFF3 annotation file for FASTA input.")
    run.add_argument("--codon-usage", type=Path, required=True, help="Codon usage CSV/TSV.")
    run.add_argument("--motifs", type=Path, required=True, help="Motif table, REBASE, or MIJAMP file.")
    run.add_argument("--topology", choices=("circular", "linear"), default="circular")
    run.add_argument("--preserve-gc", action="store_true")
    run.add_argument("--assembly-plan", action="store_true")
    run.add_argument("--mask-ranges", default="")
    run.add_argument("--protected-override-ranges", default="")
    run.add_argument("--fragment-size", type=int)
    run.add_argument("--overlap-length", type=int)
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--zip", type=Path, help="Also write a ZIP containing the output artifacts.")
    run.add_argument("--json", action="store_true", help="Print the complete JSON summary.")
    run.set_defaults(func=run_command)

    codon_bias = commands.add_parser("codon-bias", help="Build a strain-specific codon-usage table.")
    codon_input = codon_bias.add_mutually_exclusive_group(required=True)
    codon_input.add_argument("--genome", type=Path, help="Annotated GenBank genome.")
    codon_input.add_argument("--fasta", type=Path, help="Genome FASTA paired with --gff.")
    codon_bias.add_argument("--gff", type=Path, help="GFF3 annotation file for FASTA input.")
    codon_bias.add_argument("--codon-table", type=int, default=11)
    codon_bias.add_argument("--output-dir", type=Path, required=True)
    codon_bias.add_argument("--json", action="store_true", help="Print the complete JSON summary.")
    codon_bias.set_defaults(func=codon_bias_command)

    motif_finder = commands.add_parser("motif-finder", help="Find restriction motifs on both strands.")
    motif_finder.add_argument("--sequence", type=Path, required=True)
    motif_finder.add_argument("--source-type", choices=("genbank", "fasta"), default="genbank")
    motif_finder.add_argument("--gff", type=Path, help="GFF3 annotation file for FASTA input.")
    motif_finder.add_argument("--motifs", type=Path, required=True)
    motif_finder.add_argument("--topology", choices=("circular", "linear"), default="circular")
    motif_finder.add_argument("--output-dir", type=Path, required=True)
    motif_finder.add_argument("--json", action="store_true", help="Print the complete JSON summary.")
    motif_finder.set_defaults(func=motif_finder_command)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        print(f"sytogen: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())