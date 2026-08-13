"""Command-line interface for the modern SyToGen workflows."""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
import zipfile
from copy import deepcopy
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

from sytogen.api import (
    build_record_from_fasta_gff,
    motif_table_records,
    parse_motif_text,
    validate_construct_size,
)
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
    output_record.id = f"{record.id}_sytogen"
    output_record.name = f"{record.name}_sytogen"
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
    print(json.dumps(result["summary"], indent=2))
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

    if args.output == Path("-"):
        sys.stdout.write(content)
    else:
        _write_text(args.output, content)
    if errors:
        print("Warning: " + "; ".join(errors), file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="sytogen", description="Run SyToGen workflows from the command line.")
    commands = parser.add_subparsers(dest="command", required=True)

    motifs = commands.add_parser("mymotifs", help="Parse motif profile files.")
    motifs_commands = motifs.add_subparsers(dest="mymotifs_command", required=True)
    parse = motifs_commands.add_parser("parse", help="Normalize CSV, TSV, REBASE, or MIJAMP motif files.")
    parse.add_argument("inputs", type=Path, nargs="+", help="Input motif files.")
    parse.add_argument("-o", "--output", type=Path, default=Path("-"), help="Output CSV/JSON path, or - for stdout.")
    parse.add_argument("--format", choices=("csv", "json"), default="csv")
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
    run.set_defaults(func=run_command)
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