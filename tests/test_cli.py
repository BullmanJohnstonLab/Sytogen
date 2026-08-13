import csv
import json
import pytest
from pathlib import Path

from sytogen.cli import main


FIXTURES = Path(__file__).parent / "fixtures"


def test_cli_reports_package_version(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(["--version"])

    assert exit_info.value.code == 0
    assert capsys.readouterr().out.strip() == "sytogen 0.1.0"


def test_mymotifs_parse_writes_normalized_csv(tmp_path):
    output = tmp_path / "motifs.csv"

    assert main([
        "mymotifs",
        "parse",
        str(FIXTURES / "mymotif_mijamp_expectedOutput.tsv"),
        "--output",
        str(output),
    ]) == 0

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert [row["rec_seq"] for row in rows] == ["GATC", "AACNNNNNNGTGC", "CCWGG"]
    assert rows[0]["comp_meth_base"] == "3"


def test_run_writes_pipeline_artifacts(tmp_path):
    output_dir = tmp_path / "run"

    assert main([
        "run",
        "--sequence",
        str(FIXTURES / "pScout.gbk"),
        "--codon-usage",
        str(FIXTURES / "codonbias_pScout" / "codon_usage_table.csv"),
        "--motifs",
        str(FIXTURES / "motiffinder_pScout" / "motiffinder_summary.tsv"),
        "--output-dir",
        str(output_dir),
    ]) == 0

    expected = {
        "sytogen_result.fasta",
        "sytogen_result.gbk",
        "original_sequence.fasta",
        "input_sequence.gbk",
        "motifs_used.tsv",
        "motif_summary.tsv",
        "decision_matrix.tsv",
        "summary.json",
        "new_motifs_check.json",
    }
    assert {path.name for path in output_dir.iterdir()} == expected
    assert json.loads((output_dir / "summary.json").read_text())[
        "new_motifs_introduced"
    ] == 0


def test_codon_bias_writes_generated_sequence_artifacts(tmp_path):
    output_dir = tmp_path / "codon-bias"

    assert main([
        "codon-bias",
        "--genome",
        str(FIXTURES / "codonbias_pScout" / "codonbias_input.gbk"),
        "--output-dir",
        str(output_dir),
    ]) == 0

    assert {
        "codon_usage_table.csv",
        "codonbias_input.gbk",
        "codonbias_input.fasta",
        "codonbias_input.gff3",
    } == {path.name for path in output_dir.iterdir()}


def test_motif_finder_writes_hit_tables(tmp_path):
    output_dir = tmp_path / "motif-finder"

    assert main([
        "motif-finder",
        "--sequence",
        str(FIXTURES / "motiffinder_pScout" / "motiffinder_annotated.gbk"),
        "--motifs",
        str(FIXTURES / "motiffinder_pScout" / "motiffinder_summary.tsv"),
        "--topology",
        "circular",
        "--output-dir",
        str(output_dir),
    ]) == 0

    summary = json.loads((output_dir / "summary.json").read_text())
    assert summary["hits"] > 0
    assert (output_dir / "motif_hits.tsv").read_text().startswith("seqid\tposition_1based")
    assert (output_dir / "motif_hits.gff3").read_text().startswith("##gff-version 3")


def test_motif_finder_aggregates_multiple_sequence_records(tmp_path):
    sequence = tmp_path / "sequences.fasta"
    sequence.write_text(">first\nAAAATGC\n>second\nGGGATGC\n")
    motifs = tmp_path / "motifs.txt"
    motifs.write_text("ATGC\n")
    output_dir = tmp_path / "motif-finder"

    assert main([
        "motif-finder",
        "--sequence",
        str(sequence),
        "--source-type",
        "fasta",
        "--motifs",
        str(motifs),
        "--topology",
        "linear",
        "--output-dir",
        str(output_dir),
    ]) == 0

    summary = json.loads((output_dir / "summary.json").read_text())
    assert [record["sequence_id"] for record in summary["records"]] == ["first", "second"]
    assert summary["hits"] == 2
    tsv = (output_dir / "motif_hits.tsv").read_text()
    assert tsv.count("\nfirst\t") == 1
    assert tsv.count("\nsecond\t") == 1