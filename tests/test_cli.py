import csv
import json
from pathlib import Path

from sytogen.cli import main


FIXTURES = Path(__file__).parent / "fixtures"


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