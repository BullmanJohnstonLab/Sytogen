from io import BytesIO
from pathlib import Path
import sys
import base64
import csv
from zipfile import ZipFile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sytogen import create_app


FIXTURES = Path(__file__).parent / "fixtures"


def test_sytogen_page_exposes_required_workflow_controls():
    client = create_app().test_client()

    response = client.get("/sytogen")

    assert response.status_code == 200
    html = response.get_data(as_text=True)
    for token in (
        "genbank_input",
        "codon_input",
        "motif_input",
        "topology_value",
        "/api/sytogen/run",
    ):
        assert token in html


def test_sytogen_run_accepts_companion_tool_outputs():
    client = create_app().test_client()

    with (
        open(FIXTURES / "motiffinder_pEPSA5" / "motiffinder_annotated.gbk", "rb") as genbank,
        open(FIXTURES / "codonbias_pepSA5" / "codon_usage_table.csv", "rb") as codon_usage,
        open(FIXTURES / "motiffinder_pEPSA5" / "motiffinder_summary.tsv", "rb") as motif_table,
    ):
        response = client.post(
            "/api/sytogen/run",
            data={
                "genbank": (genbank, "motiffinder_annotated.gbk"),
                "codon_usage": (codon_usage, "codon_usage_table.csv"),
                "motif_table": (motif_table, "motiffinder_summary.tsv"),
                "topology": "circular",
            },
            content_type="multipart/form-data",
        )

    assert response.status_code == 200
    assert response.content_type == "application/json"

    # The API returns JSON with base64-encoded zip and plot data
    data = response.get_json()
    assert "zip_base64" in data
    assert "plot_after" in data
    assert "summary" in data

    # Decode the base64 zip and verify contents
    zip_bytes = base64.b64decode(data["zip_base64"])
    with ZipFile(BytesIO(zip_bytes)) as archive:
        assert {
            "sytogen_result.fasta",
            "sytogen_result.gbk",
            "original_sequence.fasta",
            "input_sequence.gbk",
            "motifs_used.tsv",
            "decision_matrix.tsv",
            "summary.json",
        }.issubset(set(archive.namelist()))


def test_sytogen_type_iv_motifs_are_skipped_and_marked_unchanged():
    client = create_app().test_client()

    motif_table_text = "motif\tenz_type\nATGC\t4\n"

    with (
        open(FIXTURES / "motiffinder_pEPSA5" / "motiffinder_annotated.gbk", "rb") as genbank,
        open(FIXTURES / "codonbias_pepSA5" / "codon_usage_table.csv", "rb") as codon_usage,
    ):
        response = client.post(
            "/api/sytogen/run",
            data={
                "genbank": (genbank, "motiffinder_annotated.gbk"),
                "codon_usage": (codon_usage, "codon_usage_table.csv"),
                "motif_table": (BytesIO(motif_table_text.encode("utf-8")), "motif_table.tsv"),
                "topology": "circular",
            },
            content_type="multipart/form-data",
        )

    assert response.status_code == 200
    payload = response.get_json()

    summary = payload["summary"]
    assert summary["edits_applied"] == 0
    assert summary["motifs_resolved"] == 0
    assert summary["motifs_unresolved"] == summary["motifs_input"]

    zip_bytes = base64.b64decode(payload["zip_base64"])
    with ZipFile(BytesIO(zip_bytes)) as archive:
        matrix_tsv = archive.read("decision_matrix.tsv").decode("utf-8")

    rows = list(csv.DictReader(matrix_tsv.splitlines(), delimiter="\t"))
    type_iv_rows = [r for r in rows if r.get("skip_reason") == "type_iv_skipped"]

    assert type_iv_rows, "Expected at least one Type IV skip row in decision matrix"
    assert all(r.get("chosen") in ("False", "false", "") for r in type_iv_rows)
    assert all((r.get("edit_position") or "") == "" for r in type_iv_rows)
    assert all("Type IV motif" in (r.get("reasoning") or "") for r in type_iv_rows)
