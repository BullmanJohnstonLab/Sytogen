from io import BytesIO
from pathlib import Path
import sys
import base64
import csv
from io import StringIO
from zipfile import ZipFile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sytogen import create_app


FIXTURES = Path(__file__).parent / "fixtures"


def test_sytogen_page_exposes_required_workflow_controls():
    client = create_app().test_client()

    response = client.get("/sytogen")

    assert response.status_code == 200
    html = response.get_data(as_text=True)
    for token in (
        "workflow_tab_sytogen_only",
        "workflow_tab_end_to_end",
        "workflow_panel_end_to_end",
        "sequence-input-heading",
        "codon-usage-section",
        "sytogen-run-form",
        "/mymotif",
        "/motiffinder",
        "/codon-bias",
        "genbank_input",
        "codon_input",
        "motif_input",
        "topology_value",
    ):
        assert token in html


def test_motiffinder_can_return_json_for_chained_workflows():
    client = create_app().test_client()

    with open(FIXTURES / "motiffinder_pEPSA5" / "motiffinder_annotated.gbk", "rb") as genbank:
        response = client.post(
            "/api/motiffinder/run",
            data={
                "source_type": "genbank",
                "sequence_file": (genbank, "motiffinder_annotated.gbk"),
                "motif_file": (
                    BytesIO(b"<enz_type>2<rec_seq>ATGC<meth_base>C<>") ,
                    "motifs.txt",
                ),
                "response_format": "json",
            },
            content_type="multipart/form-data",
        )

    assert response.status_code == 200
    payload = response.get_json()
    assert payload["annotated_gbk"]
    assert payload["motif_summary"]


def test_codonbias_can_return_json_for_chained_workflows():
    client = create_app().test_client()

    with open(FIXTURES / "codonbias_pepSA5" / "codonbias_input.gbk", "rb") as genome:
        response = client.post(
            "/api/codonbias/run",
            data={
                "source_type": "genbank",
                "genome_file": (genome, "codonbias_input.gbk"),
                "response_format": "json",
            },
            content_type="multipart/form-data",
        )

    assert response.status_code == 200
    payload = response.get_json()
    assert payload["codon_usage_csv"].strip()
    assert payload["zip_base64"]


def test_sytogen_rejects_constructs_over_20kb():
    record = SeqRecord(Seq("A" * 20_001), id="oversized")
    record.annotations["molecule_type"] = "DNA"
    genbank = StringIO()
    SeqIO.write(record, genbank, "genbank")

    client = create_app().test_client()
    response = client.post(
        "/api/sytogen/run",
        data={
            "genbank": (BytesIO(genbank.getvalue().encode("utf-8")), "oversized.gbk"),
            "codon_usage": (BytesIO(b"codon,fraction\nAAA,1\n"), "codons.csv"),
            "motif_table": (BytesIO(b"motif\nenz_type\n"), "motifs.tsv"),
            "topology": "circular",
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 400
    assert "20,000 bp" in response.json["error"]


def test_motiffinder_returns_compact_motif_summary():
    client = create_app().test_client()
    with open(FIXTURES / "motiffinder_pEPSA5" / "motiffinder_annotated.gbk", "rb") as genbank:
        response = client.post(
            "/api/motiffinder/run",
            data={
                "source_type": "genbank",
                "sequence_file": (genbank, "motiffinder_annotated.gbk"),
                "motif_file": (
                    BytesIO(b"<enz_type>2<rec_seq>ATGC<meth_base>C<>"),
                    "motifs.txt",
                ),
            },
            content_type="multipart/form-data",
        )

    assert response.status_code == 200
    assert response.json["motif_summary"] == [{
        "motif": "ATGC",
        "enzyme_type": "2",
        "hits": response.json["motif_summary"][0]["hits"],
        "forward_hits": response.json["motif_summary"][0]["forward_hits"],
        "reverse_hits": response.json["motif_summary"][0]["reverse_hits"],
    }]
    assert response.json["motif_summary"][0]["hits"] > 0


def test_motiffinder_uses_circular_plot_when_genbank_locus_marks_circular():
    client = create_app().test_client()

    with open(FIXTURES / "pScout 1.gbk", "rb") as genbank:
        response = client.post(
            "/api/motiffinder/run",
            data={
                "source_type": "genbank",
                "sequence_file": (genbank, "pScout 1.gbk"),
                "motif_file": (
                    BytesIO(b"<enz_type>2<rec_seq>GATC<meth_base>C<>"),
                    "motifs.txt",
                ),
            },
            content_type="multipart/form-data",
        )

    assert response.status_code == 200
    plot = response.json["plot"]
    assert plot is not None
    assert "polar" in plot["layout"]


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
    assert "motif_summary" in data

    # Decode the base64 zip and verify contents
    zip_bytes = base64.b64decode(data["zip_base64"])
    with ZipFile(BytesIO(zip_bytes)) as archive:
        assert {
            "sytogen_result.fasta",
            "sytogen_result.gbk",
            "original_sequence.fasta",
            "input_sequence.gbk",
            "motifs_used.tsv",
            "motif_summary.tsv",
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
