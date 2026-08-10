from io import BytesIO
from pathlib import Path

from sytogen import create_app


FIXTURES = Path(__file__).parent / "fixtures"


def test_mymotif_imports_csv_motif_table():
    client = create_app().test_client()

    response = client.post(
        "/api/mymotif/parse",
        data={
            "motif_files": (
                BytesIO(b'rec_seq,enz_type,meth_base\n"GATC",2,A\n'),
                "motifs.csv",
            )
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    assert response.json["file_errors"] == []
    assert response.json["motifs"] == [{
        "rec_seq": "GATC",
        "enz_type": "2",
        "meth_base": "A",
        "meth_type": "-",
        "comp_meth_base": "-",
        "comp_meth_type": "-",
    }]


def test_mymotif_imports_compact_rebase_file_without_trailing_delimiter():
    client = create_app().test_client()

    response = client.post(
        "/api/mymotif/parse",
        data={
            "motif_files": (
                BytesIO(b"<enz_type>2<rec_seq>GATCNAC<meth_base>6<meth_type>m6A"),
                "strain.rebase",
            )
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    assert response.json["motifs"][0]["rec_seq"] == "GATCNAC"
    assert response.json["motifs"][0]["enz_type"] == "2"
    assert response.json["motifs"][0]["meth_base"] == "6"
    assert response.json["motifs"][0]["meth_type"] == "m6A"
    assert response.json["motifs"][0]["comp_meth_base"] == "-"
    assert response.json["motifs"][0]["comp_meth_type"] == "-"


def test_mymotif_imports_mijamp_expected_output_file():
    client = create_app().test_client()

    with open(FIXTURES / "mymotif_mijamp_expectedOutput.tsv", "rb") as motif_file:
        response = client.post(
            "/api/mymotif/parse",
            data={
                "motif_files": (
                    motif_file,
                    "mymotif_mijamp_expectedOutput.tsv",
                )
            },
            content_type="multipart/form-data",
        )

    assert response.status_code == 200
    assert [motif["rec_seq"] for motif in response.json["motifs"]] == [
        "GATC",
        "AACNNNNNNGTGC",
        "GCACNNNNNNGTT",
        "CCWGG",
    ]
    first = response.json["motifs"][0]
    assert first["enz_type"] == "2"
    assert first["meth_base"] == "2"
    assert first["meth_type"] == "m6A"
    assert first["comp_meth_base"] == "3"
    assert first["comp_meth_type"] == "m6A"
