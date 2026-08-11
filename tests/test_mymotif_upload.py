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
        "meth_base": "2",
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


def test_mymotif_imports_mijamp_token_with_m6a_marker():
    client = create_app().test_client()

    response = client.post(
        "/api/mymotif/parse",
        data={
            "motif_files": (
                BytesIO(b"G(m6A)TC\n"),
                "motifs.tsv",
            )
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    first = response.json["motifs"][0]
    assert first["rec_seq"] == "GATC"
    assert first["enz_type"] == "2"
    assert first["meth_base"] == "2"
    assert first["meth_type"] == "m6A"
    assert first["comp_meth_base"] == "3"
    assert first["comp_meth_type"] == "m6A"


def test_mymotif_imports_mijamp_block_header_methylation_type():
    client = create_app().test_client()

    content = (
        b"#6mA modified motifs\n"
        b"#motif\tmethylCounts\ttotalCounts\t%motifsMethyl\n"
        b"G(6mA)TC\t1\t1\t100\n"
        b"A(6mA)CNNNNNNGTGC\t1\t1\t100\n"
        b"#5mC modified motifs\n"
        b"#motif\tmethylCounts\ttotalCounts\t%motifsMethyl\n"
        b"C(5mC)WGG\t1\t1\t100\n"
    )

    response = client.post(
        "/api/mymotif/parse",
        data={
            "motif_files": (
                BytesIO(content),
                "mijamp.tsv",
            )
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    motifs = response.json["motifs"]
    assert [motif["rec_seq"] for motif in motifs] == ["GATC", "AACNNNNNNGTGC", "CCWGG"]
    assert [motif["meth_type"] for motif in motifs] == ["m6A", "m6A", "m5C"]
    assert [motif["meth_base"] for motif in motifs] == ["2", "2", "2"]


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


def test_mymotif_imports_meth_base_as_position_and_unknown_type_as_unk():
    client = create_app().test_client()

    response = client.post(
        "/api/mymotif/parse",
        data={
            "motif_files": (
                BytesIO(b"<enz_type>2<rec_seq>ATGC<meth_base>C<meth_type>99"),
                "motifs.rebase",
            )
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    first = response.json["motifs"][0]
    assert first["meth_base"] == "4"
    assert first["meth_type"] == "Unk"
