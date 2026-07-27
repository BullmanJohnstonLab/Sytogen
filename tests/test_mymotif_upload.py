from io import BytesIO

from sytogen import create_app


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
