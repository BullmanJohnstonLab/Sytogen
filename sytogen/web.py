from flask import Blueprint, render_template
from flask import request, send_file, jsonify
import tempfile
import os
from sytogen.scripts.codon_bias_estimator import run_codon_bias

web = Blueprint("web", __name__)

# Main page and tool forms

@web.route("/")
def index():
    return render_template("index.html")

@web.route("/mymotif")
def mymotif_form():
    return render_template("mymotif.html")

@web.route("/motiffinder")
def motiffinder_form():
    return render_template("motiffinder.html")

@web.route("/codon-bias")
def codon_bias_form():
    return render_template("codonbias.html")

@web.route("/sytogen")
def sytogen_form():
    return render_template("sytogen.html")

# Additional informational pages

@web.route("/explained")
def explained():
    return render_template("explained.html")

@web.route("/what-is")
def whatis():
    return render_template("whatis.html")

@web.route("/user-guide")
def user_guide():
    return render_template("user_guide.html")

@web.route("/codonbias/run", methods=["POST"])
def run_codonbias():
    if "genome_file" not in request.files:
        return jsonify(error="Missing genome_file"), 400

    genome_file = request.files["genome_file"]
    codon_table = request.form.get("codon_table", "11")

    with tempfile.TemporaryDirectory() as tmpdir:
        genome_path = os.path.join(tmpdir, genome_file.filename)
        genome_file.save(genome_path)

        zip_path = run_codon_bias(
            genome_path=genome_path,
            codon_table=int(codon_table),
            output_dir=tmpdir,
        )

        return send_file(
            zip_path,
            as_attachment=True,
            download_name="codonbias_output.zip",
        )
