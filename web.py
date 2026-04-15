from flask import Blueprint, render_template

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


