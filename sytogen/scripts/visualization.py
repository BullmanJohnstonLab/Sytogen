"""
visualization.py
=================
Recreates the plasmid-map visualization legacy_sytogen.py used to produce
(circular + linear maps via dna_features_viewer/matplotlib, plus a combined
PDF), adapted to the rewrite's data structures.

Legacy wrote these to disk per-run as:
    {id}_circ_tool_repr.png / .svg
    {id}_lin_tool_repr.png  / .svg
    {id}_combined_plot_tool_repr.pdf

render_plasmid_maps() below produces the same set, in memory, as a
{filename: bytes} dict so the caller (api.py) can drop them straight into
the output zip instead of writing to a shared output folder.

Same visual language as legacy:
    - genes (CDS/ORF/Marker) all one color, with labels
    - each distinct restriction-motif recognition sequence gets its own
      color and a legend entry (legacy's "RMM-i" per-pattern features)
    - the edits SyToGen actually applied are overlaid as their own color
      (legacy didn't have this — SyToGen's rewrite tracks edits as
      first-class 'SyT' features on the output GenBank, so we surface them
      here too)

Requires 'dna_features_viewer' and 'matplotlib' — the exact same two
dependencies legacy_sytogen.py already required for this; nothing new is
being added to the dependency footprint.
"""

import io
import copy
import random

import matplotlib
matplotlib.use("Agg")  # headless rendering — must happen before pyplot import
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.backends.backend_pdf

from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


_BASE_PALETTE = [
    "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99",
]
_GENE_COLOR = "#a6cee3"   # same blue legacy used for CDS
_EDIT_COLOR = "#d62728"
_GENE_TYPES = {"CDS", "ORF", "Marker"}


def _random_hex():
    return "#" + "".join(random.choice("0123456789ABCDEF") for _ in range(6))


def _build_annotated_record(output_record, motifs):
    """
    output_record already carries the gene features (CDS/ORF/Marker) from
    the input GenBank plus the 'SyT' edit-marker features api.py adds for
    every applied mutation. Layer one feature per distinct motif *pattern*
    on top (mirroring legacy's per-pattern "RMM-i" features) so each
    recognition-site type gets its own color and legend entry.

    A motif that spans the circular origin (end >= length — see
    genome_model.Motif / sytogen_runner._parse_motifs) is split into two
    features, the same way legacy split any "join(...)" location before
    plotting it.
    """
    record = copy.deepcopy(output_record)
    length = len(record.seq)

    motif_feature_types = {}
    for i, pattern in enumerate(sorted({m.motif for m in motifs})):
        motif_feature_types[pattern] = f"Motif-{i}"

    for m in motifs:
        ftype = motif_feature_types[m.motif]
        strand = 1 if m.strand == "+" else -1
        start = m.start % length
        end = m.end % length

        if m.end < length:
            record.features.append(SeqFeature(
                FeatureLocation(start, end + 1, strand=strand),
                type=ftype, qualifiers={"note": [m.motif]},
            ))
        else:
            # origin-spanning — two pieces, same as legacy's join() split
            record.features.append(SeqFeature(
                FeatureLocation(start, length, strand=strand),
                type=ftype, qualifiers={"note": [m.motif]},
            ))
            record.features.append(SeqFeature(
                FeatureLocation(0, end + 1, strand=strand),
                type=ftype, qualifiers={"note": [m.motif]},
            ))

    return record, motif_feature_types


def _build_color_map(record, motif_feature_types):
    feature_types = sorted({f.type for f in record.features})
    palette = list(_BASE_PALETTE)
    random.shuffle(palette)

    color_map = {}
    for ftype in feature_types:
        if ftype == "SyT":
            color_map[ftype] = _EDIT_COLOR
        elif ftype in _GENE_TYPES:
            color_map[ftype] = _GENE_COLOR
        elif ftype in motif_feature_types.values():
            color_map[ftype] = palette.pop() if palette else _random_hex()
        else:
            color_map[ftype] = "#000000"
    return color_map


def render_plasmid_maps(output_record, motifs, title="", figure_width=9):
    """
    Parameters
    ----------
    output_record : Bio.SeqRecord
        The final annotated construct (same object api.py already builds
        for sytogen_result.gbk — gene features + 'SyT' edit markers).
    motifs : list[Motif]
        The input motif list from run_sytogen_pipeline's result dict.
    title : str
        Plot title (e.g. the sequence id).
    figure_width : int
        Same 'image_width' knob legacy exposed via --image_width (default 9).

    Returns
    -------
    dict[str, bytes] — filenames mapped to file contents, ready to drop
    into a zip. Any piece that fails to render (missing/incompatible
    features, a dna_features_viewer edge case, etc.) is silently skipped
    rather than failing the whole SyToGen run, mirroring legacy's own
    try/except-around-each-plot behavior.
    """
    annotated_record, motif_feature_types = _build_annotated_record(output_record, motifs)
    color_map = _build_color_map(annotated_record, motif_feature_types)
    label_types = _GENE_TYPES | set(motif_feature_types.values())

    class _Translator(BiopythonTranslator):
        def compute_feature_color(self, feature):
            return color_map.get(feature.type, "#888888")

        def compute_feature_label(self, feature):
            if feature.type not in label_types:
                return None
            return BiopythonTranslator.compute_feature_label(self, feature)

    def _legend_patches():
        patches = [mpatches.Patch(color=_GENE_COLOR, label="Gene (CDS/ORF/Marker)")]
        for pattern, ftype in motif_feature_types.items():
            patches.append(mpatches.Patch(color=color_map[ftype], label=pattern))
        if "SyT" in color_map:
            patches.append(mpatches.Patch(color=_EDIT_COLOR, label="Edit applied"))
        return patches

    outputs = {}
    ax1 = ax2 = None

    # --- circular map ---
    try:
        graphic_record = _Translator().translate_record(
            annotated_record, record_class=CircularGraphicRecord
        )
        ax1, _ = graphic_record.plot(figure_width=figure_width)
        ax1.legend(handles=_legend_patches(), loc="upper left", bbox_to_anchor=(1.05, 1.0))
        if title:
            ax1.set_title(title, fontsize=14)

        png_buf = io.BytesIO()
        ax1.figure.savefig(png_buf, bbox_inches="tight", dpi=400, format="png")
        outputs["plasmid_map_circular.png"] = png_buf.getvalue()

        svg_buf = io.BytesIO()
        ax1.figure.savefig(svg_buf, bbox_inches="tight", format="svg")
        outputs["plasmid_map_circular.svg"] = svg_buf.getvalue()
    except Exception:
        ax1 = None

    # --- linear map ---
    try:
        graphic_record = _Translator().translate_record(annotated_record)
        ax2, _ = graphic_record.plot(figure_width=figure_width)
        ax2.legend(handles=_legend_patches(), loc="upper left", bbox_to_anchor=(1.05, 1.0))
        if title:
            ax2.set_title(title, fontsize=14)

        png_buf = io.BytesIO()
        ax2.figure.savefig(png_buf, bbox_inches="tight", dpi=400, format="png")
        outputs["plasmid_map_linear.png"] = png_buf.getvalue()

        svg_buf = io.BytesIO()
        ax2.figure.savefig(svg_buf, bbox_inches="tight", format="svg")
        outputs["plasmid_map_linear.svg"] = svg_buf.getvalue()
    except Exception:
        ax2 = None

    # --- combined PDF, only if both maps rendered ---
    if ax1 is not None and ax2 is not None:
        try:
            pdf_buf = io.BytesIO()
            pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_buf)
            pdf.savefig(ax1.figure, bbox_inches="tight")
            pdf.savefig(ax2.figure, bbox_inches="tight")
            pdf.close()
            outputs["plasmid_map_combined.pdf"] = pdf_buf.getvalue()
        except Exception:
            pass

    plt.close("all")
    return outputs
