"""
visualization.py
=================
Interactive plasmid-map visualization using Plotly, replacing the earlier
dna_features_viewer/matplotlib approach (static PNG/SVG/PDF that got too
busy with labeling once a construct had more than a handful of features).

Two figures are produced:
    - "before": every motif SyToGen originally found, plotted against the
      unedited construct.
    - "after": the same motif footprint, but styled by what happened to
      each one (resolved / unresolved / newly introduced), with hover
      text pulled directly from that motif's decision-matrix row — so
      mousing over a marker shows exactly why SyToGen made the choice it
      did, without needing to cross-reference decision_matrix.tsv by hand.

Layout (both circular and linear):
    - Outermost ring/row: genes (CDS/ORF/Marker) and protected/regulatory
      sites, each their own color, in two adjacent sub-bands that
      together make up one bordered "genes" ring — genes labeled directly
      on the plot, protected sites hover-only (there are usually more of
      them, and they're typically short).
    - Inner rings/rows: one per distinct motif pattern, each its own
      color, nested progressively inward (circular) or stacked below
      (linear). Every ring/row — gene ring and each motif track — is
      bounded by its own black border and separated from its neighbors
      by a visible gap, so nothing overlaps or bleeds into the next ring.
    - Background is gray, distinct from any marker/ring color, so
      borders and light-colored markers both stay visible against it.

Both figures are returned as plotly.graph_objects.Figure objects. Callers
can embed them live on a page via fig.to_json() + Plotly.js, or export a
self-contained interactive HTML file via fig.to_html() — no image-export
dependency (e.g. kaleido) is needed for either.
"""

import plotly.graph_objects as go
import plotly.colors


GENE_TYPES = {"CDS", "ORF", "Marker"}
PROTECTED_TYPES = {"regulatory", "misc_feature", "rep_origin", "promoter", "RBS"}
MAX_PROTECTED_LENGTH = 100  # bp — mirrors sytogen_runner._parse_protected_regions

GENE_COLOR = "#7fa8c9"
PROTECTED_COLOR = "#c98f7f"
BORDER_COLOR = "black"
BORDER_WIDTH = 1.5
BACKGROUND_COLOR = "#ececec"
DONUT_CENTER_COLOR = "#1c1c1c"

RESOLVED_OPACITY = 0.30
NEW_MOTIF_COLOR = "#e63946"
NEW_MOTIF_SYMBOL = "x"


# =============================================================================
# Shared data prep
# =============================================================================

GENES_COLOR_KEY = "__GENES__"


def _assign_colors(motif_patterns):
    """
    One color per distinct motif pattern PLUS genes, sampled evenly across
    the full viridis colorscale — genes always take the first color (the
    start of the scale), then each motif pattern gets the next one in
    order, so genes and motifs share one consistent, ordered palette
    instead of genes being an unrelated fixed color.
    """
    ordered = [GENES_COLOR_KEY] + sorted(motif_patterns)
    n = len(ordered)
    if n == 1:
        sample_points = [0.0]
    else:
        sample_points = [i / (n - 1) for i in range(n)]
    colors = plotly.colors.sample_colorscale("Viridis", sample_points)
    return dict(zip(ordered, colors))


def _is_motiffinder_hit_marker(feature):
    """
    MotifFinder writes each motif it finds back into the GenBank as its own
    misc_feature entry (qualifiers include 'ID'='motif_hit_NNNN', 'motif',
    'hit_seq'). These mark motif sites, not protected regulatory regions —
    same exclusion sytogen_runner._parse_protected_regions applies, kept in
    sync here so the gene-ring protected-site band matches what the actual
    editing pipeline treated as protected.
    """
    qualifiers = feature.qualifiers
    feature_id = str(qualifiers.get("ID", [""])[0])
    return (
        "motif" in qualifiers
        or "hit_seq" in qualifiers
        or feature_id.startswith("motif_hit_")
    )


def _extract_spans(record, feature_types, length_cutoff=None):
    """
    Shared extraction for genes and protected regions alike — pulls
    (id, start, end, strand) for every feature of the given type(s),
    handling an origin-spanning CompoundLocation the same way
    sytogen_runner._parse_genes / _parse_protected_regions do. BioPython's
    naive min/max span for a wrapped join() location covers nearly the
    whole molecule, not just the feature's real footprint.
    """
    length = len(record.seq)
    spans = []
    for i, feature in enumerate(record.features):
        if feature.type not in feature_types:
            continue
        if feature_types is PROTECTED_TYPES and _is_motiffinder_hit_marker(feature):
            continue

        parts = getattr(feature.location, "parts", [feature.location])
        naive_start = int(feature.location.start)
        naive_end = int(feature.location.end) - 1
        declared_length = sum(len(part) for part in parts)

        if len(parts) > 1 and (naive_end - naive_start + 1) > declared_length:
            start = int(parts[0].start)
            end = length + int(parts[-1].end) - 1  # raw, may exceed length
        else:
            start = naive_start
            end = naive_end

        if length_cutoff is not None and (end - start + 1) > length_cutoff:
            continue

        strand = "+" if feature.location.strand >= 0 else "-"
        default_label = f"gene_{i + 1}" if feature_types is GENE_TYPES else f"{feature.type}_{i + 1}"
        label = (
            feature.qualifiers.get("gene", [None])[0]
            or feature.qualifiers.get("locus_tag", [None])[0]
            or feature.qualifiers.get("sequence", [None])[0]
            or default_label
        )
        spans.append({"id": label, "start": start, "end": end, "strand": strand})
    return spans


def _extract_genes(record):
    return _extract_spans(record, GENE_TYPES)


def _extract_protected_regions(record):
    return _extract_spans(record, PROTECTED_TYPES, length_cutoff=MAX_PROTECTED_LENGTH)


# =============================================================================
# MotifFinder-native data (dicts, not Bio.SeqRecord/SeqFeature)
# =============================================================================
# MotifFinder works from its own plain dict formats (motiffinder_backend.py):
# GFF3-style feature dicts (seqid/type/start/end/strand/attrs, 1-based) and
# hit dicts (pos_0/strand/rec_seq/hit_seq/enz_type/..., 0-based). These
# extractors mirror _extract_spans/_extract_genes/_extract_protected_regions
# above but read that shape directly, so build_motiffinder_map() doesn't
# need a Bio.SeqRecord at all — MotifFinder never builds one for the
# FASTA+GFF3 input path.

def _is_motiffinder_hit_marker_dict(feature):
    attrs = feature.get("attrs", "")
    feature_id = ""
    for part in attrs.split(";"):
        part = part.strip()
        if part.lower().startswith("id="):
            feature_id = part.split("=", 1)[1]
    attrs_lower = attrs.lower()
    return (
        "motif=" in attrs_lower
        or "hit_seq=" in attrs_lower
        or feature_id.startswith("motif_hit_")
    )


def _label_from_gff3_attrs(attrs, fallback):
    for key in ("Name=", "gene=", "ID="):
        for part in attrs.split(";"):
            part = part.strip()
            if part.startswith(key):
                return part.split("=", 1)[1]
    return fallback


def _extract_spans_from_gff3_dicts(features, feature_types, length_cutoff=None):
    spans = []
    for i, f in enumerate(features):
        if f.get("type") not in feature_types:
            continue
        if feature_types is PROTECTED_TYPES and _is_motiffinder_hit_marker_dict(f):
            continue

        start = int(f["start"]) - 1  # GFF3 1-based inclusive -> 0-based
        end = int(f["end"]) - 1
        if length_cutoff is not None and (end - start + 1) > length_cutoff:
            continue

        strand = f.get("strand", "+")
        if strand not in ("+", "-"):
            strand = "+"
        default_label = f"gene_{i + 1}" if feature_types is GENE_TYPES else f"{f.get('type')}_{i + 1}"
        label = _label_from_gff3_attrs(f.get("attrs", ""), default_label)
        spans.append({"id": label, "start": start, "end": end, "strand": strand})
    return spans


def _hit_tracks(hits):
    """Group MotifFinder's hit dicts into one track per distinct recognition
    pattern, plus the gene color to use for this same run (colors are
    assigned together so genes get the first viridis color)."""
    patterns = sorted({h["rec_seq"] for h in hits})
    colors = _assign_colors(patterns)
    gene_color = colors[GENES_COLOR_KEY]

    tracks = []
    for pattern in patterns:
        points = []
        for h in hits:
            if h["rec_seq"] != pattern:
                continue
            end = h["pos_0"] + len(pattern) - 1
            hover = f"{pattern} ({h['strand']} strand)<br>position {h['pos_0']}-{end}"
            if h.get("enz_type"):
                hover += f"<br>enzyme type {h['enz_type']}"
            points.append({"position": h["pos_0"], "hover": hover})
        tracks.append({"label": pattern, "color": colors[pattern], "points": points})
    return tracks, gene_color


def build_motiffinder_map(features, hits, sequence_length, topology, title=""):
    """
    MotifFinder's own map — every motif hit it found, one track per
    distinct pattern, against the construct's genes and protected sites.
    There's no before/after distinction here (MotifFinder never edits
    anything); this is what SyToGen's "before" plot is built from the
    same way, just fed from MotifFinder's own hit-search results directly
    instead of a decision matrix.
    """
    genes = _extract_spans_from_gff3_dicts(features, GENE_TYPES)
    protected_regions = _extract_spans_from_gff3_dicts(features, PROTECTED_TYPES, length_cutoff=MAX_PROTECTED_LENGTH)
    tracks, gene_color = _hit_tracks(hits)

    builder = _build_circular_figure if topology == "circular" else _build_linear_figure
    return builder(genes, protected_regions, tracks, sequence_length, title or "Motif map", gene_color)


# =============================================================================
# SyToGen-native data (Motif objects, decision matrix)
# =============================================================================

def _reasoning_lookup(decision_matrix):
    """
    {(motif, start, end, strand): reasoning} for the row that actually
    represents each motif's outcome — the chosen candidate if it
    resolved, or the sentinel row if it didn't. Rejected-but-not-chosen
    alternative candidates are skipped; they're not what happened to the
    motif, just other options that were considered.
    """
    groups = {}
    for row in decision_matrix:
        key = (row["motif"], row["motif_start"], row["motif_end"], row["motif_strand"])
        groups.setdefault(key, []).append(row)

    lookup = {}
    for key, rows in groups.items():
        chosen = next((r for r in rows if r.get("chosen")), None)
        representative = chosen or rows[0]
        lookup[key] = representative.get("reasoning", "")
    return lookup


def _motif_status(motif, resolved_motif_keys, reasoning_lookup):
    """
    'resolved' | 'unresolved', from the authoritative resolved_motif_keys
    set built during the actual edit loop — not inferred from whether a
    decision-matrix row happens to exist, since a motif resolved as a
    side effect of a DIFFERENT edit never gets a row of its own.
    """
    key = (motif.motif, motif.start, motif.end, motif.strand)
    status = "resolved" if key in resolved_motif_keys else "unresolved"
    reasoning = reasoning_lookup.get(
        key,
        "Resolved as a side effect of a different nearby edit."
        if status == "resolved" else
        "No decision-matrix entry found for this occurrence."
    )
    return status, reasoning


# =============================================================================
# Ring/row layout
# =============================================================================

def _compute_circular_bands(n_motif_tracks, gap=0.015):
    """
    Non-overlapping radius bands, outside to inside: the combined gene
    ring (split into a genes sub-band and a protected-sites sub-band,
    distinguished by color only — no divider line between them), then
    one band per motif track — each separated from its neighbors by a
    small `gap` so borders don't touch, but kept tight rather than
    spread thin.
    """
    gene_outer, gene_split, gene_inner = 0.98, 0.90, 0.82
    track_region_outer = gene_inner - gap
    track_region_inner = 0.20

    tracks = []
    if n_motif_tracks:
        step = (track_region_outer - track_region_inner) / n_motif_tracks
        for i in range(n_motif_tracks):
            outer = track_region_outer - i * step
            inner = outer - (step - gap)
            tracks.append((outer, max(inner, track_region_inner)))

    return {
        "gene_outer": gene_outer, "gene_split": gene_split, "gene_inner": gene_inner,
        "tracks": tracks, "center_radius": track_region_inner,
    }


def _compute_linear_rows(n_motif_tracks):
    """Row *centers*, top to bottom: genes, protected sites, then one per motif track."""
    total_rows = 2 + n_motif_tracks
    return {
        "gene_row": total_rows,
        "protected_row": total_rows - 1,
        "track_rows": [total_rows - 2 - i for i in range(n_motif_tracks)],
        "total_rows": total_rows,
    }


# =============================================================================
# Circular figure
# =============================================================================

def _arc_theta_width(span, length):
    """(center_theta_degrees, width_degrees), wrap-aware."""
    full_span = span["end"] - span["start"] + 1
    center = (span["start"] + full_span / 2.0) % length
    return (center / length) * 360.0, (full_span / length) * 360.0


def _add_circular_border(fig, radius):
    thetas = list(range(0, 361, 2))
    fig.add_trace(go.Scatterpolar(
        r=[radius] * len(thetas), theta=thetas, mode="lines",
        line=dict(color=BORDER_COLOR, width=BORDER_WIDTH),
        hoverinfo="skip", showlegend=False,
    ))


def _add_circular_center_fill(fig, radius, color):
    """Solid dark disk filling the donut hole (r=0 to radius) — the
    'inside' the innermost motif ring, rather than blending into the
    same gray as the rest of the plot's background."""
    thetas = list(range(0, 361, 4))
    fig.add_trace(go.Scatterpolar(
        r=[radius] * len(thetas), theta=thetas, mode="lines",
        fill="toself", fillcolor=color,
        line=dict(color=color, width=0),
        hoverinfo="skip", showlegend=False,
    ))


def _add_arc_band(fig, spans, length, outer, inner, color, name, hover_fn, label_fn=None):
    if not spans:
        return
    thetas = [_arc_theta_width(s, length)[0] for s in spans]
    widths = [_arc_theta_width(s, length)[1] for s in spans]
    fig.add_trace(go.Barpolar(
        r=[outer - inner] * len(spans), theta=thetas, width=widths, base=inner,
        marker_color=color, marker_line_width=0, name=name,
        hovertext=[hover_fn(s) for s in spans], hoverinfo="text", opacity=0.95,
    ))
    if label_fn:
        mid_r = (outer + inner) / 2.0
        fig.add_trace(go.Scatterpolar(
            r=[mid_r] * len(spans), theta=thetas, mode="text",
            text=[label_fn(s) for s in spans], textfont=dict(size=9, color="black"),
            hoverinfo="skip", showlegend=False,
        ))


def _build_circular_figure(genes, protected_regions, motif_tracks, length, title, gene_color=GENE_COLOR):
    fig = go.Figure()
    bands = _compute_circular_bands(len(motif_tracks))

    _add_arc_band(
        fig, genes, length, bands["gene_outer"], bands["gene_split"],
        gene_color, "Genes",
        hover_fn=lambda g: f"{g['id']} ({g['strand']} strand)<br>{g['start']}-{g['end'] % length}",
        label_fn=lambda g: g["id"],
    )
    _add_arc_band(
        fig, protected_regions, length, bands["gene_split"], bands["gene_inner"],
        PROTECTED_COLOR, "Protected sites",
        hover_fn=lambda p: f"{p['id']} (protected)<br>{p['start']}-{p['end'] % length}",
    )
    # Single line at the outer edge and single line at the inner edge of
    # the combined gene+protected ring — no divider line at gene_split;
    # the color change alone is enough to tell genes from protected sites.
    _add_circular_border(fig, bands["gene_outer"])
    _add_circular_border(fig, bands["gene_inner"])

    for track, (outer, inner) in zip(motif_tracks, bands["tracks"]):
        radius = (outer + inner) / 2.0
        thetas = [(p["position"] % length) / length * 360.0 for p in track["points"]]
        fig.add_trace(go.Scatterpolar(
            r=[radius] * len(thetas), theta=thetas, mode="markers",
            marker=dict(
                color=track.get("point_color", track["color"]),
                size=track.get("size", 9), symbol=track.get("symbol", "circle"),
                opacity=track.get("opacity", 1.0), line=dict(width=1, color="white"),
            ),
            name=track["label"], hovertext=[p["hover"] for p in track["points"]],
            hoverinfo="text",
        ))
        _add_circular_border(fig, outer)

    # One shared line at the innermost boundary (the outermost track already
    # drew its own outer edge above; each subsequent track only needs its
    # own outer edge too, since that IS the previous track's inner edge —
    # avoids drawing two lines on top of each other at every internal gap).
    if motif_tracks:
        _add_circular_border(fig, bands["tracks"][-1][1])

    # Dark "donut hole" center, rather than blending into the background.
    _add_circular_center_fill(fig, bands["center_radius"], DONUT_CENTER_COLOR)

    fig.update_layout(
        title=title,
        paper_bgcolor=BACKGROUND_COLOR,
        polar=dict(
            bgcolor=BACKGROUND_COLOR,
            radialaxis=dict(visible=False, range=[0, 1.05]),
            angularaxis=dict(
                rotation=90, direction="clockwise",
                tickmode="array",
                tickvals=[0, 90, 180, 270],
                ticktext=["0", f"{length//4}", f"{length//2}", f"{3*length//4}"],
            ),
        ),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=-0.15),
        margin=dict(l=40, r=40, t=60, b=40),
    )
    return fig


# =============================================================================
# Linear figure
# =============================================================================

def _add_row_border(fig, x0, x1, y, half_height):
    fig.add_shape(
        type="rect", x0=x0, x1=x1,
        y0=y - half_height, y1=y + half_height,
        line=dict(color=BORDER_COLOR, width=BORDER_WIDTH), fillcolor="rgba(0,0,0,0)",
    )


def _build_linear_figure(genes, protected_regions, motif_tracks, length, title, gene_color=GENE_COLOR):
    fig = go.Figure()
    rows = _compute_linear_rows(len(motif_tracks))
    row_height = 0.8

    def _add_span_row(spans, y, color, name, hover_fn, label_fn=None):
        for s in spans:
            end = min(s["end"], length - 1)
            fig.add_shape(
                type="rect", x0=s["start"], x1=end,
                y0=y - row_height / 2, y1=y + row_height / 2,
                fillcolor=color, line=dict(width=0), opacity=0.95,
            )
            fig.add_trace(go.Scatter(
                x=[(s["start"] + end) / 2], y=[y], mode="markers",
                marker=dict(opacity=0), hovertext=[hover_fn(s)], hoverinfo="text",
                showlegend=False,
            ))
            if label_fn:
                fig.add_annotation(
                    x=(s["start"] + end) / 2, y=y, text=label_fn(s),
                    showarrow=False, font=dict(size=9, color="black"),
                )
        _add_row_border(fig, 0, length, y, row_height / 2)
        if spans:
            # dummy trace so the row shows up in the legend with its color
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers",
                marker=dict(color=color, size=10), name=name,
            ))

    _add_span_row(genes, rows["gene_row"], gene_color, "Genes",
                  hover_fn=lambda g: f"{g['id']} ({g['strand']} strand)<br>{g['start']}-{min(g['end'], length-1)}",
                  label_fn=lambda g: g["id"])
    _add_span_row(protected_regions, rows["protected_row"], PROTECTED_COLOR, "Protected sites",
                  hover_fn=lambda p: f"{p['id']} (protected)<br>{p['start']}-{min(p['end'], length-1)}")

    for track, y in zip(motif_tracks, rows["track_rows"]):
        xs = [p["position"] for p in track["points"]]
        fig.add_trace(go.Scatter(
            x=xs, y=[y] * len(xs), mode="markers",
            marker=dict(
                color=track.get("point_color", track["color"]),
                size=track.get("size", 10), symbol=track.get("symbol", "circle"),
                opacity=track.get("opacity", 1.0), line=dict(width=1, color="white"),
            ),
            name=track["label"], hovertext=[p["hover"] for p in track["points"]],
            hoverinfo="text",
        ))
        _add_row_border(fig, 0, length, y, row_height / 2)

    fig.update_layout(
        title=title,
        paper_bgcolor=BACKGROUND_COLOR,
        plot_bgcolor=BACKGROUND_COLOR,
        xaxis=dict(title="Position (bp)", range=[0, length]),
        yaxis=dict(visible=False, range=[0, rows["total_rows"] + 1]),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=-0.25),
        margin=dict(l=40, r=40, t=60, b=40),
    )
    return fig


# =============================================================================
# Public entry point
# =============================================================================

def build_plasmid_maps(output_record, motifs, new_motifs, decision_matrix,
                        resolved_motif_keys, sequence_length, topology, title=""):
    """
    Returns (fig_before, fig_after) — both plotly.graph_objects.Figure.

    fig_before: every originally-found motif occurrence, one track per
    distinct pattern, plotted against the unedited construct's genes and
    protected sites.

    fig_after: the same footprint, but each marker styled by outcome —
    solid = still present (unresolved), faded/hollow = resolved (gone),
    plus an 'x' marker in a shared alert color for anything
    _final_new_motif_check found that wasn't there originally. Hover text
    on fig_after is the actual decision-matrix reasoning for that motif
    (resolved_motif_keys is the authoritative status — it's built during
    the real edit loop, since a motif resolved as a side effect of a
    different edit never gets a decision-matrix row of its own to infer
    status from).

    Genes and protected sites are read directly from output_record —
    output_record is a deepcopy of the original input record with edit
    markers added on top, so its original annotations (CDS/ORF/Marker,
    regulatory/misc_feature/rep_origin/promoter/RBS) are all still there.
    """
    genes = _extract_genes(output_record)
    protected_regions = _extract_protected_regions(output_record)
    patterns = {m.motif for m in motifs} | {nm["motif"] for nm in new_motifs}
    colors = _assign_colors(patterns)
    gene_color = colors[GENES_COLOR_KEY]
    reasoning_lookup = _reasoning_lookup(decision_matrix)

    builder = _build_circular_figure if topology == "circular" else _build_linear_figure

    # ---- before ----
    before_tracks = []
    for pattern in sorted({m.motif for m in motifs}):
        points = []
        for m in motifs:
            if m.motif != pattern:
                continue
            points.append({
                "position": m.start,
                "hover": f"{pattern} ({m.strand} strand)<br>position {m.start}-{m.end}",
            })
        before_tracks.append({"label": pattern, "color": colors[pattern], "points": points})

    fig_before = builder(genes, protected_regions, before_tracks, sequence_length,
                          title or "Motifs before SyToGen", gene_color)

    # ---- after ----
    after_tracks = []
    for pattern in sorted({m.motif for m in motifs}):
        resolved_points, unresolved_points = [], []
        for m in motifs:
            if m.motif != pattern:
                continue
            status, reasoning = _motif_status(m, resolved_motif_keys, reasoning_lookup)
            hover = (
                f"{pattern} ({m.strand} strand)<br>position {m.start}-{m.end}<br>"
                f"<b>{status.upper()}</b><br>{reasoning}"
            )
            point = {"position": m.start, "hover": hover}
            if status == "resolved":
                resolved_points.append(point)
            else:
                unresolved_points.append(point)

        if unresolved_points:
            after_tracks.append({
                "label": f"{pattern} (unresolved)", "color": colors[pattern],
                "points": unresolved_points, "opacity": 1.0,
            })
        if resolved_points:
            after_tracks.append({
                "label": f"{pattern} (resolved)", "color": colors[pattern],
                "point_color": colors[pattern], "points": resolved_points,
                "opacity": RESOLVED_OPACITY, "symbol": "circle-open",
            })

    if new_motifs:
        new_points = [{
            "position": nm["start"],
            "hover": (
                f"{nm['motif']} ({nm['strand']} strand)<br>position {nm['start']}-{nm['end']}<br>"
                f"<b>NEWLY INTRODUCED</b><br>Not present in the original construct — "
                f"introduced as a side effect of editing elsewhere."
            ),
        } for nm in new_motifs]
        after_tracks.append({
            "label": "newly introduced", "color": NEW_MOTIF_COLOR,
            "point_color": NEW_MOTIF_COLOR, "points": new_points,
            "symbol": NEW_MOTIF_SYMBOL, "size": 12,
        })

    fig_after = builder(genes, protected_regions, after_tracks, sequence_length,
                         title or "Motifs after SyToGen", gene_color)

    return fig_before, fig_after
