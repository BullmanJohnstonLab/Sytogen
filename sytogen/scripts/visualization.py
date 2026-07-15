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
    - Outer ring/line: genes (CDS/ORF/Marker), one uniform color.
    - Inner rings/lines: one per distinct motif pattern, each its own
      color, nested progressively inward (circular) or stacked below
      (linear).

Both figures are returned as plotly.graph_objects.Figure objects. Callers
can embed them live on a page via fig.to_json() + Plotly.js, or export a
self-contained interactive HTML file via fig.to_html() — no image-export
dependency (e.g. kaleido) is needed for either.
"""

import plotly.graph_objects as go


GENE_TYPES = {"CDS", "ORF", "Marker"}
GENE_COLOR = "#7fa8c9"
GENE_ROW_LABEL = "Genes"

RESOLVED_OPACITY = 0.30
NEW_MOTIF_COLOR = "#e63946"
NEW_MOTIF_SYMBOL_CIRCULAR = "x"
NEW_MOTIF_SYMBOL_LINEAR = "x"

_PALETTE = [
    "#1f78b4", "#33a02c", "#ff7f00", "#6a3d9a", "#e31a1c",
    "#b15928", "#a6cee3", "#b2df8a", "#fdbf6f", "#cab2d6",
    "#ffff99", "#fb9a99",
]


# =============================================================================
# Shared data prep
# =============================================================================

def _assign_pattern_colors(patterns):
    """One color per distinct motif pattern, deterministic (sorted) order."""
    ordered = sorted(patterns)
    return {p: _PALETTE[i % len(_PALETTE)] for i, p in enumerate(ordered)}


def _extract_genes(record):
    """
    Pull (id, start, end, strand) for every gene feature, handling an
    origin-spanning CompoundLocation the same way sytogen_runner._parse_genes
    does — BioPython's naive min/max span for a wrapped join() location
    covers nearly the whole molecule, not just the gene's real footprint.
    """
    length = len(record.seq)
    genes = []
    for i, feature in enumerate(record.features):
        if feature.type not in GENE_TYPES:
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

        strand = "+" if feature.location.strand >= 0 else "-"
        gene_id = (
            feature.qualifiers.get("gene", [None])[0]
            or feature.qualifiers.get("locus_tag", [None])[0]
            or feature.qualifiers.get("sequence", [None])[0]
            or f"{feature.type}_{i}"
        )
        genes.append({"id": gene_id, "start": start, "end": end, "strand": strand})
    return genes


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
# Circular figure
# =============================================================================

def _gene_arc_theta_width(gene, length):
    """(center_theta_degrees, width_degrees), wrap-aware."""
    start, end = gene["start"] % length, gene["end"]
    span = gene["end"] - gene["start"] + 1
    center = (gene["start"] + span / 2.0) % length
    return (center / length) * 360.0, (span / length) * 360.0


def _build_circular_figure(genes, motif_tracks, length, title):
    fig = go.Figure()

    gene_theta = [_gene_arc_theta_width(g, length)[0] for g in genes]
    gene_width = [_gene_arc_theta_width(g, length)[1] for g in genes]
    gene_hover = [f"{g['id']} ({g['strand']} strand)<br>{g['start']}-{g['end'] % length}" for g in genes]

    n_tracks = len(motif_tracks)
    outer_r, inner_r = 0.98, 0.55
    gene_r = 1.0

    if genes:
        fig.add_trace(go.Barpolar(
            r=[0.06] * len(genes),
            theta=gene_theta,
            width=gene_width,
            base=gene_r - 0.06,
            marker_color=GENE_COLOR,
            marker_line_width=0,
            name=GENE_ROW_LABEL,
            hovertext=gene_hover,
            hoverinfo="text",
            opacity=0.9,
        ))

    track_step = (outer_r - inner_r) / max(n_tracks, 1)
    for i, track in enumerate(motif_tracks):
        radius = outer_r - i * track_step
        thetas = [(p["position"] % length) / length * 360.0 for p in track["points"]]
        radii = [radius] * len(thetas)
        fig.add_trace(go.Scatterpolar(
            r=radii,
            theta=thetas,
            mode="markers",
            marker=dict(
                color=track.get("point_color", track["color"]),
                size=track.get("size", 9),
                symbol=track.get("symbol", "circle"),
                opacity=track.get("opacity", 1.0),
                line=dict(width=1, color="white"),
            ),
            name=track["label"],
            hovertext=[p["hover"] for p in track["points"]],
            hoverinfo="text",
        ))

    fig.update_layout(
        title=title,
        polar=dict(
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

def _build_linear_figure(genes, motif_tracks, length, title):
    fig = go.Figure()

    n_tracks = len(motif_tracks)
    gene_y = n_tracks + 1
    row_height = 0.8

    for g in genes:
        end = min(g["end"], length - 1)
        fig.add_shape(
            type="rect",
            x0=g["start"], x1=end,
            y0=gene_y - row_height / 2, y1=gene_y + row_height / 2,
            fillcolor=GENE_COLOR, line=dict(width=0), opacity=0.9,
        )
        fig.add_trace(go.Scatter(
            x=[(g["start"] + end) / 2], y=[gene_y],
            mode="markers", marker=dict(opacity=0),
            hovertext=[f"{g['id']} ({g['strand']} strand)<br>{g['start']}-{end}"],
            hoverinfo="text", showlegend=False,
        ))

    for i, track in enumerate(motif_tracks):
        y = n_tracks - i
        xs = [p["position"] for p in track["points"]]
        fig.add_trace(go.Scatter(
            x=xs, y=[y] * len(xs),
            mode="markers",
            marker=dict(
                color=track.get("point_color", track["color"]),
                size=track.get("size", 10),
                symbol=track.get("symbol", "circle"),
                opacity=track.get("opacity", 1.0),
                line=dict(width=1, color="white"),
            ),
            name=track["label"],
            hovertext=[p["hover"] for p in track["points"]],
            hoverinfo="text",
        ))

    fig.update_layout(
        title=title,
        xaxis=dict(title="Position (bp)", range=[0, length]),
        yaxis=dict(visible=False, range=[0, gene_y + 1]),
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
    distinct pattern, plotted against the unedited construct's genes.

    fig_after: the same footprint, but each marker styled by outcome —
    solid = still present (unresolved), faded/hollow = resolved (gone),
    plus an 'x' marker in a shared alert color for anything
    _final_new_motif_check found that wasn't there originally. Hover text
    on fig_after is the actual decision-matrix reasoning for that motif
    (resolved_motif_keys is the authoritative status — it's built during
    the real edit loop, since a motif resolved as a side effect of a
    different edit never gets a decision-matrix row of its own to infer
    status from).
    """
    genes = _extract_genes(output_record)
    patterns = {m.motif for m in motifs} | {nm["motif"] for nm in new_motifs}
    colors = _assign_pattern_colors(patterns)
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

    fig_before = builder(genes, before_tracks, sequence_length,
                          title or "Motifs before SyToGen")

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
            "symbol": NEW_MOTIF_SYMBOL_CIRCULAR if topology == "circular" else NEW_MOTIF_SYMBOL_LINEAR,
            "size": 12,
        })

    fig_after = builder(genes, after_tracks, sequence_length,
                         title or "Motifs after SyToGen")

    return fig_before, fig_after
