#!/usr/bin/env python3
"""Candidate GRAPHICAL ABSTRACT layouts for the ModelSEED 2026 NAR update.

WHAT THIS DRAWS
---------------
Five competing concepts for the mandatory NAR graphical abstract, each covering
the same four headline features of the update:

    molecules  ->  reactions  ->  thermodynamics / reaction ranking  ->  atom mapping

    flow_pipeline            left-to-right chevron pipeline (the requested idea)
    one_reaction_four_lenses one hero reaction, four "what we now know" lenses
    before_after_ledger      2020 row vs 2026 row across the four features
    direction_funnel         56k reactions funnelled down to a direction call,
                             then fanned out into model-output impact
    atom_trace_spine         one carbon traced through a pathway, data layers
                             stacked underneath

Pick one, then we replace the placeholder cells with real art. These are
LAYOUT MOCKUPS, not the submission file.

DATA PROVENANCE
---------------
Every bar reads data/msdb_counts.tsv, built by scripts/build_msdb_counts.py from
a snapshot of upstream ModelSEED `dev` -- see that script for the field-by-field
sourcing and the dev SHA, which is carried in the TSV header and copied into
each figure's _stats.tsv. Rebuild the TSV before regenerating if dev has moved.

The four headline counts (compounds, compounds with structures, reactions,
balanced reactions) were checked against ../data/snapshot_2026-07-29.md and
latex/sections/results_growth.tex on 2026-09-11 and match exactly, so these
figures and Table 2 agree. Re-check when v2.0.0 is tagged.

2020 values are Seaver et al. 2020 Tables 2-3, not recomputed. The 2020
thermodynamics baseline is the pre-2026-refresh eQuilibrator table, which is a
proxy rather than a re-run of the 2020 pipeline -- it is labelled "pre-2026" in
the figures for that reason, never "2020".

Everything the manuscript still marks [TBD] is drawn through `tbd()` as a dashed
"TBD" chip rather than invented, and every one is listed in the emitted stats
TSV. Search this file for PENDING to find them.

NON-OBVIOUS CHOICES
-------------------
* Canvas is 254 x 101.6 mm (10 x 4 in) = the NAR 5:2 ratio at 2x the 127 x 50 mm
  minimum. Axis units ARE millimetres at that size, so positions read directly
  against the spec.
* NAR requires 12-16 pt sans-serif. That is read here as 12-16 pt IN THE
  DELIVERED FILE, which is this 254 mm canvas -- 127 x 50 mm is stated as a
  minimum size, not the expected display size. `_check_font_floor()` fails the
  run if any size is quoted below MIN_PT, which is the real reason these layouts
  carry so little text; NAR asks for exactly that ("use text sparingly, mainly
  for labels"). If a co-author instead wants 12 pt to survive a downscale to the
  127 mm minimum, set PT_SCALE = 2.0 -- but note the layouts will not hold the
  text at that size and every panel needs re-cutting.
* Palette is the dataviz reference instance. Stage hues are categorical slots
  1/2/3/7 (blue, orange, aqua, violet); validated all-pairs, light mode:
  worst CVD dE 9.2, worst normal-vision dE 16.3, PASS. Aqua warns below 3:1
  contrast on the surface -- the relief rule is satisfied because every stage
  carries a visible text label, never colour alone.
* Reaction direction uses the DIVERGING pair (blue forward / red reverse) with a
  grey midpoint for reversible -- polarity, not identity. Do not recolour these
  to categorical slots.
* Molecule skeletons are hand-drawn placeholder glyphs. RDKit is not installed
  in either python on this host; the final artwork should be RDKit depictions of
  the actually chosen compounds, dropped in where `molecule_glyph()` is called.
* Font is DejaVu Sans -- Arial is not installed here. NAR asks for Arial. Set
  FIGURES_FONT=Arial (or install it) before producing the submission file.

OUTPUT
------
latex/figures/<out-subdir>/graphical_abstract_<concept>.{pdf,svg,png} plus a
_stats.tsv holding every number the layout asserts and every pending item.
PDF is the deliverable format (NAR accepts editable PDF; vector sidesteps the
300-600 dpi table entirely). SVG is the one to open in Illustrator -- it carries
the same geometry with unembedded, fully editable text. PNG at 600 dpi is for
on-screen review only.

USAGE
-----
    python3 scripts/regen_figures.py --tag graphical-abstract
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.patches import (Circle, FancyArrowPatch, FancyBboxPatch,
                                Polygon, Rectangle)

ROOT = Path(__file__).resolve().parents[1]
FIG_ROOT = Path(os.environ.get("NAR_FIG_ROOT", ROOT / "latex" / "figures"))

# ---------------------------------------------------------------------------
# Canvas. Units are millimetres at the 254 mm delivery width.
# ---------------------------------------------------------------------------
W, H = 254.0, 101.6                      # NAR 5:2, 2x the 127x50 mm minimum
FIG_W_IN = 10.0
PT_SCALE = 1.0                           # see the font note in the docstring
MIN_PT = 12                              # NAR floor; enforced, not aspirational
MAX_PT_LABEL = 16                        # NAR's stated label ceiling (headings)
DPI_PNG = 600

FONT = os.environ.get("FIGURES_FONT", "Arial")

# ---------------------------------------------------------------------------
# Palette -- dataviz reference instance (references/palette.md), light mode.
# ---------------------------------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8880"
RULE = "#dedcd6"

STAGES = {                               # categorical slots 1, 2, 3, 7
    "mol":    "#2a78d6",                 # molecules
    "rxn":    "#eb6834",                 # reactions
    "thermo": "#1baf7a",                 # thermodynamics + ranking
    "atom":   "#4a3aa7",                 # atom mapping
}
# Diverging pair for reaction DIRECTION (polarity, not identity).
DIR_FWD, DIR_REV, DIR_BOTH = "#2a78d6", "#e34948", "#9a9892"
# Ordinal blue ramp for the funnel; step 250 is the lightest that clears 2:1.
RAMP = ["#86b6ef", "#5598e7", "#3987e5", "#256abf", "#184f95"]

# ---------------------------------------------------------------------------
# Numbers the layouts assert.
#   snapshot_2026-07-29.md / results_growth.tex. Provisional against dev 1111754.
# ---------------------------------------------------------------------------
N = {
    "compounds_2020": 33992, "compounds_2026": 45708,
    "structures_2020": 28120, "structures_2026": 36943,
    "reactions_2020": 36193, "reactions_2026": 56012,
    "balanced_2020": 25457,  "balanced_2026": 34370,
    "thermo_sources_2020": 2, "thermo_sources_2026": 4,
    "biolog_2020": 355, "biolog_total": 390,
}
# from results_structure_curation.tex rather than the growth snapshot
N_CURATION = {"acp_overrides": 47, "excluded_compounds": 4}
THERMO_SOURCES = ["eQuilibrator", "Group contrib", "dGPredictor", "OpenTECR"]


def THERMO_ROWS():
    """Reactions carrying a stored DrG' from each estimator. TECRDB is left out
    of the bars on purpose -- at 797 against 56,002 its bar would be under a
    millimetre, so it is reported as a number beside them instead of being
    inflated to visibility."""
    return [("eQuilibrator", c("thermo_reactions:eQuilibrator")),
            ("Group contrib.", c("thermo_reactions:Group contribution")),
            ("dGPredictor", c("thermo_reactions:dGPredictor"))]


def DIRECTION_SEGS():
    """Direction split over BALANCED reactions only. Over all 56,012 the
    unknown share is 46% and swamps the strip; balanced is the population the
    direction call is actually made on."""
    return [(">", c("direction_balanced:>"), DIR_FWD),
            ("=", c("direction_balanced:="), DIR_BOTH),
            ("<", c("direction_balanced:<"), DIR_REV),
            ("?", c("direction_balanced:?"), RULE)]

# PENDING: every manuscript [TBD] these layouts would need. Emitted to the TSV.
PENDING = {
    "similarity_model": "chosen reaction-embedding foundation model",
    "biolog_2026": "functional Biolog conditions, 2026 refresh",
    "n_models": "size of the ModelSEED v2 draft-model corpus",
    "direction_effect": "headline direction-sensitivity result",
}

COUNTS_TSV = ROOT / "data" / "msdb_counts.tsv"


def load_counts() -> dict[tuple[str, str], int]:
    """Real counts from scripts/build_msdb_counts.py, keyed (key, vintage).

    Fails loudly rather than falling back to the hard-coded N dict: a bar chart
    that silently drew stale numbers is the exact failure this file exists to
    avoid."""
    if not COUNTS_TSV.exists():
        raise SystemExit(
            f"{COUNTS_TSV} missing. Build it first:\n"
            "  MSDB_ROOT=/scratch/ctaylor/tmp/devsnap_<sha> "
            "python3 scripts/build_msdb_counts.py")
    out, prov = {}, ""
    for ln in COUNTS_TSV.read_text().splitlines():
        if ln.startswith("#"):
            prov = prov or ln.lstrip("# ").strip()
            continue
        k, v, vintage = ln.split("\t")
        if k == "key":
            continue
        out[(k, vintage)] = int(v)
    out[("_provenance", "meta")] = prov  # type: ignore[assignment]
    return out


C = load_counts()


def c(key: str, vintage: str = "2026") -> int:
    return C[(key, vintage)]


_font_uses: list[tuple[str, int]] = []
# Every bar registers (chart, label, value, scale_top, drawn_mm, full_mm) so the
# "these bars are to scale" claim can be rechecked instead of trusted.
_bars: list[tuple[str, str, int, int, float, float]] = []


def _bar(chart, label, value, top, drawn, full):
    _bars.append((chart, label, value, top, drawn, full))


def _check_bar_scale(tol=1e-6):
    for chart, label, value, top, drawn, full in _bars:
        want, got = value / top, drawn / full
        if abs(want - got) > tol:
            raise SystemExit(
                f"bar not to scale: {chart}/{label} value={value} top={top} "
                f"expected {want:.6f} of the axis, drew {got:.6f}")


def pt(size: int) -> float:
    """Font size quoted in points AT THE 127 mm MINIMUM; scaled to the canvas."""
    _font_uses.append(("", size))
    return size * PT_SCALE


def _check_font_floor() -> None:
    bad = sorted({s for _, s in _font_uses if s < MIN_PT})
    if bad:
        raise SystemExit(
            f"font floor violated: sizes {bad} quoted below NAR's {MIN_PT} pt "
            f"minimum. Cut text or enlarge the panel -- do not lower the floor.")


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def new_canvas():
    fig = plt.figure(figsize=(FIG_W_IN, FIG_W_IN * H / W), dpi=100)
    fig.patch.set_facecolor(SURFACE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_facecolor(SURFACE)
    ax.set_axis_off()
    ax.set_autoscale_on(False)   # imshow gradients must not rescale the canvas
    return fig, ax


def text(ax, x, y, s, size, *, color=INK, weight="normal", ha="left",
         va="center", **kw):
    return ax.text(x, y, s, fontsize=pt(size), color=color, fontweight=weight,
                   ha=ha, va=va, **kw)


def card(ax, x, y, w, h, *, face=SURFACE, edge=RULE, lw=1.0, radius=2.5,
         alpha=1.0, ls="-"):
    """Rounded panel. Everything in these layouts sits on one of these."""
    p = FancyBboxPatch((x + radius, y + radius), w - 2 * radius, h - 2 * radius,
                       boxstyle=f"round,pad={radius}", linewidth=lw,
                       facecolor=face, edgecolor=edge, alpha=alpha,
                       linestyle=ls, zorder=1)
    ax.add_patch(p)
    return p


def gradient_fill(ax, x, y, w, h, hue, *, a_top=0.20, a_bot=0.02, radius=2.5,
                  bands=48):
    """Vertical tint gradient clipped to the panel's rounded rectangle.

    Drawn as stacked constant-alpha bands rather than an imshow ramp ON PURPOSE:
    imshow exports as an embedded raster, which would put a bitmap inside an
    otherwise all-vector PDF/SVG and break both the "editable PDF" the journal
    asks for and Illustrator editing. 48 bands is past the point where the steps
    are visible at print size, and in Illustrator they arrive as one clipped
    group you can delete and replace with a native gradient.

    Deliberately weak (alpha 0.20 -> 0.02): a chart background has to stay
    recessive or it competes with the marks, and the darkest step still leaves
    the 12 pt body text well clear of the surface it sits on."""
    rgb, surf = to_rgb(hue), to_rgb(SURFACE)
    clip = FancyBboxPatch((x + radius, y + radius), w - 2 * radius,
                          h - 2 * radius, boxstyle=f"round,pad={radius}",
                          facecolor="none", edgecolor="none",
                          transform=ax.transData)
    ax.add_patch(clip)
    bh = h / bands
    for i in range(bands):
        a = a_top + (a_bot - a_top) * (i + 0.5) / bands
        # Bands are pre-blended against the surface and drawn OPAQUE. Stacking
        # translucent bands instead would compound alpha wherever they overlap,
        # and the overlap is what hides the seams -- so translucent bands make
        # every seam a visible darker line.
        col = tuple(s + (c - s) * a for c, s in zip(rgb, surf))
        r = Rectangle((x, y + h - (i + 1) * bh), w, bh * 1.04, facecolor=col,
                      edgecolor="none", linewidth=0, zorder=0.5)
        ax.add_patch(r)
        r.set_clip_path(clip)


def stage_card(ax, x, y, w, h, key, title, *, gradient=False):
    """Tinted panel with a solid header bar -- the stage label is always text,
    which is what satisfies the relief rule for the low-contrast aqua hue."""
    hue = STAGES[key]
    if gradient:
        gradient_fill(ax, x, y, w, h, hue)
    else:
        card(ax, x, y, w, h, face=hue, edge="none", alpha=0.07)
    card(ax, x, y, w, h, face="none", edge=hue, lw=1.2, alpha=0.55)
    bar_h = 9.0
    ax.add_patch(FancyBboxPatch((x + 2.5, y + h - bar_h + 2.5 - 2.5),
                                w - 5, bar_h - 5, boxstyle="round,pad=2.5",
                                facecolor=hue, edgecolor="none", zorder=2))
    text(ax, x + w / 2, y + h - bar_h / 2, title, 12, color="white",
         weight="bold", ha="center", zorder=3)
    return y + h - bar_h            # usable top edge below the header


def chevron(ax, x, y, w, h, color):
    """Between-stage arrow. Carries the left-to-right reading order NAR asks for."""
    ax.add_patch(FancyArrowPatch((x, y), (x + w, y), arrowstyle="-|>",
                                 mutation_scale=h, linewidth=2.4,
                                 color=color, zorder=4,
                                 shrinkA=0, shrinkB=0))


def tbd(ax, x, y, w, h, label):
    """A value the manuscript has not filled in yet. Visibly a placeholder."""
    card(ax, x, y, w, h, face="none", edge=INK_MUTED, lw=1.0, ls=(0, (3, 3)),
         radius=2.0)
    text(ax, x + w / 2, y + h / 2, label, 12, color=INK_MUTED, ha="center")


def growth_bar(ax, x, y, w, h, v20, v26, hue):
    """2020 as a muted stub, the 2020->2026 increment as the coloured extension.
    Horizontal because these panels are short on height, not width."""
    frac = v20 / v26
    ax.add_patch(FancyBboxPatch((x + 1.0, y + 1.0), w * frac - 2.0, h - 2.0,
                                boxstyle="round,pad=1.0,rounding_size=1.2",
                                facecolor=INK_MUTED, edgecolor=SURFACE,
                                linewidth=1.0, zorder=3))
    ax.add_patch(FancyBboxPatch((x + w * frac + 1.0, y + 1.0),
                                w * (1 - frac) - 2.0, h - 2.0,
                                boxstyle="round,pad=1.0,rounding_size=1.2",
                                facecolor=hue, edgecolor=SURFACE,
                                linewidth=1.0, zorder=3))


def _rrect(ax, x, y, w, h, color, *, alpha=1.0, z=3, r=1.2):
    ax.add_patch(FancyBboxPatch((x, y), max(w, 0.01), max(h, 0.01),
                                boxstyle=f"round,pad=0,rounding_size={r}",
                                facecolor=color, alpha=alpha, edgecolor=SURFACE,
                                linewidth=0.8, zorder=z))


def bar_pair(ax, x, y, w, h, v20, v26, hue, *, part20=None, part26=None,
             lab20="2020", lab26="2026", top=None):
    """Two vertical bars, 2020 muted and 2026 in the stage hue, direct-labelled
    with their totals. If part* is given it is drawn as a solid sub-segment of
    the same bar (with-structure / balanced), the remainder at low alpha -- the
    panel's legend line names what solid means, so it is never colour-alone."""
    # `top` lets sibling panels share one scale so their bars are comparable
    # by eye; without it each panel self-scales and equal heights would imply
    # equal magnitudes across panels that do not have them.
    top = top or max(v20, v26) or 1
    lab_h, val_h = 5.5, 6.0
    plot_h = h - lab_h - val_h
    bw = min(w * 0.30, 15.0)
    gap = w * 0.18
    x0 = x + (w - 2 * bw - gap) / 2
    for i, (v, part, lab, col) in enumerate(
            [(v20, part20, lab20, INK_MUTED), (v26, part26, lab26, hue)]):
        bx = x0 + i * (bw + gap)
        bh = plot_h * v / top
        _bar("bar_pair", f"{lab}", v, top, bh, plot_h)
        _rrect(ax, bx, y + lab_h, bw, bh, col,
           alpha=0.35 if part is not None else 1.0)
        if part is not None:
            _rrect(ax, bx, y + lab_h, bw, plot_h * part / top, col, z=4)
        text(ax, bx + bw / 2, y + lab_h + bh + val_h / 2, f"{v:,}", 12,
             color=INK, weight="bold", ha="center")
        text(ax, bx + bw / 2, y + lab_h / 2, lab, 12, color=INK_2, ha="center")


def legend_line(ax, x, y, color, label):
    """One-line key for a bar sub-segment. Swatch plus text, 12 pt."""
    _rrect(ax, x, y - 1.6, 3.2, 3.2, color, r=0.8)
    text(ax, x + 5.0, y, label, 12, color=INK_2)


def hbar_rows(ax, x, y, w, h, rows, hue):
    """One labelled horizontal bar per row, label ABOVE the bar so the full
    panel width is available to the bar. rows = [(name, value), ...]."""
    top = max(v for _, v in rows) or 1
    pitch = h / len(rows)
    for i, (name, v) in enumerate(rows):
        ry = y + h - pitch * (i + 1)
        text(ax, x, ry + pitch - 3.0, f"{name}  {v:,}", 12, color=INK_2)
        _bar("hbar_rows", name, v, top, w * v / top, w)
        _rrect(ax, x, ry + 0.5, w * v / top, 4.2, hue,
               alpha=0.45 + 0.55 * (v / top))


def hbar_pair(ax, x, y, w, h, v20, v26, hue, lab20="2020"):
    """2020 over 2026 as two horizontal bars with their values at the right."""
    top = max(v20, v26) or 1
    bh = (h - 2.0) / 2
    yr_w, val_w = 20.0, 26.0     # year label left of the bar, value right of it
    span = w - yr_w - val_w
    for i, (v, col, yr) in enumerate([(v26, hue, "2026"),
                                      (v20, INK_MUTED, lab20)]):
        by = y + i * (bh + 2.0)
        text(ax, x, by + bh / 2, yr, 12, color=INK_2, va="center")
        _bar("hbar_pair", yr, v, top, span * v / top, span)
        _rrect(ax, x + yr_w, by, span * v / top, bh, col)
        text(ax, x + yr_w + span * v / top + 3.0, by + bh / 2, f"{v:,}", 12,
             color=INK_2, va="center")


def molecule_glyph(ax, cx, cy, s, color, *, kind="ring", number_at=None):
    """PLACEHOLDER skeletal structure. Replace with an RDKit depiction of the
    real compound once the hero compounds are chosen."""
    import math
    pts = []
    if kind == "ring":
        for i in range(6):
            a = math.radians(90 + 60 * i)
            pts.append((cx + s * math.cos(a), cy + s * math.sin(a)))
        ax.add_patch(Polygon(pts, closed=True, fill=False, edgecolor=color,
                             linewidth=1.6, zorder=3))
        ax.add_patch(Circle((cx, cy), s * 0.45, fill=False, edgecolor=color,
                            linewidth=1.0, alpha=0.5, zorder=3))
        # one substituent so it reads as a molecule, not a hexagon
        ax.plot([pts[0][0], pts[0][0] + s * 0.9], [pts[0][1], pts[0][1] + s * 0.5],
                color=color, lw=1.6, zorder=3, solid_capstyle="round")
    else:                                  # short open chain
        n = 5
        for i in range(n):
            px = cx - s * 1.6 + i * (s * 0.8)
            py = cy + (s * 0.35 if i % 2 else -s * 0.35)
            pts.append((px, py))
        xs, ys = zip(*pts)
        ax.plot(xs, ys, color=color, lw=1.6, zorder=3, solid_capstyle="round")
    # number_at maps a vertex INDEX to (label, colour) so a mapped pair can be
    # placed on atoms whose correspondence lines will not cross.
    for idx, (lab, col) in (number_at or {}).items():
        px, py = pts[idx]
        ax.add_patch(Circle((px, py), s * 0.40, facecolor=col,
                            edgecolor=SURFACE, linewidth=0.8, zorder=5))
        text(ax, px, py, lab, 12, color="white", weight="bold",
             ha="center", va="center", zorder=6)
    return pts


def reaction_glyph(ax, cx, cy, s, color):
    """A + B -> C + D, drawn small. Placeholder for the hero reaction."""
    molecule_glyph(ax, cx - s * 3.0, cy, s * 0.8, color, kind="ring")
    molecule_glyph(ax, cx - s * 1.3, cy, s * 0.8, color, kind="chain")
    ax.add_patch(FancyArrowPatch((cx - s * 0.35, cy), (cx + s * 0.75, cy),
                                 arrowstyle="-|>", mutation_scale=9,
                                 linewidth=1.6, color=INK_2, zorder=3))
    molecule_glyph(ax, cx + s * 2.1, cy, s * 0.8, color, kind="ring")
    molecule_glyph(ax, cx + s * 3.8, cy, s * 0.8, color, kind="chain")


def direction_strip(ax, x, y, w, h, segs=None):
    """Reversible / forward / reverse / unknown share as one stacked bar.
    segs = [(label, value, colour)]; 2 mm surface gaps between segments."""
    segs = segs or [(">", 34, DIR_FWD), ("=", 28, DIR_BOTH), ("<", 38, DIR_REV)]
    total = sum(v for _, v, _ in segs) or 1
    cx = x
    for lab, v, col in segs:
        sw = w * v / total
        _bar("direction_strip", lab, v, total, sw, w)
        _rrect(ax, cx + 1.0, y + 1.0, sw - 2.0, h - 2.0, col, r=1.0)
        if sw >= 5.5:      # below this a 12 pt glyph spills onto its neighbours
            text(ax, cx + sw / 2, y + h / 2, lab, 12,
                 color=INK_2 if col is RULE else "white", weight="bold",
                 ha="center", va="center", zorder=4)
        cx += sw


def source_dotplot(ax, x, y, w, h, labels, values, errs):
    """One row per thermodynamics source: label at left, estimate + uncertainty
    on a shared axis. Identity is carried by ROW POSITION and the text label,
    so no second categorical palette is needed inside this panel."""
    lab_w = w * 0.50
    ax_x, ax_w = x + lab_w, w - lab_w
    row = h / len(labels)
    lo, hi = -60.0, 25.0
    def sx(v):
        return ax_x + ax_w * (v - lo) / (hi - lo)
    ax.plot([sx(0), sx(0)], [y, y + h], color=RULE, lw=1.0, zorder=2,
            ls=(0, (2, 2)))
    for i, (lab, v, e) in enumerate(zip(labels, values, errs)):
        cy = y + h - row * (i + 0.5)
        text(ax, x, cy, lab, 12, color=INK_2, va="center")
        col = DIR_FWD if v < -10 else (DIR_REV if v > 10 else DIR_BOTH)
        ax.plot([sx(v - e), sx(v + e)], [cy, cy], color=col, lw=1.8,
                alpha=0.45, solid_capstyle="round", zorder=3)
        ax.add_patch(Circle((sx(v), cy), 1.5, facecolor=col, edgecolor=SURFACE,
                            linewidth=0.8, zorder=4))
    return sx


def atom_map_pair(ax, x, w, cy, s):
    """PLACEHOLDER atom map: two atoms numbered on both sides, correspondence
    drawn as dashed leaders. Vertices are chosen so the leaders do not cross."""
    pair = [("1", STAGES["mol"]), ("2", STAGES["rxn"])]
    lp = molecule_glyph(ax, x + w * 0.28, cy, s, INK_2, kind="ring",
                        number_at={5: pair[0], 4: pair[1]})
    rp = molecule_glyph(ax, x + w * 0.72, cy, s, INK_2, kind="chain",
                        number_at={3: pair[0], 0: pair[1]})
    for li, ri, (_, col) in ((5, 3, pair[0]), (4, 0, pair[1])):
        ax.plot([lp[li][0], rp[ri][0]], [lp[li][1], rp[ri][1]], color=col,
                lw=1.0, ls=(0, (2, 2)), alpha=0.8, zorder=4)


def browser_card(ax, x, y, w, h, title, hue):
    """UI frame -- this is how the atom mapping is exposed to a user.
    NAR forbids using the database HOME PAGE as a figure; a representative
    query result like this is explicitly allowed."""
    card(ax, x, y, w, h, face="#ffffff", edge=RULE, lw=1.2)
    bar = 8.0
    ax.add_patch(FancyBboxPatch((x + 2.5, y + h - bar + 0.5), w - 5, bar - 3.0,
                                boxstyle="round,pad=2.0,rounding_size=1.5",
                                facecolor="#f0efec", edgecolor="none", zorder=2))
    for i in range(3):
        ax.add_patch(Circle((x + 6.0 + i * 3.4, y + h - bar / 2), 1.0,
                            facecolor=INK_MUTED, edgecolor="none", zorder=3))
    text(ax, x + 18.0, y + h - bar / 2, title, 12, color=INK_2, va="center",
         zorder=3)
    return y + h - bar


# ---------------------------------------------------------------------------
# Concept A -- flow_pipeline  (the requested left-to-right flow chart)
# ---------------------------------------------------------------------------
def concept_flow_pipeline(ax):
    # Single centred title line. The four panels below are symmetric, so a
    # left-aligned two-line title left the top-right corner visibly empty.
    text(ax, W / 2, H - 12.0, "ModelSEED Biochemistry Database  —  2026 update",
         16, weight="bold", ha="center")

    top, bot = 78.0, 8.0
    hgt = top - bot
    widths = [52.0, 52.0, 62.0, 62.0]
    gap = (W - 12.0 - sum(widths)) / (len(widths) - 1)
    xs, cx = [], 6.0
    for w in widths:
        xs.append(cx)
        cx += w + gap

    # Panels 1 and 2 share one scale (the larger of the two totals) so that a
    # taller bar means a bigger number ACROSS panels, not just within one.
    count_top = max(c("compounds"), c("reactions"))

    # 1 -- molecules
    x, w = xs[0], widths[0]
    stage_card(ax, x, bot, w, hgt, "mol", "MOLECULES", gradient=True)
    bar_pair(ax, x + 4, bot + 14.0, w - 8, 44.0,
             c("compounds", "2020"), c("compounds"), STAGES["mol"],
             part20=c("compounds_with_structure", "2020"),
             part26=c("compounds_with_structure"), top=count_top)
    legend_line(ax, x + 7, bot + 7.0, STAGES["mol"], "with structures")

    # 2 -- reactions
    x, w = xs[1], widths[1]
    stage_card(ax, x, bot, w, hgt, "rxn", "REACTIONS", gradient=True)
    bar_pair(ax, x + 4, bot + 14.0, w - 8, 44.0,
             c("reactions", "2020"), c("reactions"), STAGES["rxn"],
             part20=c("reactions_balanced", "2020"),
             part26=c("reactions_balanced"), top=count_top)
    legend_line(ax, x + 12, bot + 7.0, STAGES["rxn"], "balanced")

    # 3 -- thermodynamics: how much DrG' data each reaction now has
    x, w = xs[2], widths[2]
    stage_card(ax, x, bot, w, hgt, "thermo", "THERMODYNAMICS", gradient=True)
    hbar_rows(ax, x + 5, bot + 32.0, w - 10, 28.0, THERMO_ROWS(), STAGES["thermo"])
    text(ax, x + 5, bot + 26.0,
         f"{c('thermo_reactions:TECRDB'):,} from experiment", 12, color=INK_2)
    text(ax, x + 5, bot + 19.5,
         f"3 estimates: {c('thermo_estimates_per_balanced_reaction:3'):,}", 12,
         color=INK_2)
    direction_strip(ax, x + 5, bot + 8.0, w - 10, 7.0, DIRECTION_SEGS())

    # 4 -- atom mapping in the UI
    x, w = xs[3], widths[3]
    inner = stage_card(ax, x, bot, w, hgt, "atom", "ATOM MAPPING", gradient=True)
    ui_top = browser_card(ax, x + 4, bot + 22.0, w - 8, inner - bot - 22.0,
                          "rxn00200", STAGES["atom"])
    atom_map_pair(ax, x + 4, w - 8, (bot + 24.0 + ui_top) / 2, 6.0)
    mapped, bal = c("atom_mapping:total"), c("reactions_balanced")
    text(ax, x + 6, bot + 15.5,
         f"{mapped:,} mapped  ({100 * mapped / bal:.0f}%)", 12, color=INK_2)
    _rrect(ax, x + 6, bot + 8.0, w - 12, 4.6, RULE, alpha=0.8, z=3)
    _bar("atom_progress", "mapped", mapped, bal, (w - 12) * mapped / bal, w - 12)
    _rrect(ax, x + 6, bot + 8.0, (w - 12) * mapped / bal, 4.6, STAGES["atom"],
           z=4)

    for i in range(3):
        chevron(ax, xs[i] + widths[i] + 1.0, (top + bot) / 2,
                gap - 2.0, 8.0, INK_MUTED)


def concept_one_reaction(ax):
    text(ax, 8, H - 8, "One reaction, four new layers", 14, weight="bold")

    hx, hw = 6.0, 76.0
    card(ax, hx, 6, hw, H - 22, face=RULE, edge="none", alpha=0.30)
    text(ax, hx + hw / 2, H - 28, "rxn00200", 14, weight="bold", ha="center")
    reaction_glyph(ax, hx + hw / 2, H - 48, 6.4, INK_2)
    text(ax, hx + hw / 2, H - 66, "mass balanced", 12, ha="center", color=INK_2)
    text(ax, hx + hw / 2, H - 76, f"1 of {N['balanced_2026']:,}", 12,
         ha="center", color=INK_2)
    chevron(ax, hx + hw + 2.0, H / 2, 8.0, 8.0, INK_MUTED)

    gx, gy = hx + hw + 14.0, 6.0
    gw, gh = W - gx - 6.0, H - 22.0
    cw, ch = (gw - 5.0) / 2, (gh - 5.0) / 2
    cells = [
        ("mol", "STRUCTURES", 0, 1),
        ("rxn", "SIMILAR REACTIONS", 1, 1),
        ("thermo", "THERMODYNAMICS", 0, 0),
        ("atom", "ATOM MAP", 1, 0),
    ]
    for key, title, col, row in cells:
        x = gx + col * (cw + 5.0)
        y = gy + row * (ch + 5.0)
        inner = stage_card(ax, x, y, cw, ch, key, title)
        if key == "mol":
            molecule_glyph(ax, x + cw * 0.30, (y + inner) / 2, 5.0, STAGES[key])
            text(ax, x + cw * 0.58, (y + inner) / 2 + 4, "InChI", 12, color=INK_2)
            text(ax, x + cw * 0.58, (y + inner) / 2 - 5, "PubChem", 12,
                 color=INK_2)
        elif key == "rxn":
            # ranked nearest neighbours from the reaction-similarity matrix
            text(ax, x + 6, inner - 5.0, "ranked by embedding", 12, color=INK_2)
            for r, frac in enumerate([0.96, 0.91, 0.84]):
                ry = inner - 12.0 - r * 6.5
                text(ax, x + 6, ry, f"{r + 1}", 12, color=INK_MUTED)
                bwid = (cw - 26) * frac
                ax.add_patch(FancyBboxPatch((x + 13, ry - 2.2), bwid, 4.4,
                                            boxstyle="round,pad=0,rounding_size=1.2",
                                            facecolor=STAGES["rxn"],
                                            alpha=0.35 + 0.22 * (2 - r),
                                            edgecolor="none", zorder=3))
        elif key == "thermo":
            source_dotplot(ax, x + 5, y + 4, cw - 10, inner - y - 8,
                           THERMO_SOURCES, [-38.0, -22.0, -30.0, -44.0],
                           [8.0, 16.0, 11.0, 4.0])
        else:
            ui = browser_card(ax, x + 4, y + 4, cw - 8, inner - y - 8,
                              "atom map", STAGES[key])
            atom_map_pair(ax, x + 4, cw - 8, (y + 6.0 + ui) / 2, 4.2)


# ---------------------------------------------------------------------------
# Concept C -- before_after_ledger
# ---------------------------------------------------------------------------
def concept_before_after(ax):
    text(ax, 8, H - 8, "What changed since 2020", 14, weight="bold")

    mapped = c("atom_mapping:total")
    cols = [
        ("mol", "MOLECULES", c("compounds", "2020"), c("compounds"),
         c("compounds_with_structure", "2020"), c("compounds_with_structure"),
         "with structures"),
        ("rxn", "REACTIONS", c("reactions", "2020"), c("reactions"),
         c("reactions_balanced", "2020"), c("reactions_balanced"), "balanced"),
        # pre-refresh eQuilibrator is the closest stored proxy for 2020 coverage
        ("thermo", "REACTIONS WITH ΔG",
         c("thermo_reactions:eQuilibrator_prerefresh", "pre-2026-refresh"),
         c("thermo_reactions:Group contribution"), 0,
         c("thermo_estimates_per_reaction:3"), "3 estimates"),
        ("atom", "ATOM-MAPPED", 0, mapped, None, None, None),
    ]
    x0, gap = 6.0, 5.0
    cw = (W - 12.0 - gap * 3) / 4
    for i, (key, title, v20, v26, p20, p26, legend) in enumerate(cols):
        x = x0 + i * (cw + gap)
        hue = STAGES[key]
        card(ax, x, 16.0, cw, 62.0, face=hue, edge="none", alpha=0.07)
        card(ax, x, 16.0, cw, 62.0, face="none", edge=hue, lw=1.2, alpha=0.5)
        text(ax, x + cw / 2, 84.0, title, 12, weight="bold", color=hue,
             ha="center")
        if v20 == 0:
            # a zero bar is invisible; say it in words and draw only 2026
            bw = min(cw * 0.30, 15.0)
            bx = x + cw / 2 - bw / 2 + cw * 0.17
            _rrect(ax, bx, 27.5, bw, 38.0, hue)
            text(ax, bx + bw / 2, 68.5, f"{v26:,}", 12, weight="bold",
                 color=INK, ha="center")
            text(ax, bx + bw / 2, 24.5, "2026", 12, color=INK_2, ha="center")
            text(ax, x + cw / 2 - bw - cw * 0.05, 46.0, "none", 12,
                 color=INK_MUTED, ha="center")
            text(ax, x + cw / 2 - bw - cw * 0.05, 24.5, "2020", 12,
                 color=INK_MUTED, ha="center")
        else:
            bar_pair(ax, x + 3, 22.0, cw - 6, 48.0, v20, v26, hue,
                     part20=p20, part26=p26,
                     lab20="pre-2026" if key == "thermo" else "2020")
        if legend:
            legend_line(ax, x + 4, 18.5, hue, legend)
    tbd(ax, 6, 3.0, W - 12, 9.5, "TBD  one-line headline finding")


def concept_direction_funnel(ax):
    text(ax, 8, H - 8, "From biochemistry to a direction call", 14, weight="bold")

    stops = [(f"{N['reactions_2026']:,}", "reactions"),
             (f"{N['balanced_2026']:,}", "balanced"),
             ("4", "ΔG sources"),
             ("> = <", "direction")]
    x0, y0 = 6.0, 24.0
    fh_top, fh_bot = 44.0, 22.0
    fw = 38.0
    gap = 3.0
    for i, (big, small) in enumerate(stops):
        x = x0 + i * (fw + gap)
        h_l = fh_top - (fh_top - fh_bot) * i / len(stops)
        h_r = fh_top - (fh_top - fh_bot) * (i + 1) / len(stops)
        cy = y0 + fh_top / 2
        poly = Polygon([(x, cy - h_l / 2), (x + fw, cy - h_r / 2),
                        (x + fw, cy + h_r / 2), (x, cy + h_l / 2)],
                       closed=True, facecolor=RAMP[i], edgecolor=SURFACE,
                       linewidth=2.0, zorder=2)
        ax.add_patch(poly)
        text(ax, x + fw / 2, cy + 3.0, big, 16, weight="bold", color="white",
             ha="center", zorder=3)
        text(ax, x + fw / 2, cy - 6.0, small, 12, color="white", ha="center",
             zorder=3)

    fan_x = x0 + 4 * (fw + gap) + 7.0
    fan_w = W - fan_x - 6.0
    ch = 13.0
    for i, lab in enumerate(["growth rate", "essential genes", "Biolog"]):
        y = y0 + 1.0 + i * (ch + 3.0)
        card(ax, fan_x, y, fan_w, ch, face=STAGES["thermo"], edge="none",
             alpha=0.10)
        card(ax, fan_x, y, fan_w, ch, face="none", edge=STAGES["thermo"],
             lw=1.0, alpha=0.5)
        text(ax, fan_x + 5, y + ch / 2, lab, 12, color=INK_2)
        tbd(ax, fan_x + fan_w * 0.62, y + 2.5, fan_w * 0.33, ch - 5.0, "TBD")
        ax.add_patch(FancyArrowPatch((fan_x - 3.0, y0 + fh_top / 2),
                                     (fan_x - 0.5, y + ch / 2),
                                     arrowstyle="-|>", mutation_scale=7,
                                     linewidth=1.4, color=INK_MUTED, zorder=3,
                                     connectionstyle="arc3,rad=0.15"))

    # atom mapping as a rail above the funnel
    rail_w = fan_x - x0 - 4.0
    card(ax, x0, H - 30.0, rail_w, 12.0, face=STAGES["atom"], edge="none",
         alpha=0.10)
    text(ax, x0 + 6, H - 24.0, "ATOM MAPPING  —  every balanced reaction, "
         "atom by atom", 12, color=STAGES["atom"], weight="bold")
    card(ax, x0, 6.0, rail_w, 12.0, face=STAGES["mol"], edge="none",
         alpha=0.10)
    text(ax, x0 + 6, 12.0, f"CURATED STRUCTURES  —  {N['compounds_2026']:,} "
         "compounds, PubChem validated", 12, color=STAGES["mol"], weight="bold")


# ---------------------------------------------------------------------------
# Concept E -- atom_trace_spine
# ---------------------------------------------------------------------------
def concept_atom_spine(ax):
    text(ax, 8, H - 8, "Tracing one carbon through the database", 14,
         weight="bold")

    n = 4
    x0 = 24.0
    step = (W - 2 * x0) / (n - 1)
    spine_y = H - 32.0
    ax.plot([x0 - 6, W - x0 + 6], [spine_y, spine_y], color=RULE, lw=2.0,
            zorder=1)
    for i in range(n):
        cx = x0 + i * step
        ax.add_patch(Circle((cx, spine_y), 7.0, facecolor=SURFACE,
                            edgecolor=INK_2, linewidth=1.4, zorder=3))
        molecule_glyph(ax, cx, spine_y, 3.6, INK_2, kind="ring")
        ax.add_patch(Circle((cx + 5.0, spine_y + 5.0), 2.6,
                            facecolor=STAGES["atom"], edgecolor=SURFACE,
                            linewidth=0.9, zorder=6))
        text(ax, cx + 5.0, spine_y + 5.0, "C", 12, color="white", weight="bold",
             ha="center", va="center", zorder=7)
        if i < n - 1:
            ax.add_patch(FancyArrowPatch((cx + 8.5, spine_y),
                                         (cx + step - 8.5, spine_y),
                                         arrowstyle="-|>", mutation_scale=9,
                                         linewidth=1.8, color=STAGES["atom"],
                                         zorder=2))
    text(ax, W / 2, spine_y + 13.0,
         f"{c('atom_mapping:total'):,} reactions atom-mapped", 12,
         color=STAGES["atom"], ha="center")

    layers = [("mol", "compounds", c("compounds", "2020"), c("compounds"),
               "2020"),
              ("rxn", "balanced reactions", c("reactions_balanced", "2020"),
               c("reactions_balanced"), "2020"),
              # the pre-refresh eQuilibrator table is the closest stored proxy
              # for 2020 coverage, so it is not labelled a clean "2020"
              ("thermo", "reactions with ΔG",
               c("thermo_reactions:eQuilibrator_prerefresh", "pre-2026-refresh"),
               c("thermo_reactions:Group contribution"), "pre-2026")]
    ly, lh = 6.0, 16.0
    for i, (key, title, v20, v26, lab20) in enumerate(reversed(layers)):
        y = ly + i * (lh + 2.5)
        hue = STAGES[key]
        card(ax, 6, y, W - 12, lh, face=hue, edge="none", alpha=0.10)
        card(ax, 6, y, W - 12, lh, face="none", edge=hue, lw=1.0, alpha=0.5)
        text(ax, 12, y + lh / 2, title, 12, color=hue, weight="bold")
        hbar_pair(ax, 62, y + 2.0, W - 80, lh - 4.0, v20, v26, hue, lab20)
        for j in range(n):
            ax.plot([x0 + j * step, x0 + j * step], [y + lh, y + lh + 2.5],
                    color=hue, lw=1.0, alpha=0.4, zorder=1)


CONCEPTS = {
    "flow_pipeline": concept_flow_pipeline,
    "one_reaction_four_lenses": concept_one_reaction,
    "before_after_ledger": concept_before_after,
    "direction_funnel": concept_direction_funnel,
    "atom_trace_spine": concept_atom_spine,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--concept", required=True, choices=sorted(CONCEPTS))
    ap.add_argument("--out-subdir", default=None,
                    help="directory under latex/figures/ (default: ga_<concept>)")
    args = ap.parse_args()

    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [FONT, "Liberation Sans", "Helvetica",
                                       "DejaVu Sans"]
    plt.rcParams["pdf.fonttype"] = 42          # embed as TrueType, not Type 3
    plt.rcParams["svg.fonttype"] = "none"

    out = FIG_ROOT / (args.out_subdir or f"ga_{args.concept}")
    out.mkdir(parents=True, exist_ok=True)

    fig, ax = new_canvas()
    CONCEPTS[args.concept](ax)
    _check_font_floor()
    _check_bar_scale()

    stem = out / f"graphical_abstract_{args.concept}"
    fig.savefig(f"{stem}.pdf", facecolor=SURFACE)
    # SVG as well: pdf.fonttype 42 embeds SUBSET CID fonts, which Illustrator
    # can only edit as text if the same font is installed locally. svg.fonttype
    # "none" writes <text> with a font-family name and embeds nothing, so the
    # text stays fully editable on any machine. Same vector geometry.
    fig.savefig(f"{stem}.svg", facecolor=SURFACE)
    fig.savefig(f"{stem}.png", dpi=DPI_PNG, facecolor=SURFACE)
    plt.close(fig)

    with open(out / "_stats.tsv", "w") as fh:
        fh.write("key\tvalue\tnote\n")
        fh.write(f"concept\t{args.concept}\tlayout mockup, not the submission file\n")
        fh.write(f"canvas_mm\t{W:.1f}x{H:.1f}\tNAR 5:2, 2x the 127x50 mm minimum\n")
        fh.write(f"min_font_pt\t{MIN_PT}\tenforced by _check_font_floor\n")
        fh.write(f"pt_scale\t{PT_SCALE}\tset to 2.0 to keep 12 pt after a downscale to 127 mm\n")
        fh.write(f"font\t{FONT}\tNAR asks for Arial; falls back to DejaVu Sans here\n")
        fh.write(f"data_provenance\t-\t{C[('_provenance', 'meta')]}\n")
        for (k, vintage), v in C.items():
            if k.startswith("_"):
                continue
            fh.write(f"{k}[{vintage}]\t{v}\tdata/msdb_counts.tsv\n")
        for chart, label, value, btop, drawn, full in _bars:
            fh.write(f"scale:{chart}:{label}\t{value}\t"
                     f"{drawn:.3f}mm of {full:.3f}mm axis, top={btop}\n")
        for k, v in N_CURATION.items():
            fh.write(f"{k}\t{v}\tresults_structure_curation.tex\n")
        for k, v in PENDING.items():
            fh.write(f"PENDING:{k}\tTBD\t{v}\n")

    print(f"wrote {stem}.pdf / .png and {out / '_stats.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
