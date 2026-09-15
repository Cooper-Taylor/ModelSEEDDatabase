#!/usr/bin/env python3
"""Candidate GRAPHICAL ABSTRACT layouts for the ModelSEED 2026 NAR update.

WHAT THIS DRAWS
---------------
Six competing concepts for the mandatory NAR graphical abstract, each covering
the same four headline features of the update:

    molecules  ->  reactions  ->  thermodynamics / reaction ranking  ->  atom mapping

    flow_pipeline            left-to-right chevron pipeline (the requested idea)
    one_reaction_four_lenses one hero reaction, four "what we now know" lenses
    before_after_ledger      2020 row vs 2026 row across the four features
    direction_funnel         56k reactions funnelled down to a direction call,
                             then fanned out into model-output impact
    atom_trace_spine         one carbon traced through a pathway, data layers
                             stacked underneath
    server_hub               the database as a hub: molecules/reactions/
                             structures in from the left, atom mapping out the
                             bottom, four DrG' sources out the right merging
                             into one prediction that is graded gold/silver/
                             bronze. Drawn from the hand sketch in
                             assets/sketch_server_hub_IMG_6387.jpeg.

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
import math
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.path import Path as MplPath
import matplotlib.image as mpimg
import numpy as np
from matplotlib.patches import (Circle, Ellipse, FancyArrowPatch,
                                FancyBboxPatch,
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
DPI_VECTOR = 1200      # resampling target for images inside the PDF/SVG

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
# Evidence tiers drawn as the metals they are named after. Validated light,
# all pairs: worst CVD dE 14.3, worst normal-vision dE 16.0 -- both clear.
# Two known, accepted deviations from the palette rules:
#   * silver fails the chroma floor (it IS grey) -- unavoidable given the
#     semantics, and safe here because every tier carries a text label, so
#     colour is never the only channel;
#   * gold sits at 2.05:1 against the surface, so the relief rule applies --
#     again satisfied by the direct labels.
GRADE_RAMP = ["#d4af37", "#8e9296", "#a9601f"]

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

# UI screen capture used in the atom-mapping panel. NAR forbids using the
# database HOME PAGE as a figure but explicitly allows "a representative screen
# dump" of a query result, which a reaction page is.
#
# CROP is (x0, x1, y0, y1) as FRACTIONS of the source image, applied before the
# capture is drawn. It exists for two reasons beyond framing:
#   * the site chrome carries the ModelSEED LOGO, and NAR bans trademarked
#     logos from the graphical abstract outright ("the text UniProt is fine,
#     but not the logo") -- so the nav bar must be cropped away, not just
#     shrunk;
#   * the same chrome shows the signed-in username, which should not ship.
# The default keeps the first two structures and the reaction arrow only.
SCREENSHOT = Path(os.environ.get(
    "NAR_UI_SHOT",
    ROOT / "latex" / "figures" / "ga_flow_pipeline" / "screenshot.png"))
# Crop measured against the committed 3420x2214 capture: keeps the WHOLE
# equation -- the (2) coefficient, GTP, the reversible arrow, PPi, H+ and
# GppppG, with their names and compound IDs. All four participants are needed
# for the atom mapping to make sense. Everything the crop removes is removed on
# purpose --
#   * the ModelSEED wordmark/logo in the site header (NAR ban),
#   * the signed-in username "ctaylor",
#   * the browser chrome, which shows personal bookmarks and a
#     staging.modelseed.org URL that must not appear in a published figure.
# Re-measure these if the capture is ever replaced; they are fractions of the
# source image, not pixels, but they are specific to this framing.
SCREENSHOT_CROP = (0.3140, 0.7300, 0.3784, 0.5521)

# server_hub reserves a slot for an atom-mapping illustration lifted from the
# web UI. Point this at a file and it is placed and reported; leave it absent
# and a dashed "paste here" slot is drawn instead, with the gap recorded in
# _stats.tsv. Never invent the artwork.
ATOM_MAP_IMAGE = Path(os.environ.get(
    "NAR_ATOM_MAP_IMAGE",
    ROOT / "assets" / "atom_mapping_capture.png"))
# server_hub also reserves a slot in the Computational box for a small bar
# chart of per-estimator coverage.
BAR_CHART_IMAGE = Path(os.environ.get(
    "NAR_BAR_CHART_IMAGE",
    ROOT / "assets" / "computational_bars.png"))
_atom_map_dpi: list[float] = []

# A CARTOON of the eQuilibrator reported-uncertainty histogram -- panel C, first
# chart, of Papers/NAR_Update_2026/figures/main_figures_draft.pdf (origin/dev).
# Traced from that PDF's vector path data: 30 bars over a 0-2 kcal/mol axis,
# downsampled here to 15 by taking the taller of each pair so the silhouette
# survives, then normalised to its own peak.
#
# THIS IS NOT DATA. It carries no axis, no counts and no scale, and it exists to
# say "most reactions have a small reported uncertainty, with a long tail". The
# real numbers are in that figure. _stats.tsv records it as a cartoon so nobody
# later mistakes it for a plotted series.
EQ_SIGMA_CARTOON = [0.74, 0.55, 1.00, 0.93, 0.58, 0.40, 0.28, 0.20,
                    0.14, 0.13, 0.09, 0.07, 0.03, 0.04, 0.02]
# The whole equation is 3.70:1 -- in a ~49 mm panel that is 11 mm tall and each
# structure lands at 8 mm. So the four participants are cut out INDIVIDUALLY and
# re-laid as a wrapped two-line equation, which more than doubles each structure
# at the same panel width. Each tile keeps its own compound name and ID.
# Fractions of the source image; all four share one y band.
MOLECULE_Y = (0.3799, 0.5498)
# Centres, not edges, with ONE shared half-width: hand-measured edges differed
# by a few pixels, which made the tiles different heights and left the row
# visibly unaligned.
# Half-width trimmed from 0.04125 to 0.0380 and each centre re-derived from the
# measured INK bounding box rather than the page's box edges. Same tile size on
# the page, less dead margin inside it, so the structures render larger.
# Checked: no crop reaches the reaction arrow or either "+" on the source page.
MOLECULE_HALF_W = 0.0380
MOLECULES = [("GTP", 0.36462), ("PPi", 0.48509),
             ("H+", 0.58436), ("GppppG", 0.68333)]


def molecule_crop(i):
    _, cx = MOLECULES[i]
    return (cx - MOLECULE_HALF_W, cx + MOLECULE_HALF_W, *MOLECULE_Y)


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


_panels: list[tuple[str, float, float, float, float]] = []


def _panel(name, x, y, w, h):
    _panels.append((name, x, y, w, h))


def _check_text_overlaps(fig, ax, pad=0.5, min_gap=1.5):
    """Fail the run if any two text labels overlap, or if one leaves the canvas.

    "Looks fine to me" does not survive a layout change three edits later, and
    at a 12 pt floor in 57 mm panels collisions are the default failure mode.
    Text over a BAR is fine and not checked -- only text against text."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    items = []
    for a in ax.texts:
        if id(a) in _nocheck:
            continue
        bb = a.get_window_extent(renderer=r).expanded(1.0, 1.0)
        items.append((a.get_text(), bb))
    bad = []
    for i in range(len(items)):
        ti, bi = items[i]
        for j in range(i + 1, len(items)):
            tj, bj = items[j]
            ov = (min(bi.x1, bj.x1) - max(bi.x0, bj.x0),
                  min(bi.y1, bj.y1) - max(bi.y0, bj.y0))
            # -min_gap, not 0: labels that merely fail to overlap still read
            # as one run-on string. Require visible daylight between them.
            if ov[0] > -min_gap and ov[1] > -min_gap:
                bad.append(f"{ti!r} too close to {tj!r} "
                           f"(gap {-max(ov):.1f}px, need {min_gap})")
    fw, fh = fig.canvas.get_width_height()
    for txt, bb in items:
        if bb.x0 < -pad or bb.y0 < -pad or bb.x1 > fw + pad or bb.y1 > fh + pad:
            bad.append(f"{txt!r} runs off the canvas")
    # A label may not leave the panel it sits in. Text-vs-text passes happily
    # while a caption spills into the neighbouring panel, so check this too.
    for name, px, py, pw, ph in _panels:
        p0 = ax.transData.transform((px, py))
        p1 = ax.transData.transform((px + pw, py + ph))
        for txt, bb in items:
            cx, cy = (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2
            if not (p0[0] <= cx <= p1[0] and p0[1] <= cy <= p1[1]):
                continue
            if bb.x0 < p0[0] - pad or bb.x1 > p1[0] + pad:
                over = max(p0[0] - bb.x0, bb.x1 - p1[0])
                bad.append(f"{txt!r} overflows panel {name} "
                           f"horizontally by {over:.1f}px")
            if bb.y0 < p0[1] - pad or bb.y1 > p1[1] + pad:
                over = max(p0[1] - bb.y0, bb.y1 - p1[1])
                bad.append(f"{txt!r} overflows panel {name} "
                           f"vertically by {over:.1f}px")
    if bad:
        raise SystemExit("text layout problems:\n  " + "\n  ".join(bad))



_containers: list[tuple] = []     # (x, y, w, h, name) in mm; see contain()


def contain(x, y, w, h, name, ellipse=False):
    """Register a box that text placed inside it must not escape.

    _check_text_overlaps() polices text against TEXT. Nothing policed text
    against the BOX it sits in, and that is the failure this figure kept
    producing: 'STRUCTURES' is 32.0 mm at 12 pt, not the 25 mm it looks like,
    and 'ΔG PREDICTIONS' is 40.6 mm. Both silently overhung their panels.
    Measure, do not estimate -- and then let the run fail if it regresses."""
    _containers.append((x, y, w, h, name, ellipse))


def _check_text_in_containers(fig, ax, pad=1.2):
    if not _containers:
        return
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bad = []
    for a in ax.texts:
        bb = a.get_window_extent(renderer=r)
        (x0, y0), (x1, y1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        for bx, by, bw, bh, name, is_ell in _containers:
            if not (bx <= cx <= bx + bw and by <= cy <= by + bh):
                continue                       # not this box's text
            if is_ell:
                ex, ey = bx + bw / 2, by + bh / 2
                ra, rb = bw / 2 - pad, bh / 2 - pad
                out = [(px, py) for px in (x0, x1) for py in (y0, y1)
                       if ((px - ex) / ra) ** 2 + ((py - ey) / rb) ** 2 > 1.0]
                if out:
                    bad.append(f"{a.get_text()!r} escapes {name}: "
                               f"{len(out)} of 4 corners outside the ellipse")
                continue
            if x0 < bx + pad or x1 > bx + bw - pad:
                over = max(bx + pad - x0, x1 - (bx + bw - pad))
                bad.append(f"{a.get_text()!r} overflows {name} horizontally "
                           f"by {over:.1f} mm")
            if y0 < by + pad or y1 > by + bh - pad:
                over = max(by + pad - y0, y1 - (by + bh - pad))
                bad.append(f"{a.get_text()!r} overflows {name} vertically "
                           f"by {over:.1f} mm")
    if bad:
        raise SystemExit("text escapes its box:\n  " + "\n  ".join(bad))


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


_capture_dpi: list[float] = []   # effective dpi of the placed UI capture
_nocheck: set[int] = set()


def text(ax, x, y, s, size, *, color=INK, weight="normal", ha="left",
         va="center", check=True, **kw):
    """check=False for labels drawn INSIDE a mark (atom numbers on their
    circles, > = < on their strip segments). Those are meant to sit close
    together and are separated by the mark, not by whitespace, so the
    collision audit must not police them."""
    a = ax.text(x, y, s, fontsize=pt(size), color=color, fontweight=weight,
                ha=ha, va=va, **kw)
    if not check:
        _nocheck.add(id(a))
    return a


def card(ax, x, y, w, h, *, face=SURFACE, edge=RULE, lw=1.0, radius=2.5,
         alpha=1.0, ls="-"):
    """Rounded panel. Everything in these layouts sits on one of these."""
    p = FancyBboxPatch((x + radius, y + radius), w - 2 * radius, h - 2 * radius,
                       boxstyle=f"round,pad={radius}", linewidth=lw,
                       facecolor=face, edgecolor=edge, alpha=alpha,
                       linestyle=ls, zorder=1)
    ax.add_patch(p)
    return p


def flask_glyph(ax, cx, cy, h, hue, *, lw=1.2, z=4):
    """An Erlenmeyer flask, h tall and 0.92h wide, centred on (cx, cy).

    Outline drawn as an OPEN polyline so the mouth stays open, with the liquid
    as a separate translucent fill: a closed outline would put a lid on it."""
    def P(pts):
        return [(cx + px * h, cy + py * h) for px, py in pts]
    ax.add_patch(Polygon(
        P([(-0.20, 0.50), (-0.20, 0.08), (-0.46, -0.38), (-0.40, -0.50),
           (0.40, -0.50), (0.46, -0.38), (0.20, 0.08), (0.20, 0.50)]),
        closed=False, facecolor="none", edgecolor=hue, linewidth=lw,
        joinstyle="round", zorder=z))
    # The liquid's top corners sit ON the sloping wall: the wall runs from
    # (0.20, 0.08) to (0.46, -0.38), so at y = -0.16 it has reached x = 0.336.
    ax.add_patch(Polygon(
        P([(-0.336, -0.16), (-0.46, -0.38), (-0.40, -0.50), (0.40, -0.50),
           (0.46, -0.38), (0.336, -0.16)]),
        closed=True, facecolor=hue, edgecolor="none", alpha=0.45, zorder=z))
    ax.plot(*zip(*P([(-0.29, 0.50), (0.29, 0.50)])), color=hue, lw=lw,
            solid_capstyle="round", zorder=z)


def chip_glyph(ax, cx, cy, h, hue, *, lw=1.2, z=4):
    """A processor die with legs, h wide and h tall, centred on (cx, cy)."""
    die, pin, core = 0.34, 0.16, 0.15
    for f in (-0.18, 0.0, 0.18):
        ax.plot([cx - (die + pin) * h, cx + (die + pin) * h],
                [cy + f * h, cy + f * h], color=hue, lw=lw,
                solid_capstyle="round", zorder=z)
        ax.plot([cx + f * h, cx + f * h],
                [cy - (die + pin) * h, cy + (die + pin) * h], color=hue,
                lw=lw, solid_capstyle="round", zorder=z)
    ax.add_patch(Rectangle((cx - die * h, cy - die * h), 2 * die * h,
                           2 * die * h, facecolor=_blend(hue, 0.30),
                           edgecolor=hue, linewidth=lw, zorder=z))
    ax.add_patch(Rectangle((cx - core * h, cy - core * h), 2 * core * h,
                           2 * core * h, facecolor="none", edgecolor=hue,
                           linewidth=lw * 0.8, zorder=z))


def confluence_arrow(ax, x0, rows, mid, xb, xj, xtip, color, *, shaft=2.8,
                     head=4.5, head_hw=6.0, alpha=1.0, z=2):
    """A straight middle shaft that the outer shafts curve into.

    Each outer arm runs horizontally to xb, then crosses to the middle line by
    xj on a smoothstep (3t^2 - 2t^3). That curve has HORIZONTAL tangents at
    both ends, so an arm leaves its own row and arrives on the trunk flush --
    no corner marks where the ribbon was bent, and no join to hide.

    Drawn as separate solid shapes in one OPAQUE colour rather than as a single
    traced outline. They union seamlessly because the colour is flat, which
    keeps the confluence a matter of geometry instead of of tracing a polygon
    around three ribbons that merge.
    """
    xh = xtip - head
    def box(x, y, w, h):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=color, edgecolor="none",
                               alpha=alpha, zorder=z))
    # Abutting shapes leave a hairline where their antialiased edges meet, so
    # every piece overlaps the next by 0.3 mm. The overlap is invisible: the
    # colour is flat and opaque, and the shaft is far narrower than the head.
    box(x0, mid - shaft / 2, xh - x0 + 0.3, shaft)    # the middle arm, straight
    t = np.linspace(0.0, 1.0, 120)
    for y in rows:
        if abs(y - mid) < 1e-9:
            continue
        box(x0, y - shaft / 2, xb - x0 + 0.3, shaft)
        c = np.stack([xb + t * (xj - xb),
                      y + (mid - y) * (3 * t ** 2 - 2 * t ** 3)], axis=1)
        d = np.gradient(c, axis=0)
        n = np.stack([-d[:, 1], d[:, 0]], axis=1)
        n /= np.linalg.norm(n, axis=1, keepdims=True)
        rail = np.concatenate([c + n * shaft / 2, (c - n * shaft / 2)[::-1]])
        ax.add_patch(Polygon(list(map(tuple, rail)), closed=True,
                             facecolor=color, edgecolor="none", alpha=alpha,
                             zorder=z))
    ax.add_patch(Polygon([(xh, mid + head_hw), (xtip, mid),
                          (xh, mid - head_hw)], closed=True, facecolor=color,
                         edgecolor="none", alpha=alpha, zorder=z))


def card_down_arrow(ax, x, y, w, h, hue, *, arrow=True, radius=2.5, lw=1.4,
                    shaft=4.0, shaft_l=1.8, head=9.0, head_l=3.4,
                    fill_alpha=None, z=1):
    """A rounded panel, optionally with a block arrow extruded from its BOTTOM
    edge, as ONE closed outline.

    One path rather than a card with an arrow parked underneath it: the bottom
    border runs along to the shaft, turns down and round the arrowhead, and
    comes back up to carry on, so the arrow is unmistakably part of the box and
    there is no join to hide. Filled opaque and pre-blended, for the same
    reason the cards are -- a translucent shape over another translucent one
    darkens the overlap.
    """
    # FILL_A and _blend are defined further down the module, so the default
    # is resolved at call time rather than at def time.
    fill_alpha = FILL_A if fill_alpha is None else fill_alpha
    r = radius
    cx = x + w / 2

    def corner(ox, oy, a0, a1, n=12):
        t = np.linspace(a0, a1, n)
        return [(ox + r * math.cos(v), oy + r * math.sin(v)) for v in t]

    p = [(x + r, y)]
    if arrow:
        p += [(cx - shaft / 2, y),
              (cx - shaft / 2, y - shaft_l),
              (cx - head / 2, y - shaft_l),
              (cx, y - shaft_l - head_l),
              (cx + head / 2, y - shaft_l),
              (cx + shaft / 2, y - shaft_l),
              (cx + shaft / 2, y)]
    p += [(x + w - r, y)]
    p += corner(x + w - r, y + r, -math.pi / 2, 0.0)
    p += [(x + w, y + h - r)]
    p += corner(x + w - r, y + h - r, 0.0, math.pi / 2)
    p += [(x + r, y + h)]
    p += corner(x + r, y + h - r, math.pi / 2, math.pi)
    p += [(x, y + r)]
    p += corner(x + r, y + r, math.pi, 1.5 * math.pi)
    ax.add_patch(Polygon(p, closed=True, facecolor=_blend(hue, fill_alpha),
                         edgecolor=hue, linewidth=lw, joinstyle="round",
                         zorder=z))


def gradient_fill(ax, x, y, w, h, hue, *, a_in=0.24, a_out=0.02, radius=2.5,
                  bands=72, cx_f=0.5, cy_f=0.45):
    """Radial tint gradient clipped to the panel's rounded rectangle: strongest
    at the centre, fading to almost nothing at the corners.

    Drawn as concentric OPAQUE pre-blended ellipses, largest first, rather than
    an imshow ramp. Two reasons: imshow exports as an embedded raster, which
    would put a bitmap in an otherwise all-vector PDF; and stacking translucent
    shapes compounds alpha at every boundary, turning each one into a visible
    ring. Opaque overdraw has neither problem.

    Deliberately weak (0.24 -> 0.02): a chart background has to stay recessive
    or it competes with the marks."""
    rgb, surf = to_rgb(hue), to_rgb(SURFACE)
    clip = FancyBboxPatch((x + radius, y + radius), w - 2 * radius,
                          h - 2 * radius, boxstyle=f"round,pad={radius}",
                          facecolor="none", edgecolor="none",
                          transform=ax.transData)
    ax.add_patch(clip)
    ex, ey = x + w * cx_f, y + h * cy_f
    # outermost ellipse has to cover the corners, hence the sqrt(2) overshoot
    rx, ry = w * 1.45, h * 1.45
    for i in range(bands):
        f = i / (bands - 1)              # 0 at the rim, 1 at the centre
        a = a_out + (a_in - a_out) * f
        col = tuple(s + (c - s) * a for c, s in zip(rgb, surf))
        e = Ellipse((ex, ey), rx * (1 - f) + 0.001, ry * (1 - f) + 0.001,
                    facecolor=col, edgecolor="none", linewidth=0, zorder=0.5)
        ax.add_patch(e)
        e.set_clip_path(clip)


def stage_card(ax, x, y, w, h, key, title, *, gradient=False, radius=2.5):
    """Tinted panel with a header bar FLUSH to the card's top and sides.

    The header is a plain rectangle clipped to the card's rounded-rect path, so
    its top corners pick up the card radius and its bottom edge stays square.
    Insetting it instead (the obvious FancyBboxPatch approach) leaves a strip of
    panel background above and beside the header, which reads as a misprint.

    The stage label is always text, which is what satisfies the relief rule for
    the low-contrast aqua hue."""
    _panel(title, x, y, w, h)
    hue = STAGES[key]
    clip = FancyBboxPatch((x + radius, y + radius), w - 2 * radius,
                          h - 2 * radius, boxstyle=f"round,pad={radius}",
                          facecolor="none", edgecolor="none",
                          transform=ax.transData)
    ax.add_patch(clip)
    if gradient:
        gradient_fill(ax, x, y, w, h, hue, radius=radius)
    else:
        card(ax, x, y, w, h, face=hue, edge="none", alpha=0.07)
    card(ax, x, y, w, h, face="none", edge=hue, lw=1.2, alpha=0.55)
    bar_h = 9.0
    hdr = Rectangle((x, y + h - bar_h), w, bar_h, facecolor=hue,
                    edgecolor="none", linewidth=0, zorder=2)
    ax.add_patch(hdr)
    hdr.set_clip_path(clip)
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


def _rrect(ax, x, y, w, h, color, *, alpha=1.0, z=3, r=1.2, clip=None,
           edge=SURFACE, lw=0.8):
    p = FancyBboxPatch((x, y), max(w, 0.01), max(h, 0.01),
                       boxstyle=f"round,pad=0,rounding_size={r}",
                       facecolor=color, alpha=alpha, edgecolor=edge,
                       linewidth=lw, zorder=z)
    ax.add_patch(p)
    if clip is not None:
        p.set_clip_path(clip)
    return p


def tick_label(v):
    return "0" if v == 0 else f"{v // 1000}k"


def chart_frame(ax, x, y, w, h, *, top, ticks, horizontal=False,
                tick_labels=True):
    """Plot area for a bar chart: its own light surface inset from the panel
    gradient, gridlines at `ticks`, tick labels outside the box, and a stronger
    baseline at zero. The gridlines are what let a reader read a value off the
    chart instead of relying on the printed number."""
    _rrect(ax, x, y, w, h, "#ffffff", z=0.8, r=1.5)
    for v in ticks:
        f = v / top
        if horizontal:
            gx = x + w * f
            ax.plot([gx, gx], [y, y + h], color=RULE, lw=0.6, zorder=1.0)
            if tick_labels:
                text(ax, gx, y - 3.6, tick_label(v), 12, color=INK_MUTED,
                     ha="center", va="center")
        else:
            gy = y + h * f
            ax.plot([x, x + w], [gy, gy], color=RULE, lw=0.6, zorder=1.0)
            if tick_labels:
                text(ax, x - 1.8, gy, tick_label(v), 12, color=INK_MUTED,
                     ha="right", va="center")
    # zorder 5 puts the axis ABOVE the bars (z 3-4) so it reads as one
    # unbroken line through the chart instead of being chopped up by each bar
    # and its surface-coloured edge stroke.
    if horizontal:
        ax.plot([x, x], [y, y + h], color=INK_MUTED, lw=1.2, zorder=5)
    else:
        ax.plot([x, x + w], [y, y], color=INK_MUTED, lw=1.2, zorder=5)
    # Invisible clip for the bars: they are drawn overshooting the baseline by
    # their corner radius and clipped back to it, which squares off the end that
    # sits on the axis while leaving the data end rounded.
    clip = Rectangle((x, y), w, h, facecolor="none", edgecolor="none",
                     transform=ax.transData)
    ax.add_patch(clip)
    return clip


def bar_chart(ax, x, y, w, h, values, labels, hue, *, top, ticks, parts=None,
              val_y=None, xlab_y=None, tick_labels=True):
    """Vertical bars inside a chart_frame. values[0] is the reference series and
    is drawn muted; the rest take the panel hue. `parts` draws a solid
    sub-segment of the same bar (with-structures / balanced / clean)."""
    clip = chart_frame(ax, x, y, w, h, top=top, ticks=ticks,
                       tick_labels=tick_labels)
    n = len(values)
    # One equal slice per bar, bar centred in its slice. This maximises the
    # centre-to-centre distance, which is what the value labels above the bars
    # actually need -- packing the bars tighter collides the labels.
    slice_w = w / n
    bw = min(slice_w * 0.52, 14.0)
    for i, v in enumerate(values):
        bx = x + slice_w * (i + 0.5) - bw / 2
        bh = h * v / top
        col = INK_MUTED if i == 0 else hue
        part = (parts or [None] * n)[i]
        _bar("bar_chart", labels[i], v, top, bh, h)
        R = 1.2
        _rrect(ax, bx, y - R, bw, bh + R, col,
               alpha=0.35 if part is not None else 1.0, z=3, r=R, clip=clip)
        if part is not None:
            _rrect(ax, bx, y - R, bw, h * part / top + R, col, z=4, r=R,
                   clip=clip)
        if val_y is not None:
            text(ax, bx + bw / 2, val_y, f"{v:,}", 12, color=INK,
                 weight="bold", ha="center")
        if xlab_y is not None:
            text(ax, bx + bw / 2, xlab_y, labels[i], 12, color=INK_2,
                 ha="center")


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


def legend_line(ax, x, y, entries):
    """Key for bar sub-segments. entries = [(colour, alpha, label), ...], laid
    out left to right; returns the x it finished at."""
    for cols, alpha, label in entries:
        cols = cols if isinstance(cols, (list, tuple)) else [cols]
        # Several swatches share one label where the same meaning is drawn in
        # more than one colour -- the solid part of BOTH the 2020 and the 2026
        # bar means "with structures", so the key has to show both.
        for col in cols:
            _rrect(ax, x, y - 1.6, 3.2, 3.2, col, alpha=alpha, r=0.8)
            x += 4.2
        text(ax, x + 1.0, y, label, 12, color=INK_2)
        x += 1.0 + 2.35 * len(label) + 6.0
    return x


def hbar_rows(ax, x, y, w, h, rows, hue, *, top=None, ticks=(),
              tick_labels=True, colors=None):
    """One labelled horizontal bar per row inside a chart_frame, label ABOVE the
    bar so the full box width is available to the bar itself.

    `colors` gives an explicit per-row colour, used by the ranking panel where
    the rows are ORDINAL tiers and take a one-hue ramp. Without it, rows share
    the panel hue and vary only in alpha with their own value."""
    top = top or max(v for _, v in rows) or 1
    clip = chart_frame(ax, x, y, w, h, top=top, ticks=ticks, horizontal=True,
                       tick_labels=tick_labels)
    pitch = h / len(rows)
    for i, (name, v) in enumerate(rows):
        ry = y + h - pitch * (i + 1)
        text(ax, x + 1.5, ry + pitch - 3.2, f"{name} {v:,}", 12, color=INK_2)
        _bar("hbar_rows", name, v, top, w * v / top, w)
        col = hue if colors is None else colors[i]
        bar_h = min(5.0, pitch * 0.38)
        R = 1.2
        _rrect(ax, x - R, ry + 0.8, w * v / top + R, bar_h, col,
               alpha=1.0 if colors is not None else 0.45 + 0.55 * (v / top),
               z=3, r=R, clip=clip)


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
             ha="center", va="center", zorder=6, check=False)
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
                 ha="center", va="center", zorder=4, check=False)
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


def capture_aspect(crop=None, default=3.70):
    """Width/height of a cropped region, for laying out around it."""
    if not SCREENSHOT.exists():
        return default
    img = mpimg.imread(SCREENSHOT)
    ih, iw = img.shape[0], img.shape[1]
    fx0, fx1, fy0, fy1 = crop or SCREENSHOT_CROP
    return ((fx1 - fx0) * iw) / ((fy1 - fy0) * ih)


def ui_capture(ax, x, y, w, hue, *, radius=1.5, crop=None,
                border=True):
    """Draw the cropped UI screen capture, scaled to width `w`, top-aligned at
    `y`. Returns the height it consumed.

    If the file is absent it draws a labelled placeholder instead of failing --
    the rest of the figure is still worth looking at, and the missing asset is
    reported in _stats.tsv rather than silently omitted.

    NOTE: a screen capture is a raster. This is the one element that puts a
    bitmap into an otherwise all-vector PDF; see the module docstring."""
    if not SCREENSHOT.exists():
        h = w / 1.87
        card(ax, x, y - h, w, h, face="none", edge=INK_MUTED, lw=1.0,
             ls=(0, (3, 3)), radius=2.0)
        text(ax, x + w / 2, y - h / 2 + 2.2, "UI capture", 12, color=INK_MUTED,
             ha="center")
        text(ax, x + w / 2, y - h / 2 - 3.0, "not found", 12, color=INK_MUTED,
             ha="center")
        return h
    img = mpimg.imread(SCREENSHOT)
    ih, iw = img.shape[0], img.shape[1]
    fx0, fx1, fy0, fy1 = crop or SCREENSHOT_CROP
    img = img[int(fy0 * ih):int(fy1 * ih), int(fx0 * iw):int(fx1 * iw)]
    h = w * img.shape[0] / img.shape[1]
    # Effective resolution of the SOURCE pixels at the placed size. The export
    # dpi can upsample past this, but this is the real detail available, and
    # it is what has to clear NAR's 300 dpi floor for colour half-tones.
    _capture_dpi.append(img.shape[1] / w * 25.4)
    im = ax.imshow(img, extent=(x, x + w, y - h, y), aspect="auto",
                   interpolation="none", zorder=3)
    clip = FancyBboxPatch((x + radius, y - h + radius), w - 2 * radius,
                          h - 2 * radius, boxstyle=f"round,pad={radius}",
                          facecolor="none", edgecolor="none",
                          transform=ax.transData)
    ax.add_patch(clip)
    im.set_clip_path(clip)
    if border:
        card(ax, x, y - h, w, h, face="none", edge=hue, lw=0.8, radius=radius,
             alpha=0.5)
    return h


def equilibrium(ax, cx, cy, w=4.4, color=INK_2, *, lw=1.0, ms=5, sep=0.7):
    """Reversible-reaction symbol as two offset half-arrows.

    Defaults are the small inline size used inside the input chips; pass lw/ms/
    sep up for a standalone symbol in an equation."""
    for dy, x0, x1 in ((sep, -1, 1), (-sep, 1, -1)):
        ax.add_patch(FancyArrowPatch((cx + x0 * w / 2, cy + dy),
                                     (cx + x1 * w / 2, cy + dy),
                                     arrowstyle="-|>", mutation_scale=ms,
                                     linewidth=lw, color=color, zorder=4,
                                     shrinkA=0, shrinkB=0))


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
    text(ax, W / 2, H - 9.0, "ModelSEED Biochemistry Database  —  2026 update",
         16, weight="bold", ha="center")

    top, bot = 84.0, 5.0
    hgt = top - bot
    # Panels are sized to their content, not equally: the thermodynamics source
    # labels are the longest strings in the figure and set that panel's width.
    # Widths are set from MEASURED label widths at 12 pt, not guessed:
    # "Group contrib  56,002" is 47.4 mm and sets panel 3; a pair of 16.8 mm
    # value labels sets panels 1-2; "by ΔG evidence" sets panel 4.
    widths = [45.0, 43.0, 56.0, 37.0, 49.0]
    gap = (W - 12.0 - sum(widths)) / (len(widths) - 1)
    xs, cx = [], 6.0
    for w in widths:
        xs.append(cx)
        cx += w + gap

    BOX_Y, BOX_H = 26.0, 40.0
    # LEG_Y is the shared key row: the blue/orange swatch keys and the green
    # panels' caption lines all sit on it, so the bottom of every panel reads
    # as one line across the figure.
    VAL_Y, XLAB_Y, LEG_Y = 70.5, 20.0, 11.0
    COUNT_TOP = 60000
    COUNT_TICKS = (0, 20000, 40000, 60000)

    # 1 -- molecules  (carries the y axis for itself and panel 2)
    x, w = xs[0], widths[0]
    stage_card(ax, x, bot, w, hgt, "mol", "MOLECULES", gradient=True)
    bar_chart(ax, x + 4, BOX_Y, w - 8, BOX_H,
              [c("compounds", "2020"), c("compounds")], ["2020", "2026"],
              STAGES["mol"], top=COUNT_TOP, ticks=COUNT_TICKS,
              parts=[c("compounds_with_structure", "2020"),
                     c("compounds_with_structure")],
              val_y=VAL_Y, xlab_y=XLAB_Y, tick_labels=False)
    legend_line(ax, x + 2, LEG_Y,
                [([INK_MUTED, STAGES["mol"]], 1.0, "with structures")])

    # 2 -- reactions, same 0-60k scale as panel 1
    x, w = xs[1], widths[1]
    stage_card(ax, x, bot, w, hgt, "rxn", "REACTIONS", gradient=True)
    bar_chart(ax, x + 4, BOX_Y, w - 8, BOX_H,
              [c("reactions", "2020"), c("reactions")], ["2020", "2026"],
              STAGES["rxn"], top=COUNT_TOP, ticks=COUNT_TICKS,
              parts=[c("reactions_balanced", "2020"),
                     c("reactions_balanced")],
              val_y=VAL_Y, xlab_y=XLAB_Y, tick_labels=False)
    legend_line(ax, x + 5, LEG_Y,
                [([INK_MUTED, STAGES["rxn"]], 1.0, "balanced")])

    # 3 -- thermodynamics: how much DrG' data each reaction has
    x, w = xs[2], widths[2]
    stage_card(ax, x, bot, w, hgt, "thermo", "THERMODYNAMICS", gradient=True)
    hbar_rows(ax, x + 4, BOX_Y, w - 8, BOX_H, THERMO_ROWS(), STAGES["thermo"],
              top=COUNT_TOP, ticks=COUNT_TICKS)
    text(ax, x + 4, LEG_Y, "reactions with a ΔG", 12, color=INK_2)

    # 4 -- ranking: the same reactions ordered by how good that DrG' is.
    # Same hue as panel 3 on purpose -- this ranks the thermodynamics, it is
    # not a fifth independent topic.
    x, w = xs[3], widths[3]
    stage_card(ax, x, bot, w, hgt, "thermo", "RANKING", gradient=True)
    tiers = [("gold", c("thermo_evidence_grade:gold")),
             ("silver", c("thermo_evidence_grade:silver")),
             ("bronze", c("thermo_evidence_grade:bronze"))]
    hbar_rows(ax, x + 4, BOX_Y, w - 8, BOX_H, tiers, STAGES["thermo"],
              top=max(v for _, v in tiers), ticks=(0, 10000),
              colors=GRADE_RAMP)
    text(ax, x + 4, LEG_Y, "ΔG evidence", 12, color=INK_2)

    # 5 -- atom mapping: one headline quantity and a split, stated not plotted
    x, w = xs[4], widths[4]
    inner5 = stage_card(ax, x, bot, w, hgt, "atom", "ATOM MAPPING",
                        gradient=True)
    # The four structures, cut out individually and re-laid as a wrapped
    # equation:   2 [GTP] <=> [PPi]
    #               + [H+]  +  [GppppG]
    # Operators are drawn here rather than taken from the capture, so they sit
    # on the layout grid instead of wherever the page happened to put them.
    inner_w = w - 6
    OP2, OPEQ, OPPLUS = 2.8, 5.0, 3.6        # widths reserved for 2, <=>, +
    tile_w = (inner_w - max(OP2 + OPEQ, 2 * OPPLUS)) / 2
    tile_h = tile_w / capture_aspect(molecule_crop(0))
    row_gap = 2.0
    block_top = (bot + inner5) / 2 + (2 * tile_h + row_gap) / 2

    def tile(i, tx, ty):
        ui_capture(ax, tx, ty, tile_w, STAGES["atom"], border=False,
                   crop=molecule_crop(i))
        return tx + tile_w

    # structures sit in the upper ~2/3 of each tile; operators align to them
    def op_y(ty):
        return ty - tile_h * 0.34

    tx = x + 3
    text(ax, tx, op_y(block_top), "2", 12, color=INK_2, va="center")
    tx += OP2
    tx = tile(0, tx, block_top)
    equilibrium(ax, tx + OPEQ / 2, op_y(block_top), w=3.8)
    tx += OPEQ
    tile(1, tx, block_top)

    row2 = block_top - tile_h - row_gap
    tx = x + 3
    text(ax, tx, op_y(row2), "+", 12, color=INK_2, va="center")
    tx += OPPLUS
    tx = tile(2, tx, row2)
    text(ax, tx + OPPLUS / 2, op_y(row2), "+", 12, color=INK_2, ha="center",
         va="center")
    tx += OPPLUS
    tile(3, tx, row2)


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
            legend_line(ax, x + 4, 18.5, [(hue, 1.0, legend)])
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
             ha="center", va="center", zorder=7, check=False)
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



# --------------------------------------------------------------- server_hub
SERVER_BODY = "#2f2f33"      # near-black chassis; the sketch asks for gray/black
SERVER_UNIT = "#54545b"      # the individual rack units
SERVER_EDGE = "#1d1d21"
LED_ON = "#17c07f"           # every light is lit; a saturated step of STAGES["thermo"]
# Structures are a property OF molecules, so they take the dark end of the same
# blue ramp rather than a new hue -- related, but not mistakable for the
# molecules chip now that the counts are gone.
STRUCT_HUE = RAMP[4]
FILL_A = 0.16                # chip fills; was 0.10 and read as washed out
ARROW_A = 1.0                # arrows at full strength


def server_rack(ax, x, y, w, h, *, title="ModelSEED", units=5):
    """A rack of servers, near-black, with one green status LED.

    Drawn rather than clip-arted so it stays vector and recolourable. The units
    are horizontal because that is what makes a box read as a rack; the sketch's
    vertical divisions read as a filing cabinet at small sizes."""
    # card() hardcodes zorder=1; the chassis has to sit above the inbound
    # arrows, so raise the returned patch rather than widening the shared helper.
    card(ax, x, y, w, h, face=SERVER_BODY, edge=SERVER_EDGE, lw=1.2,
         radius=2.5).set_zorder(3)
    hdr_h = 10.0
    ax.add_patch(FancyBboxPatch(
        (x + 1.4, y + h - hdr_h + 1.4), w - 2.8, hdr_h - 2.8,
        boxstyle="round,pad=1.4", facecolor=SERVER_EDGE, edgecolor="none",
        zorder=4))
    text(ax, x + w / 2, y + h - hdr_h / 2, title, 13, color="white",
         weight="bold", ha="center", va="center", zorder=5)

    body_top, body_bot = y + h - hdr_h - 2.0, y + 4.0
    pitch = (body_top - body_bot) / units
    pad, gap = 0.9, 1.6
    for i in range(units):
        # The DRAWN unit rect, derived once. FancyBboxPatch's pad expands the
        # box it is given, so the rect is inset by `pad` on every side of the
        # call below -- everything inside the unit is centred on THIS, which is
        # what the LEDs were previously 0.8 mm out from.
        ux0, ux1 = x + 3.0, x + w - 3.0
        uy0 = body_bot + i * pitch
        uy1 = uy0 + pitch - gap
        ucy = (uy0 + uy1) / 2.0
        ax.add_patch(FancyBboxPatch(
            (ux0 + pad, uy0 + pad), (ux1 - ux0) - 2 * pad, (uy1 - uy0) - 2 * pad,
            boxstyle=f"round,pad={pad}", facecolor=SERVER_UNIT, edgecolor="none",
            zorder=4))
        slot_h = (uy1 - uy0) * 0.40
        for k in range(3):
            ax.add_patch(Rectangle((ux1 - 10.0 + k * 3.0, ucy - slot_h / 2),
                                   1.5, slot_h, facecolor=SERVER_EDGE,
                                   edgecolor="none", alpha=0.8, zorder=5))
        ax.add_patch(Circle((ux0 + 3.6, ucy), 1.25, facecolor=LED_ON,
                            edgecolor="none", zorder=5))
    # No subtitle: the figure title already carries "Biochemistry Database",
    # and anything at the chassis foot lands on the lowest source row.


def block_arrow(ax, x0, y0, x1, y1, color, *, shaft=5.0, head=9.0, head_w=None,
                alpha=1.0, z=2):
    """A fat block arrow. Horizontal or vertical only -- the layouts never need
    a diagonal one, and a general implementation would be harder to reason
    about than two cases.

    `head` is the head LENGTH along the axis; `head_w` its WIDTH across it.
    They default to the same value, which is fine while shaft << head, but a
    shaft of 4.6 against a head of 5.0 produced a rectangle with a faint point
    rather than an arrow. Set head_w explicitly for short, fat arrows."""
    head_w = head if head_w is None else head_w
    if abs(y1 - y0) < 1e-9:                      # horizontal
        d = 1.0 if x1 > x0 else -1.0
        bx = x1 - d * head
        ax.add_patch(Polygon([(x0, y0 - shaft / 2), (bx, y0 - shaft / 2),
                              (bx, y0 - head_w / 2), (x1, y0),
                              (bx, y0 + head_w / 2), (bx, y0 + shaft / 2),
                              (x0, y0 + shaft / 2)],
                             closed=True, facecolor=color, edgecolor="none",
                             alpha=alpha, zorder=z))
    else:                                        # vertical
        d = 1.0 if y1 > y0 else -1.0
        by = y1 - d * head
        ax.add_patch(Polygon([(x0 - shaft / 2, y0), (x0 - shaft / 2, by),
                              (x0 - head_w / 2, by), (x0, y1),
                              (x0 + head_w / 2, by), (x0 + shaft / 2, by),
                              (x0 + shaft / 2, y0)],
                             closed=True, facecolor=color, edgecolor="none",
                             alpha=alpha, zorder=z))


def paste_slot(ax, x, y, w, h, hue, label, img=None):
    """A deliberately empty, labelled region for artwork we do not have yet.

    Distinct from tbd(): tbd() marks a NUMBER the manuscript has not settled,
    this marks a PICTURE someone still has to drop in. If ATOM_MAP_IMAGE exists
    the picture is placed instead and the slot disappears."""
    src = ATOM_MAP_IMAGE if img is None else img
    if src.exists():
        img = mpimg.imread(src)
        ih, iw = img.shape[0], img.shape[1]
        scale = min(w / iw, h / ih)
        dw, dh = iw * scale, ih * scale
        cx, cy = x + w / 2, y + h / 2
        _atom_map_dpi.append(iw / dw * 25.4)
        ax.imshow(img, extent=(cx - dw / 2, cx + dw / 2, cy - dh / 2, cy + dh / 2),
                  aspect="auto", interpolation="none", zorder=4)
        card(ax, cx - dw / 2, cy - dh / 2, dw, dh, face="none", edge=hue,
             lw=0.8, radius=1.5, alpha=0.5)
        return
    card(ax, x, y, w, h, face="none", edge=INK_MUTED, lw=1.0, ls=(0, (3, 3)),
         radius=2.0)
    text(ax, x + w / 2, y + h / 2, label, 12, color=INK_MUTED, ha="center",
         va="center")


def _blend(hue, alpha, base=SURFACE):
    """The flat RGB you get by painting `hue` at `alpha` over `base`.

    Needed because the join patch has to repaint the chip's interior exactly,
    and a second translucent layer would darken it."""
    h, b = to_rgb(hue), to_rgb(base)
    return tuple(b[k] + (h[k] - b[k]) * alpha for k in range(3))


# Source captures live in assets/, not in a figure OUTPUT directory --
# latex/figures/<row>/ is rewritten by scripts/regen_figures.py.
MOLECULE_DIR = Path(os.environ.get("NAR_MOLECULE_DIR",
                                   ROOT / "assets" / "molecules"))
_molecules_used: list[str] = []


# The mapping colour of glyoxylate's carboxylate carbon, read off that capture.
# CO2 was captured on its own, so its carbon came out plain black -- restating
# it in this colour is what makes the dashed trace's two ends read as one atom.
TRACED_C = "#1225cd"


def _has_colour(img, hexc, tol=0.02):
    return bool(np.all(np.abs(img[..., :3] - to_rgb(hexc)) < tol, axis=-1).any())


def _molecule(name, *, recolour=None):
    """Load a molecule capture and key out its background.

    The captures are opaque -- white surround, light-grey plate -- which would
    read as three grey tiles sitting on the figure surface. Everything above a
    luminance threshold becomes transparent, leaving only the bonds and atom
    labels. The atom colouring IS the atom mapping, so nothing that carries
    meaning is near that threshold.

    The keyed image is then cropped to its ink. The captures are square but the
    drawings are not -- glyoxylate fills 52% of its tile's height and CO2 only
    12% -- so an uncropped tile spends most of its allotted space on nothing.
    Cropping makes the drawn structure, rather than the screenshot, the thing
    that gets sized.

    `recolour` is ((fx0, fx1), colour) and restates one atom label in a
    different hue: every NEUTRAL pixel in that fractional column band is
    redrawn in `colour` at its original coverage, so antialiasing survives and
    the coloured bonds, which are not neutral, are left alone."""
    img = mpimg.imread(MOLECULE_DIR / f"{name}.png")
    if img.shape[2] == 3:
        img = np.dstack([img, np.ones(img.shape[:2])])
    img = img.copy()
    lum = img[..., :3].mean(axis=2)
    img[..., 3] = np.where(lum > 0.93, 0.0, img[..., 3])
    rows = np.where(img[..., 3].any(axis=1))[0]
    cols = np.where(img[..., 3].any(axis=0))[0]
    img = img[rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1]
    if recolour is not None:
        (fx0, fx1), colour = recolour
        w = img.shape[1]
        band = np.zeros(w, bool)
        band[int(fx0 * w):int(fx1 * w)] = True
        rgb = np.asarray(to_rgb(colour))
        m = (np.ptp(img[..., :3], axis=2) < 0.06) & band[None, :] \
            & (img[..., 3] > 0)
        cover = (1.0 - img[..., :3].mean(axis=2))[m][:, None]
        img[m, :3] = rgb * cover + (1.0 - cover)
    _molecules_used.append(name)
    return img


def molecule_image(ax, img, cx, cy, h, *, z=4):
    """Place a keyed molecule centred on (cx, cy) at height h. Returns width."""
    w = h * img.shape[1] / img.shape[0]
    ax.imshow(img, extent=(cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2),
              aspect="auto", interpolation="antialiased", zorder=z)
    return w


def reaction_equation(ax, x0, x1, cy, h, items, *, color=INK):
    """Lay out `items` left to right, centred in [x0, x1] at height cy.

    Items are ("mol", image), ("txt", string) or ("eq", None). Widths are
    measured first and the whole run is centred, so adding a species does not
    silently push the equation off its panel.

    The molecules share ONE scale in mm per source pixel, so a C=O in CO2 is
    drawn the same size as a C=O in tartronate semialdehyde -- the captures are
    all rendered at the same bond length, and giving each cropped drawing the
    same height instead would inflate the small species. `h` is the height
    budget for the tallest molecule; the scale is whichever of the height and
    width budgets binds first.

    Returns the placed boxes as [(kind, cx, cy, w, h), ...] so callers can
    anchor annotations to a specific ATOM inside a specific molecule without
    re-deriving the layout."""
    GAPS = {"mol": 2.6, "txt": 3.0, "eq": 3.2}
    EQ_W = 11.0
    gaps = [GAPS[k] for k, _ in items[:-1]]
    mols = [v for k, v in items if k == "mol"]
    fixed = sum(EQ_W if k == "eq" else len(v) * 2.6
                for k, v in items if k != "mol") + sum(gaps)
    scale = h / max(m.shape[0] for m in mols)
    if sum(m.shape[1] for m in mols) * scale > x1 - x0 - fixed:
        scale = (x1 - x0 - fixed) / sum(m.shape[1] for m in mols)
    heights, widths = [], []
    for kind, val in items:
        if kind == "mol":
            widths.append(val.shape[1] * scale)
            heights.append(val.shape[0] * scale)
        else:
            widths.append(EQ_W if kind == "eq" else len(val) * 2.6)
            heights.append(h)
    total = sum(widths) + sum(gaps)
    cx = (x0 + x1) / 2 - total / 2
    placed = []
    for i, (kind, val) in enumerate(items):
        w, hi = widths[i], heights[i]
        if kind == "mol":
            molecule_image(ax, val, cx + w / 2, cy, hi)
        elif kind == "eq":
            equilibrium(ax, cx + w / 2, cy, w=EQ_W, color=color, lw=2.4, ms=13,
                        sep=1.7)
        else:
            text(ax, cx + w / 2, cy, val, 14, weight="bold", ha="center",
                 va="center", color=color, check=False)
        placed.append((kind, cx + w / 2, cy, w, hi))
        cx += w + (gaps[i] if i < len(gaps) else 0.0)
    return placed


def cartoon_histogram(ax, x, y, w, h, hue, values, *, gap=0.22):
    """A shape, not a chart: bars only, no axis, no ticks, no labels.

    Used where the figure needs to say "this distribution is skewed with a long
    tail" without asserting a number. Anything quantitative belongs in a real
    chart with a scale."""
    n = len(values)
    bw = (w - gap * (n - 1)) / n
    peak = max(values) or 1.0
    for i, v in enumerate(values):
        bh = max(h * v / peak, 0.35)
        _rrect(ax, x + i * (bw + gap), y, bw, bh, hue, alpha=0.85, z=4, r=0.35)
    ax.plot([x, x + w], [y, y], color=hue, lw=0.9, alpha=0.55, zorder=4,
            solid_capstyle="butt")


def concept_server_hub(ax):
    """The hand-sketched layout: the database as a hub.

    NO COUNTS -- a structure diagram; the counts live in the other concepts.

    SPACING IS DERIVED. Every arrow-to-object clearance is GAP and every row
    group is centred on HUB_CY. Hand-typed offsets had drifted to four
    different values before this was pulled out.

    COLOUR. Only the four validated stage hues are used (slots 1/2/3/7, all
    pairs PASS). The three inputs deliberately SHARE one hue: they are one
    family -- structures compose compounds, which compose reactions -- and the
    semicircular arcs carry that sequence. Giving them three hues needed a
    second blue that failed the normal-vision floor against the first
    (ΔE 14.7), and it also consumed the orange that the experimental arrow
    needs to stand apart from the three computational ones.
    """
    GAP = 2.0
    HUB_CY = 66.0
    IN_HUE = STAGES["mol"]            # the three database-content inputs
    EXP_HUE = STAGES["rxn"]           # OpenTECR, and the Experimental box
    COMP_HUE = STAGES["thermo"]       # the three estimators, and Computational
    SRV_X, SRV_W, SRV_H = 58.0, 34.0, 42.0
    SRV_Y = HUB_CY - SRV_H / 2

    def rows(n, pitch):
        return [HUB_CY + (n - 1) / 2 * pitch - i * pitch for i in range(n)]

    # ---- left: structures -> compounds -> reactions
    # Each box extrudes an arrow from its BOTTOM edge into the box below, so
    # the "composes" sequence is carried by the boxes themselves. Pitch 17 on
    # an 11 mm box leaves a 6 mm gap, which the 1.8 + 3.4 mm arrow clears with
    # 0.8 to spare. The top box then reaches y 88.5, just under the title.
    CHIP_X, CHIP_W, CHIP_H = 2.0, 47.0, 11.0
    in_rows = [("STRUCTURES", "chain"), ("COMPOUNDS", "ring"),
               ("REACTIONS", "rxn")]
    ys = rows(3, 17.0)
    for i, ((lab, glyph), cy) in enumerate(zip(in_rows, ys)):
        contain(CHIP_X, cy - CHIP_H / 2, CHIP_W, CHIP_H, f"{lab} chip")
        card_down_arrow(ax, CHIP_X, cy - CHIP_H / 2, CHIP_W, CHIP_H, IN_HUE,
                        arrow=i < len(in_rows) - 1)
        gx = CHIP_X + 6.4
        if glyph == "ring":
            molecule_glyph(ax, gx, cy, 2.7, IN_HUE, kind="ring")
        elif glyph == "chain":
            molecule_glyph(ax, gx, cy, 2.3, IN_HUE, kind="chain")
        else:
            # Laid out left to right by each glyph's REAL extent, not on a
            # constant pitch. For the same nominal size the chain reaches
            # 1.6*s either side of its centre while the ring reaches 0.87*s
            # left and 0.90*s right (its substituent), so an even pitch buried
            # the equilibrium's right arrowhead 0.9 mm inside the chain.
            s_r, s_c, eq_w, sep = 1.4, 1.1, 2.4, 0.9
            rx = CHIP_X + 1.4 + 0.866 * s_r
            ex = rx + 0.9 * s_r + sep + eq_w / 2
            nx = ex + eq_w / 2 + sep + 1.6 * s_c
            molecule_glyph(ax, rx, cy, s_r, IN_HUE, kind="ring")
            equilibrium(ax, ex, cy, w=eq_w)
            molecule_glyph(ax, nx, cy, s_c, IN_HUE, kind="chain")
        text(ax, CHIP_X + 13.0, cy, lab, 12, weight="bold")
        block_arrow(ax, CHIP_X + CHIP_W + GAP, cy, SRV_X - GAP, cy, IN_HUE,
                    shaft=4.2, head=4.6, head_w=7.0, alpha=ARROW_A)

    # ---- the hub
    server_rack(ax, SRV_X, SRV_Y, SRV_W, SRV_H)

    # ---- bottom: atom mapping
    ATOM_X, ATOM_W, ATOM_Y, ATOM_H = 17.0, 182.0, 1.0, 31.0
    block_arrow(ax, SRV_X + SRV_W / 2, SRV_Y - GAP, SRV_X + SRV_W / 2,
                ATOM_Y + ATOM_H + GAP, STAGES["atom"], shaft=4.0, head=5.0,
                head_w=11.0, alpha=ARROW_A)
    contain(ATOM_X, ATOM_Y, ATOM_W, ATOM_H, "atom-mapping box")
    card(ax, ATOM_X, ATOM_Y, ATOM_W, ATOM_H, face=STAGES["atom"], edge="none",
         alpha=0.09, radius=2.5)
    card(ax, ATOM_X, ATOM_Y, ATOM_W, ATOM_H, face="none", edge=STAGES["atom"],
         lw=1.4, radius=2.5)
    # Label on the LEFT so the equation gets the panel's full height; the
    # molecule captures are square and were being squeezed under a header.
    text(ax, ATOM_X + 3.5, ATOM_Y + ATOM_H / 2, "ATOM\nMAPPING", 12,
         weight="bold", color=STAGES["atom"], va="center", linespacing=1.25)
    # Separate the name from the illustration with a rule, centred in the gap
    # between them: "MAPPING" is the wider of the two lines at 22.4 mm so the
    # name ends at 45.9, and the centred equation's ink starts at about 52.
    ax.add_patch(Rectangle((ATOM_X + 28.5, ATOM_Y + 4.0), 0.7, ATOM_H - 8.0,
                           facecolor=STAGES["atom"], edgecolor="none",
                           alpha=0.45, zorder=3))
    if ATOM_MAP_IMAGE.exists():
        paste_slot(ax, ATOM_X + 30.0, ATOM_Y + 3.0, ATOM_W - 34.0,
                   ATOM_H - 6.0, STAGES["atom"], "atom mapping")
    else:
        # glyoxylate <=> CO2 + tartronate semialdehyde. The atom colouring in
        # the captures is the mapping itself, which is the point of the panel.
        # The stoichiometric 2 is dropped: this panel is about which atom goes
        # where, and the coefficient was reading as part of the panel's title.
        glyoxylate = _molecule("Glyoxalate")
        if not _has_colour(glyoxylate, TRACED_C):
            raise SystemExit(
                f"the glyoxylate capture no longer contains {TRACED_C}: "
                "re-read the traced carbon's colour before trusting the CO2 "
                "recolour, or the two ends of the arrow will disagree.")
        # CO2's C glyph is columns 111..136 of a 248-wide crop, with empty
        # columns either side of it, so this band takes the letter and no bond.
        placed = reaction_equation(
            ax, ATOM_X + 30.0, ATOM_X + ATOM_W - 4.0, ATOM_Y + ATOM_H / 2,
            28.0,
            [("mol", glyoxylate), ("eq", None),
             ("mol", _molecule("CO2", recolour=((0.443, 0.557), TRACED_C))),
             ("txt", "+"),
             ("mol", _molecule("Tartronate Semialdehyde"))])
        # Trace one atom across the reaction: the carboxylate carbon of
        # glyoxylate is the one released as CO2. Anchors are fractions of each
        # capture, read off the artwork and then mapped onto the cropped tiles
        # -- glyoxylate's carboxylate C at (0.63, 0.37), CO2's C near the middle
        # of its one-line drawing. The matching dark-green oxygens in both
        # captures confirm the pair.
        def _atom(idx, fx, fy):
            _k, mx, my, mw, mh = placed[idx]
            return mx + (fx - 0.5) * mw, my + (fy - 0.5) * mh
        a = _atom(0, 0.630, 0.368)
        b = _atom(2, 0.508, 0.457)
        # Solve the arc3 rad rather than hand-tuning it. The two anchors sit
        # at different heights inside their molecules, so a fixed rad puts the
        # crown wherever the chord happens to fall -- at -0.42 it grazed
        # glyoxylate's carbonyl oxygen. arc3 places its control point at
        # (mid + rad*dy, mid - rad*dx), which puts the curve's apex at
        # midy - rad*dx/2, so pinning the apex a fixed clearance under the
        # panel's top edge inverts to one expression. The trace then passes
        # cleanly over both structures and over the equilibrium arrows.
        apex = ATOM_Y + ATOM_H - 3.4
        rad = 2.0 * ((a[1] + b[1]) / 2 - apex) / (b[0] - a[0])
        conn = f"arc3,rad={rad:.4f}"
        # Shaft and head are two patches over the same curve. A dashed
        # FancyArrowPatch applies its dash pattern to the ARROWHEAD's outline
        # as well as the shaft, which chews the triangle into a blob with a
        # spur hanging off one side; the head has to be its own solid,
        # unstroked patch to come out crisp. The shaft then stops at the head's
        # base, which for this style is mutation_scale * 0.4 points long.
        head_pt = 15 * 0.4
        ax.add_patch(FancyArrowPatch(
            a, b, connectionstyle=conn, arrowstyle="-", linewidth=1.7,
            linestyle=(0, (2.9, 2.0)), color=STAGES["atom"], zorder=6,
            shrinkA=8.0, shrinkB=11.0 + head_pt, capstyle="round"))
        ax.add_patch(FancyArrowPatch(
            a, b, connectionstyle=conn, arrowstyle="-|>", mutation_scale=15,
            linewidth=0.0, facecolor=STAGES["atom"], edgecolor="none",
            zorder=6, shrinkA=8.0, shrinkB=11.0))

    # ---- right: one experimental source, three computational ones, merged
    # Every estimator arrow has to leave the rack ALONGSIDE it, so all four
    # rows live inside the server's own 45..87 span -- the bottom one used to
    # sit at 42.4, below the rack's floor, and read as coming from nowhere.
    # "Group Contrib." is 35.6 mm at 12 pt and still sets the box start: the
    # labels run 104..139.6, so BOX_X cannot come in past ~141.
    LBL_X = SRV_X + SRV_W + GAP
    BOX_X, BOX_W = 140.0, 44.0
    EXP_Y, EXP_H = 75.0, 12.0
    COMP_Y, COMP_H = 54.5, 12.0        # same size as Experimental
    PAN_X, PAN_W, PAN_Y, PAN_H = 192.0, 21.0, 45.0, 42.0
    PAN_HUE = INK_2

    block_arrow(ax, LBL_X, EXP_Y + EXP_H / 2, BOX_X - GAP, EXP_Y + EXP_H / 2,
                EXP_HUE, shaft=2.8, head=5.2, alpha=ARROW_A)
    text(ax, LBL_X, EXP_Y + EXP_H / 2 + 4.4, "OpenTECR", 12,
         weight="bold", color=EXP_HUE)

    # The outer two estimators curve into the middle one, which stays straight
    # and full length, so the box they enter can be the same small size as
    # Experimental instead of being stretched to catch three heads.
    #
    # Label placement is what makes the curves possible. Each curve crosses the
    # band between its own row and the middle row, so that band has to be clear
    # of type from where the curve starts. Putting "Group Contrib." BELOW its
    # shaft empties the lower band, which lets both curves start at the same x
    # (124, just past "dGPredictor" at 123.0) and stay symmetric. Each label is
    # still 1.0 mm off its own shaft and at least 1.4 from any other.
    comp_rows = [70.0, 60.5, 51.0]
    for lab, cy, dy in zip(["eQuilibrator", "dGPredictor", "Group Contrib."],
                           comp_rows, (4.56, 4.56, -4.56)):
        text(ax, LBL_X, cy + dy, lab, 12, weight="bold")
    confluence_arrow(ax, LBL_X, comp_rows, COMP_Y + COMP_H / 2, 124.0, 133.0,
                     BOX_X - GAP, COMP_HUE, alpha=ARROW_A)

    # Each kind gets its own box, and each box then feeds the ONE panel where
    # the estimates are pooled -- the grades are assigned off the merged
    # distribution, not off either source alone.
    # Icon then label, both LEFT aligned rather than centred as a group, so
    # the two rows line up with each other even though "Experimental" is 3 mm
    # shorter than "Computational". 1.5 pad + 3.8 icon + 1.6 gap + 35.6 label
    # puts the label's right edge 1.5 mm inside a 44 mm box, which clears the
    # container check's 1.2. That is the floor: the label cannot shrink.
    ICON_H, ICON_GAP, ICON_PAD = 3.8, 1.6, 1.5
    for by, bh, lab, hue, glyph in [
            (EXP_Y, EXP_H, "Experimental", EXP_HUE, flask_glyph),
            (COMP_Y, COMP_H, "Computational", COMP_HUE, chip_glyph)]:
        contain(BOX_X, by, BOX_W, bh, f"{lab} box")
        card(ax, BOX_X, by, BOX_W, bh, face=hue, edge="none", alpha=FILL_A,
             radius=2.5)
        card(ax, BOX_X, by, BOX_W, bh, face="none", edge=hue, lw=1.5,
             radius=2.5)
        glyph(ax, BOX_X + ICON_PAD + ICON_H / 2, by + bh / 2, ICON_H, hue)
        text(ax, BOX_X + ICON_PAD + ICON_H + ICON_GAP, by + bh / 2, lab, 12,
             weight="bold", va="center", color=hue)
        block_arrow(ax, BOX_X + BOX_W + GAP, by + bh / 2, PAN_X - GAP,
                    by + bh / 2, hue, shaft=2.4, head=3.2, head_w=5.6,
                    alpha=ARROW_A)

    # ---- the merged panel: the pooled uncertainty distribution.
    # No tint behind the bars -- a grey wash made the chart look switched off,
    # and the panel does not need a fill to read as a container when it has a
    # border, a tag and content.
    contain(PAN_X, PAN_Y, PAN_W, PAN_H, "dG panel")
    card(ax, PAN_X, PAN_Y, PAN_W, PAN_H, face=SURFACE, edge=PAN_HUE, lw=1.5,
         radius=2.5)
    TAG_H = 8.5
    tag_y = PAN_Y + PAN_H - 3.0 - TAG_H
    contain(PAN_X + 3.0, tag_y, PAN_W - 6.0, TAG_H, "dG tag")
    card(ax, PAN_X + 3.0, tag_y, PAN_W - 6.0, TAG_H, face=SURFACE,
         edge=PAN_HUE, lw=1.2, radius=2.0)
    text(ax, PAN_X + PAN_W / 2, tag_y + TAG_H / 2, "\u0394G", 12,
         weight="bold", ha="center", va="center", color=INK)
    # Cartoon of eQuilibrator's reported-uncertainty distribution. A real
    # image at BAR_CHART_IMAGE overrides it.
    if BAR_CHART_IMAGE.exists():
        paste_slot(ax, PAN_X + 3.0, PAN_Y + 3.0, PAN_W - 6.0,
                   PAN_H - 9.0 - TAG_H, PAN_HUE, "bar chart",
                   img=BAR_CHART_IMAGE)
    else:
        cartoon_histogram(ax, PAN_X + 3.0, PAN_Y + 3.0, PAN_W - 6.0,
                          PAN_H - 9.0 - TAG_H, COMP_HUE, EQ_SIGMA_CARTOON)

    # ---- classification
    GR_X, GR_W, GR_H = 224.0, 23.0, 13.0
    grades = [("GOLD", GRADE_RAMP[0]), ("SILVER", GRADE_RAMP[1]),
              ("BRONZE", GRADE_RAMP[2])]
    stack_cy = PAN_Y + PAN_H / 2
    for (lab, col), cy in zip(grades,
                              [stack_cy + 14.0, stack_cy, stack_cy - 14.0]):
        block_arrow(ax, PAN_X + PAN_W + GAP, cy, GR_X - GAP, cy, col,
                    shaft=2.6, head=4.6, alpha=ARROW_A)
        contain(GR_X, cy - GR_H / 2, GR_W, GR_H, f"{lab} chip")
        card(ax, GR_X, cy - GR_H / 2, GR_W, GR_H, face=col, edge="none",
             alpha=0.20, radius=2.5)
        card(ax, GR_X, cy - GR_H / 2, GR_W, GR_H, face="none", edge=col,
             lw=1.5, radius=2.5)
        text(ax, GR_X + GR_W / 2, cy, lab, 12, weight="bold", color=col,
             ha="center")
    # Centred on the column. "Classification" is 32.3 mm against a 23 mm
    # column, so it overhangs both sides by 4.6; the stack sits at 224 rather
    # than 229 precisely so the right overhang stays on the canvas.
    text(ax, GR_X + GR_W / 2, stack_cy + 14.0 + GR_H / 2 + 4.0,
         "Classification", 12, weight="bold", ha="center", color=INK)

CONCEPTS = {
    "flow_pipeline": concept_flow_pipeline,
    "one_reaction_four_lenses": concept_one_reaction,
    "before_after_ledger": concept_before_after,
    "direction_funnel": concept_direction_funnel,
    "atom_trace_spine": concept_atom_spine,
    "server_hub": concept_server_hub,
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
    _check_text_overlaps(fig, ax)
    _check_text_in_containers(fig, ax)

    stem = out / f"graphical_abstract_{args.concept}"
    # dpi matters even for the vector formats: matplotlib resamples embedded
    # IMAGES to the output device resolution, and at the default 100 dpi the UI
    # capture came out at 100 dpi -- under NAR's 300 dpi floor for colour
    # half-tones. Vector geometry is unaffected by this.
    fig.savefig(f"{stem}.pdf", facecolor=SURFACE, dpi=DPI_VECTOR)
    # SVG as well: pdf.fonttype 42 embeds SUBSET CID fonts, which Illustrator
    # can only edit as text if the same font is installed locally. svg.fonttype
    # "none" writes <text> with a font-family name and embeds nothing, so the
    # text stays fully editable on any machine. Same vector geometry.
    fig.savefig(f"{stem}.svg", facecolor=SURFACE, dpi=DPI_VECTOR)
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
        fh.write(f"ui_capture\t{'present' if SCREENSHOT.exists() else 'MISSING'}"
                 f"\t{SCREENSHOT}\n")
        if _capture_dpi:
            d = _capture_dpi[0]
            fh.write(f"ui_capture_dpi\t{d:.0f}\tsource pixels at the placed "
                     f"size; {'clears' if d >= 300 else 'BELOW'} NAR's 300 dpi "
                     "floor for colour half-tones\n")
        fh.write(f"atom_mapping_reaction\t2 glyoxylate <=> CO2 + tartronate "
                 f"semialdehyde\tmolecule captures: "
                 f"{', '.join(_molecules_used) if _molecules_used else 'none'}\n")
        fh.write(f"eq_sigma_cartoon\t{len(EQ_SIGMA_CARTOON)} bars\tCARTOON, not "
                 "data: silhouette traced from panel C chart 1 of "
                 "figures/main_figures_draft.pdf, no axis or scale drawn\n")
        fh.write(f"atom_map_image\t{'present' if ATOM_MAP_IMAGE.exists() else 'EMPTY SLOT'}"
                 f"\t{ATOM_MAP_IMAGE} -- server_hub reserves a paste slot for this\n")
        if _atom_map_dpi:
            d = _atom_map_dpi[0]
            fh.write(f"atom_map_dpi\t{d:.0f}\tsource pixels at the placed size; "
                     f"{'clears' if d >= 300 else 'BELOW'} NAR's 300 dpi floor\n")
        fh.write(f"ui_capture_crop\t{SCREENSHOT_CROP}\t"
                 "fractions x0,x1,y0,y1 -- crops away the nav bar, which "
                 "carries the logo NAR bans and the signed-in username\n")
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
