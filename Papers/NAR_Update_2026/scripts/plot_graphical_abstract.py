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
import matplotlib.image as mpimg
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


def equilibrium(ax, cx, cy, w=4.4, color=INK_2):
    """Reversible-reaction symbol as two offset half-arrows."""
    for dy, x0, x1 in ((0.7, -1, 1), (-0.7, 1, -1)):
        ax.add_patch(FancyArrowPatch((cx + x0 * w / 2, cy + dy),
                                     (cx + x1 * w / 2, cy + dy),
                                     arrowstyle="-|>", mutation_scale=5,
                                     linewidth=1.0, color=color, zorder=4,
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
    text(ax, W / 2, H - 9.0, "ModelSEED Biochemistry Database  —  2026 update",
         16, weight="bold", ha="center")

    GAP = 2.0
    HUB_CY = 62.0
    IN_HUE = STAGES["mol"]            # the three database-content inputs
    EXP_HUE = STAGES["rxn"]           # OpenTECR, and the Experimental box
    COMP_HUE = STAGES["thermo"]       # the three estimators, and Computational
    SRV_X, SRV_W, SRV_H = 78.0, 34.0, 44.0
    SRV_Y = HUB_CY - SRV_H / 2

    def rows(n, pitch):
        return [HUB_CY + (n - 1) / 2 * pitch - i * pitch for i in range(n)]

    # ---- left: structures -> compounds -> reactions
    CHIP_X, CHIP_W, CHIP_H = 10.0, 48.0, 14.0
    in_rows = [("STRUCTURES", "chain"), ("COMPOUNDS", "ring"),
               ("REACTIONS", "rxn")]
    ys = rows(3, 14.0)
    for (lab, glyph), cy in zip(in_rows, ys):
        contain(CHIP_X, cy - CHIP_H / 2, CHIP_W, CHIP_H, f"{lab} chip")
        card(ax, CHIP_X, cy - CHIP_H / 2, CHIP_W, CHIP_H, face=IN_HUE,
             edge="none", alpha=FILL_A, radius=2.5)
        card(ax, CHIP_X, cy - CHIP_H / 2, CHIP_W, CHIP_H, face="none",
             edge=IN_HUE, lw=1.4, radius=2.5)
        gx = CHIP_X + 6.4
        if glyph == "ring":
            molecule_glyph(ax, gx, cy, 2.7, IN_HUE, kind="ring")
        elif glyph == "chain":
            molecule_glyph(ax, gx, cy, 2.3, IN_HUE, kind="chain")
        else:
            molecule_glyph(ax, gx - 2.8, cy, 1.6, IN_HUE, kind="ring")
            equilibrium(ax, gx, cy, w=2.2)
            molecule_glyph(ax, gx + 2.8, cy, 1.6, IN_HUE, kind="chain")
        text(ax, CHIP_X + 13.0, cy, lab, 12, weight="bold")
        block_arrow(ax, CHIP_X + CHIP_W + GAP, cy, SRV_X - GAP, cy, IN_HUE,
                    shaft=4.2, head=7.0, alpha=ARROW_A)

    # Semicircular influence arcs down the LEFT edge: structures compose
    # compounds, compounds compose reactions. They bow out into the 10 mm
    # margin and land on the upper third of the box below, so they read as
    # "feeds into" rather than as a second data flow.
    for cy_up, cy_dn in zip(ys[:-1], ys[1:]):
        ax.add_patch(FancyArrowPatch(
            (CHIP_X, cy_up - CHIP_H / 4), (CHIP_X, cy_dn + CHIP_H / 6),
            connectionstyle="arc3,rad=0.95", arrowstyle="-|>",
            mutation_scale=13, linewidth=1.8, color=INK_2, zorder=5,
            shrinkA=1.0, shrinkB=1.0))

    # ---- the hub
    server_rack(ax, SRV_X, SRV_Y, SRV_W, SRV_H)

    # ---- bottom: atom mapping
    ATOM_X, ATOM_W, ATOM_Y, ATOM_H = 45.0, 157.0, 2.0, 23.0
    block_arrow(ax, SRV_X + SRV_W / 2, SRV_Y - GAP, SRV_X + SRV_W / 2,
                ATOM_Y + ATOM_H + GAP, STAGES["atom"], shaft=5.0, head=5.5,
                head_w=13.0, alpha=ARROW_A)
    contain(ATOM_X, ATOM_Y, ATOM_W, ATOM_H, "atom-mapping box")
    card(ax, ATOM_X, ATOM_Y, ATOM_W, ATOM_H, face=STAGES["atom"], edge="none",
         alpha=0.09, radius=2.5)
    card(ax, ATOM_X, ATOM_Y, ATOM_W, ATOM_H, face="none", edge=STAGES["atom"],
         lw=1.4, radius=2.5)
    text(ax, ATOM_X + 5.0, ATOM_Y + ATOM_H - 4.8, "ATOM MAPPING", 12,
         weight="bold", color=STAGES["atom"])
    paste_slot(ax, ATOM_X + 5.0, ATOM_Y + 3.0, ATOM_W - 10.0, ATOM_H - 11.5,
               STAGES["atom"], "paste atom-mapping capture here")

    # ---- right: one experimental source, three computational ones
    # "Group contribution" is 46.2 mm at 12 pt and sets where the boxes start.
    BOX_X, BOX_W = 162.0, 40.0
    EXP_Y, EXP_H = 71.0, 14.0
    COMP_Y, COMP_H = 34.0, 35.0

    block_arrow(ax, SRV_X + SRV_W + GAP, EXP_Y + EXP_H / 2, BOX_X - GAP,
                EXP_Y + EXP_H / 2, EXP_HUE, shaft=3.6, head=6.4, alpha=ARROW_A)
    text(ax, SRV_X + SRV_W + GAP, EXP_Y + EXP_H / 2 + 4.6, "OpenTECR", 12,
         weight="bold", color=EXP_HUE)

    comp_rows = [COMP_Y + COMP_H * f for f in (0.86, 0.55, 0.24)]
    for lab, cy in zip(["eQuilibrator", "dGPredictor", "Group contribution"],
                       comp_rows):
        block_arrow(ax, SRV_X + SRV_W + GAP, cy, BOX_X - GAP, cy, COMP_HUE,
                    shaft=3.6, head=6.4, alpha=ARROW_A)
        text(ax, SRV_X + SRV_W + GAP, cy + 4.6, lab, 12, weight="bold")

    for (bx, by, bw, bh, lab, hue) in [
            (BOX_X, EXP_Y, BOX_W, EXP_H, "Experimental", EXP_HUE),
            (BOX_X, COMP_Y, BOX_W, COMP_H, "Computational", COMP_HUE)]:
        contain(bx, by, bw, bh, f"{lab} box")
        card(ax, bx, by, bw, bh, face=hue, edge="none", alpha=FILL_A,
             radius=2.5)
        card(ax, bx, by, bw, bh, face="none", edge=hue, lw=1.5, radius=2.5)
        text(ax, bx + bw / 2, by + bh - 5.0, lab, 12, weight="bold",
             ha="center", color=hue)
    # room for a bar chart of per-estimator coverage, pasted later
    paste_slot(ax, BOX_X + 4.0, COMP_Y + 3.0, BOX_W - 8.0, COMP_H - 12.0,
               COMP_HUE, "bar chart", img=BAR_CHART_IMAGE)

    # A rail so all three grades read as drawing on BOTH boxes, rather than
    # gold appearing to come from the experimental box it happens to sit beside.
    RAIL_X = BOX_X + BOX_W + GAP
    card(ax, RAIL_X, COMP_Y + 2.0, 1.8, (EXP_Y + EXP_H) - COMP_Y - 4.0,
         face=INK_MUTED, edge="none", radius=0.9, alpha=0.55)

    # ---- classification
    GR_X, GR_W, GR_H = 218.0, 32.0, 14.0
    grades = [("GOLD", GRADE_RAMP[0]), ("SILVER", GRADE_RAMP[1]),
              ("BRONZE", GRADE_RAMP[2])]
    stack_cy = (COMP_Y + EXP_Y + EXP_H) / 2
    for (lab, col), cy in zip(grades, [stack_cy + 14.0, stack_cy, stack_cy - 14.0]):
        block_arrow(ax, RAIL_X + 1.8 + GAP, cy, GR_X - GAP, cy, col,
                    shaft=3.2, head=5.2, alpha=ARROW_A)
        contain(GR_X, cy - GR_H / 2, GR_W, GR_H, f"{lab} chip")
        card(ax, GR_X, cy - GR_H / 2, GR_W, GR_H, face=col, edge="none",
             alpha=0.20, radius=2.5)
        card(ax, GR_X, cy - GR_H / 2, GR_W, GR_H, face="none", edge=col,
             lw=1.5, radius=2.5)
        text(ax, GR_X + GR_W / 2, cy, lab, 12, weight="bold", color=col,
             ha="center")
    text(ax, GR_X + GR_W / 2, stack_cy + 14.0 + GR_H / 2 + 5.0,
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
