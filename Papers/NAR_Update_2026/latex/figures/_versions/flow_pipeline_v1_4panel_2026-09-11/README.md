# flow_pipeline v1 — four panels, frozen 2026-09-11

The state the flow concept was approved at, before the RANKING panel and the
radial backgrounds were added. Kept so the two can be compared side by side.

Four panels: MOLECULES / REACTIONS / THERMODYNAMICS / ATOM MAPPING, flush
headers, vertical linear gradients, no connecting arrows, atom mapping as a stat
block. All values real, from ModelSEED dev @ d64fdc6f.

`plot_graphical_abstract.py.snapshot` is the exact generator that produced these
files. To rebuild them:

    python3 <this dir>/plot_graphical_abstract.py.snapshot --concept flow_pipeline \
        --out-subdir _versions/flow_pipeline_v1_4panel_2026-09-11/rebuilt

The live generator has moved on; do not expect `regen_figures.py` to reproduce
this directory.
