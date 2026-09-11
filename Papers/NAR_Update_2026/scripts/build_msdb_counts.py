#!/usr/bin/env python3
"""Extract the counts the graphical-abstract bar charts need, from a ModelSEED
snapshot, into data/msdb_counts.tsv.

WHY A SEPARATE SCRIPT
---------------------
Scanning 56k reaction records and 46k compound records takes a few seconds. The
figure generator runs five times (once per concept), so the scan is done once
here and cached as a TSV that carries its own provenance header.

DATA
----
Reads a snapshot of upstream ModelSEED `dev`, never the working checkout --
`/scratch/ctaylor/ModelSEEDDatabase` sits on a feature branch that carries local
edits. Point MSDB_ROOT at the snapshot; see the fresh-msdb-data procedure for
how to cut one:

    git archive origin/dev Biochemistry Scripts | tar x -C /scratch/ctaylor/tmp/devsnap_<sha>

Every count is a STORED field, never recomputed:
  * balanced        -- "OK" in reaction["status"]
  * structure       -- compound has a smiles or inchikey
  * thermodynamics  -- reaction["thermodynamics"] is {source: [dG, err, op, ...]}
  * evidence        -- reaction["thermo-evidence"] is {assessment, grade, source}
  * atom mapping    -- reaction["atom_mapping"] is {confidence, data}
  * direction       -- reaction["reversibility"] in > < = ?

2020 BASELINES
--------------
Growth figures come from Seaver et al. 2020 Tables 2-3 as transcribed in
../data/snapshot_2026-07-29.md. The eQuilibrator baseline is read from
Thermodynamics/eQuilibrator/eQuilibrator-2020_reactions.tsv, whose own header
describes it as the values "stored in ModelSEED before the 2026-08
regeneration". That is a proxy for the 2020 release, not a re-run of the 2020
pipeline, and it is labelled eq_prerefresh_reactions rather than "2020" so
nothing downstream can quietly promote it.

USAGE
-----
    MSDB_ROOT=/scratch/ctaylor/tmp/devsnap_<sha> python3 scripts/build_msdb_counts.py
"""
from __future__ import annotations

import collections
import glob
import json
import os
import subprocess
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MSDB = Path(os.environ["MSDB_ROOT"])
OUT = ROOT / "data" / "msdb_counts.tsv"

# Seaver et al. 2020, Tables 2-3 (via ../data/snapshot_2026-07-29.md).
BASE_2020 = {
    "compounds": 33992, "compounds_with_structure": 28120,
    "reactions": 36193, "reactions_balanced": 25457,
    "reversible": 18399, "functional_reactions": 21403,
    "biolog_functional": 355, "biolog_total": 390,
}


def dev_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", "/scratch/ctaylor/ModelSEEDDatabase",
             "rev-parse", "--short", "origin/dev"], text=True).strip()
    except Exception:
        return "unknown"


def main() -> int:
    rows: list[tuple[str, int, str]] = []
    thermo = collections.Counter()
    n_est = collections.Counter()
    n_est_bal = collections.Counter()
    grade = collections.Counter()
    direction = collections.Counter()
    direction_bal = collections.Counter()
    atom = collections.Counter()
    nrxn = balanced = 0

    for f in sorted(glob.glob(str(MSDB / "Biochemistry" / "reaction_*.json"))):
        for r in json.load(open(f)):
            nrxn += 1
            is_bal = "OK" in (r.get("status") or "")
            balanced += is_bal
            est = 0
            for src, val in (r.get("thermodynamics") or {}).items():
                if val and val[0] is not None:
                    thermo[src] += 1
                    est += 1
            n_est[est] += 1
            if is_bal:
                n_est_bal[est] += 1
            direction[r.get("reversibility") or "?"] += 1
            if is_bal:
                direction_bal[r.get("reversibility") or "?"] += 1
            a = r.get("atom_mapping") or {}
            if a.get("data"):
                atom[a.get("confidence") or "unlabelled"] += 1
            e = r.get("thermo-evidence") or {}
            if e.get("grade"):
                grade[e["grade"]] += 1
            if e.get("source"):
                thermo.setdefault("TECRDB", 0)
                if e["source"] == "TECRDB":
                    thermo["TECRDB"] += 1

    ncpd = withstruct = 0
    for f in sorted(glob.glob(str(MSDB / "Biochemistry" / "compound_*.json"))):
        for c in json.load(open(f)):
            ncpd += 1
            if c.get("smiles") or c.get("inchikey"):
                withstruct += 1

    eq2020 = MSDB / "Biochemistry/Thermodynamics/eQuilibrator/eQuilibrator-2020_reactions.tsv"
    n_eq_pre = 0
    if eq2020.exists():
        with open(eq2020) as fh:
            n_eq_pre = sum(1 for ln in fh
                           if ln.strip() and not ln.startswith(("#", "reaction_id")))

    rows += [("compounds", ncpd, "2026"),
             ("compounds_with_structure", withstruct, "2026"),
             ("reactions", nrxn, "2026"),
             ("reactions_balanced", balanced, "2026")]
    for k, v in BASE_2020.items():
        rows.append((k, v, "2020"))
    for src, v in sorted(thermo.items()):
        rows.append((f"thermo_reactions:{src}", v, "2026"))
    rows.append(("thermo_reactions:eQuilibrator_prerefresh", n_eq_pre,
                 "pre-2026-refresh"))
    for k in sorted(n_est):
        rows.append((f"thermo_estimates_per_reaction:{k}", n_est[k], "2026"))
    for k in sorted(n_est_bal):
        rows.append((f"thermo_estimates_per_balanced_reaction:{k}",
                     n_est_bal[k], "2026"))
    for k, v in sorted(grade.items()):
        rows.append((f"thermo_evidence_grade:{k}", v, "2026"))
    for k in (">", "<", "=", "?"):
        rows.append((f"direction:{k}", direction.get(k, 0), "2026"))
        rows.append((f"direction_balanced:{k}", direction_bal.get(k, 0), "2026"))
    for k, v in sorted(atom.items()):
        rows.append((f"atom_mapping:{k}", v, "2026"))
    rows.append(("atom_mapping:total", sum(atom.values()), "2026"))

    sha = dev_sha()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        fh.write(f"# ModelSEED dev @ {sha} ({MSDB}), read {date.today()}\n")
        fh.write("# 2020 rows are Seaver et al. 2020 Tables 2-3, not recomputed.\n")
        fh.write("# All 2026 rows are STORED fields, not recomputed locally.\n")
        fh.write("key\tvalue\tvintage\n")
        for k, v, vintage in rows:
            fh.write(f"{k}\t{v}\t{vintage}\n")
    print(f"ModelSEED dev @ {sha}")
    print(f"wrote {OUT}  ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
