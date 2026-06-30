#!/usr/bin/env python
"""Run EVERY reversibility heuristic on EVERY reaction and dump the result.

The normal cascade (``Estimate_Reaction_Reversibility.py``) short-circuits at the
first heuristic that fires, so the stored report only shows the *winning* rule.
This script instead evaluates **all** heuristics on **every** reaction (no
short-circuit) via ``reversibility_heuristics.evaluate_all_heuristics`` and saves,
per reaction, what each heuristic returned -- to both CSV and JSON.

Usage::

    ./Dump_Reaction_Heuristic_Outputs.py            # EQ-level energy (default)
    ./Dump_Reaction_Heuristic_Outputs.py GC         # GC-level
    ./Dump_Reaction_Heuristic_Outputs.py DGP        # dGPredictor-level
    ./Dump_Reaction_Heuristic_Outputs.py ''         # unfiltered top-level deltag

Energy values come from the top-level ``deltag``/``deltagerr`` gated by the chosen
level's eligibility -- the same source the cascade uses. Structure-only heuristics
(ATP synthase, ABC transporter, default) are evaluated even when a reaction has no
ΔG at the chosen level; the ΔG-dependent heuristics are recorded as ``NA`` there.

Outputs (written next to the reversibility reports):
  * ``Reaction_Heuristic_Outputs[_<LEVEL>].csv``  -- one row per reaction, one
    column per heuristic (the operator it returned, ``-`` if it abstained, ``NA``
    if it could not run), plus the cascade's winning rule + operator.
  * ``Reaction_Heuristic_Outputs[_<LEVEL>].json`` -- richer per-reaction detail
    including each heuristic's status label.
"""
import csv
import json
import sys

sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions

from reversibility_heuristics import (
    DEFAULT_HEURISTICS, evaluate_all_heuristics, top_level_energy,
)

# Column order = cascade order.
HEURISTIC_NAMES = [h.__name__ for h in DEFAULT_HEURISTICS]


def _parse_db_level(argv):
    if len(argv) > 1 and argv[1] in ('EQ', 'GC', 'DGP', ''):
        return argv[1]
    return 'EQ'


def _cell(result):
    """Flatten one heuristic result to a CSV cell.

    ``(status, op)`` -> the operator; ``None`` -> ``-`` (abstained);
    ``"no-energy"`` / ``"empty"`` -> ``NA``.
    """
    if isinstance(result, tuple):
        return result[1]
    if result is None:
        return '-'
    return 'NA'  # "no-energy" or "empty"


def _cascade_winner(results):
    """First heuristic (in cascade order) that fired -> (name, operator)."""
    for name in HEURISTIC_NAMES:
        r = results.get(name)
        if isinstance(r, tuple):
            return name, r[1]
    return None, None


def main():
    db_level = _parse_db_level(sys.argv)
    suffix = ('_' + db_level) if db_level else ''
    energy_source = top_level_energy(db_level)

    helper = Reactions()
    reactions_dict = helper.loadReactions()

    csv_path = "Reaction_Heuristic_Outputs%s.csv" % suffix
    json_path = "Reaction_Heuristic_Outputs%s.json" % suffix

    detail = {}
    n_energy = 0
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["rxn_id", "has_energy", "energy_source"]
            + HEURISTIC_NAMES
            + ["cascade_winner", "cascade_operator"])

        for rxn in sorted(reactions_dict.keys()):
            rxn_entry = reactions_dict[rxn]
            source_label, has_energy, results = evaluate_all_heuristics(
                rxn_entry, energy_source, DEFAULT_HEURISTICS)
            if has_energy:
                n_energy += 1
            winner, winner_op = _cascade_winner(results)

            writer.writerow(
                [rxn, int(bool(has_energy)), source_label or ""]
                + [_cell(results[name]) for name in HEURISTIC_NAMES]
                + [winner or "", winner_op or ""])

            detail[rxn] = {
                "has_energy": bool(has_energy),
                "energy_source": source_label,
                "cascade": {"winner": winner, "operator": winner_op},
                "heuristics": {
                    name: (
                        {"status": results[name][0], "operator": results[name][1]}
                        if isinstance(results[name], tuple)
                        else (None if results[name] is None else results[name])
                    )
                    for name in HEURISTIC_NAMES
                },
            }

    with open(json_path, "w") as fh:
        json.dump(detail, fh, separators=(",", ":"), sort_keys=True)

    print("Heuristics (cascade order): %s" % ", ".join(HEURISTIC_NAMES))
    print("Reactions: %d  (with %s energy: %d)"
          % (len(reactions_dict),
             (db_level + "-level") if db_level else "top-level", n_energy))
    print("Wrote %s and %s" % (csv_path, json_path))


if __name__ == "__main__":
    main()
