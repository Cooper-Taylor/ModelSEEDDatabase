#!/usr/bin/env python3
"""One-command regression test for the Thermodynamics reaction-direction pipeline.

What it does, in order:

1. Bootstraps a read-only baseline by extracting ``Biochemistry/reaction_*.json``
   from a reference git branch (default ``dev``) into ``Scripts/Tests/dev_baseline/``.
   The baseline is gitignored; it is only re-extracted when missing or when
   ``--refresh-baseline`` is passed.
2. Runs the Thermodynamics pipeline (the same six commands as
   ``Scripts/Thermodynamics/Rerun_Thermodynamics.sh``). This mutates
   ``Biochemistry/reaction_*.json`` and ``Biochemistry/compound_*.json``
   in place — exactly what those scripts are designed to do.
3. Compares the ``reversibility`` field of every reaction in the current
   working-tree files against the baseline files and prints a pass/fail
   summary.

Usage:

    ./test_reaction_direction.py                 # bootstrap if needed, run, compare
    ./test_reaction_direction.py --no-run        # skip the pipeline; just diff
    ./test_reaction_direction.py --refresh-baseline  # re-pull baseline from dev
    ./test_reaction_direction.py --baseline-ref origin/dev  # use a different ref

Exit code is 0 when every reaction's direction matches the baseline,
1 otherwise.
"""
import argparse
import json
import os
import subprocess
import sys
from glob import glob

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(THIS_DIR, '..', '..'))
BIOCHEM_DIR = os.path.join(REPO_ROOT, 'Biochemistry')
THERMO_DIR = os.path.join(REPO_ROOT, 'Scripts', 'Thermodynamics')
BASELINE_DIR = os.path.join(THIS_DIR, 'dev_baseline')

REACTION_GLOB = 'reaction_*.json'
COMPOUND_GLOB = 'compound_*.json'

# Sentinel used by the Thermodynamics pipeline to mean "no usable energy".
SENTINEL_DG = 10000000

# Column caps for the mismatch table — keep the name column from blowing the
# width. Notes and the thermo sublist are shown in full so the reader can see
# every code (GCC, EQP, ...) and the full method label.
NAME_MAX = 40

# Map from _method_used() return value to the thermodynamics key it came from.
METHOD_KEYS = {'EQ': 'eQuilibrator', 'GC': 'Group contribution'}

# Human-readable label for the method that drove the final reversibility,
# shown in the "Thermo sublist (method)" column.
METHOD_LABELS = {'EQ': 'equilibrium', 'GC': 'group contribution'}

# The six scripts that make up the Thermodynamics reaction-direction
# pipeline, in the order Rerun_Thermodynamics.sh runs them.
PIPELINE = [
    ['./Update_Compound_GroupContribution_Energies.py'],
    ['./Update_Reaction_GroupContribution_Energies.py'],
    ['./Estimate_Reaction_Reversibility.py', 'GC'],
    ['./Update_Compound_eQuilibrator_Energies.py'],
    ['./Update_Reaction_eQuilibrator_Energies.py'],
    ['./Estimate_Reaction_Reversibility.py', 'EQ'],
]


def bootstrap_baseline(ref):
    """Pull dev's reaction_*.json and compound_*.json files into ``dev_baseline/``.

    Uses ``git show <ref>:<path>`` so the working tree of the reference
    branch is never touched. Files in ``dev_baseline/`` are read-only
    inputs to the test — never written by the pipeline. Compound files
    are needed so the diff can show per-compound energies from the
    baseline alongside the newly generated values."""
    os.makedirs(BASELINE_DIR, exist_ok=True)

    listing = subprocess.run(
        ['git', 'ls-tree', '-r', '--name-only', ref, 'Biochemistry/'],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True)
    names = [n for n in listing.stdout.splitlines()
             if (n.startswith('Biochemistry/reaction_')
                 or n.startswith('Biochemistry/compound_'))
             and n.endswith('.json')]
    if not names:
        sys.exit('ERROR: no Biochemistry/reaction_*.json or compound_*.json '
                 'files found at ref ' + ref)

    print('Bootstrapping baseline from {0} ({1} files)...'.format(ref, len(names)))
    for name in names:
        blob = subprocess.run(
            ['git', 'show', '{0}:{1}'.format(ref, name)],
            cwd=REPO_ROOT, capture_output=True, check=True)
        out_path = os.path.join(BASELINE_DIR, os.path.basename(name))
        with open(out_path, 'wb') as fh:
            fh.write(blob.stdout)
    print('Baseline written to {0}'.format(BASELINE_DIR))


def baseline_present():
    return (bool(glob(os.path.join(BASELINE_DIR, REACTION_GLOB)))
            and bool(glob(os.path.join(BASELINE_DIR, COMPOUND_GLOB))))


def run_pipeline():
    """Run the six Thermodynamics scripts from the Thermodynamics directory.

    The scripts use ``sys.path.append('../../Libs/Python/')`` and other
    paths relative to their own directory, so cwd must be THERMO_DIR."""
    print('Running Thermodynamics pipeline in {0}...'.format(THERMO_DIR))
    for cmd in PIPELINE:
        print('  $ ' + ' '.join(cmd))
        subprocess.run(cmd, cwd=THERMO_DIR, check=True)


def load_reactions(directory):
    """Return ``{reaction_id: rxn_dict}`` for every reaction_*.json
    file in ``directory``."""
    out = {}
    paths = sorted(glob(os.path.join(directory, REACTION_GLOB)))
    if not paths:
        sys.exit('ERROR: no reaction_*.json files in ' + directory)
    for path in paths:
        with open(path) as fh:
            for rxn in json.load(fh):
                out[rxn['id']] = rxn
    return out


def load_compounds(directory):
    """Return ``{compound_id: cpd_dict}`` for every compound_*.json
    file in ``directory``. Returns an empty dict if no compound files
    exist (older baselines that predate compound bootstrapping)."""
    out = {}
    paths = sorted(glob(os.path.join(directory, COMPOUND_GLOB)))
    for path in paths:
        with open(path) as fh:
            for cpd in json.load(fh):
                out[cpd['id']] = cpd
    return out


def _method_used(rxn):
    """Which energy source drove the final reversibility for ``rxn``.

    Mirrors Estimate_Reaction_Reversibility's precedence: the EQ run runs
    last and overwrites the reversibility whenever the eQuilibrator sublist
    has a usable energy; otherwise the GC run's value survives."""
    thermo = rxn.get('thermodynamics') or {}
    eq = thermo.get('eQuilibrator')
    if eq and eq[0] is not None and eq[0] != SENTINEL_DG:
        return 'EQ'
    gc = thermo.get('Group contribution')
    if gc and gc[0] is not None and gc[0] != SENTINEL_DG:
        return 'GC'
    return '—'


def _format_thermo_sublist(rxn):
    """Render the thermodynamics sublist for whichever method drove the final
    reversibility (EQ wins over GC; falls back to GC if EQ has no usable
    energy). The method tag in parentheses tells the reader which key the
    sublist came from."""
    thermo = rxn.get('thermodynamics') or {}
    method = _method_used(rxn)
    sublist = thermo.get(METHOD_KEYS.get(method, ''))
    if not sublist:
        body = '—'
    else:
        body = '[' + ', '.join(_fmt_cell(x) for x in sublist) + ']'
    return '{0}  ({1})'.format(body, METHOD_LABELS.get(method, method))


def _format_notes(rxn):
    notes = rxn.get('notes')
    if not notes:
        return '—'
    return ', '.join(notes)


def _compound_energy(cpd, method):
    """Return the [deltag, deltagerr] energy for ``method`` from ``cpd``'s
    thermodynamics, or ``None`` if missing/unusable."""
    if not cpd:
        return None
    key = METHOD_KEYS.get(method)
    if not key:
        return None
    sublist = (cpd.get('thermodynamics') or {}).get(key)
    if not sublist or sublist[0] is None or sublist[0] == SENTINEL_DG:
        return None
    return sublist


def _format_compound_energies(rxn, base_cpds, cur_cpds):
    """Render per-compound energies for the method that drove ``rxn``'s
    final reversibility, in ``cpdXXXXX: base→cur`` form for each
    compound in the reaction's stoichiometry. A bare value means both
    baseline and current agree."""
    method = _method_used(rxn)
    if method not in METHOD_KEYS:
        return '—'
    stoich = rxn.get('stoichiometry') or []
    parts = []
    for entry in stoich:
        cid = entry.get('compound')
        if not cid:
            continue
        base_e = _compound_energy(base_cpds.get(cid), method)
        cur_e = _compound_energy(cur_cpds.get(cid), method)
        base_s = _fmt_deltag(base_e[0] if base_e else None)
        cur_s = _fmt_deltag(cur_e[0] if cur_e else None)
        if base_s == cur_s:
            parts.append('{0}: {1}'.format(cid, cur_s))
        else:
            parts.append('{0}: {1}→{2}'.format(cid, base_s, cur_s))
    if not parts:
        return '—'
    return ', '.join(parts)


def _fmt_cell(x):
    if isinstance(x, str):
        return "'" + x + "'"
    if isinstance(x, float):
        return '{0:g}'.format(x)
    return str(x)


def _fmt_rev(rev):
    return rev if rev else '—'


def _fmt_deltag(dg):
    if dg is None or dg == SENTINEL_DG:
        return '—'
    return '{0:g}'.format(dg)


def _truncate(s, width):
    if len(s) <= width:
        return s
    return s[:width - 1] + '…'


def _print_row(cells, widths):
    print('| ' + ' | '.join(c.ljust(w) for c, w in zip(cells, widths)) + ' |')


def _print_mismatch_table(mismatched, current, baseline,
                          cur_cpds, base_cpds, max_show):
    headers = ['Reaction', 'Name', 'dev → new', 'Notes',
               'Thermo sublist (method)', 'Top-level `deltag`',
               'Compound energies (dev → new)']
    rows = []
    for rid in mismatched[:max_show]:
        cur = current[rid]
        base = baseline[rid]
        rows.append([
            rid,
            _truncate(cur.get('name') or '—', NAME_MAX),
            '{0} → {1}'.format(_fmt_rev(base.get('reversibility')),
                               _fmt_rev(cur.get('reversibility'))),
            _format_notes(cur),
            _format_thermo_sublist(cur),
            _fmt_deltag(cur.get('deltag')),
            _format_compound_energies(cur, base_cpds, cur_cpds),
        ])

    widths = [max(len(h), max((len(r[i]) for r in rows), default=0))
              for i, h in enumerate(headers)]

    print()
    if len(mismatched) > max_show:
        print('Direction mismatches (first {0} of {1}):'.format(
            max_show, len(mismatched)))
    else:
        print('Direction mismatches ({0}):'.format(len(mismatched)))
    _print_row(headers, widths)
    print('|' + '|'.join('-' * (w + 2) for w in widths) + '|')
    for r in rows:
        _print_row(r, widths)


def compare(current, baseline, cur_cpds, base_cpds, max_show=20):
    """Print a summary diff and return True iff every reaction matches."""
    cur_ids = set(current)
    base_ids = set(baseline)

    only_current = cur_ids - base_ids
    only_baseline = base_ids - cur_ids
    shared = cur_ids & base_ids
    mismatched = [rid for rid in sorted(shared)
                  if baseline[rid].get('reversibility')
                  != current[rid].get('reversibility')]

    print()
    print('=' * 60)
    print('Reaction direction comparison')
    print('=' * 60)
    print('  Reactions in baseline : {0}'.format(len(baseline)))
    print('  Reactions in current  : {0}'.format(len(current)))
    print('  Shared reactions      : {0}'.format(len(shared)))
    print('  Only in current       : {0}'.format(len(only_current)))
    print('  Only in baseline      : {0}'.format(len(only_baseline)))
    print('  Direction mismatches  : {0}'.format(len(mismatched)))

    if only_current:
        print()
        print('Reactions only in current (first {0}):'.format(max_show))
        for rid in sorted(only_current)[:max_show]:
            print('  + {0}'.format(rid))
    if only_baseline:
        print()
        print('Reactions only in baseline (first {0}):'.format(max_show))
        for rid in sorted(only_baseline)[:max_show]:
            print('  - {0}'.format(rid))
    if mismatched:
        _print_mismatch_table(mismatched, current, baseline,
                              cur_cpds, base_cpds, max_show)

    ok = not (only_current or only_baseline or mismatched)
    print()
    print('RESULT: ' + ('PASS' if ok else 'FAIL'))
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--baseline-ref', default='dev',
                        help='git ref to pull the baseline from (default: dev)')
    parser.add_argument('--refresh-baseline', action='store_true',
                        help='re-extract the baseline from --baseline-ref')
    parser.add_argument('--no-run', action='store_true',
                        help='skip the Thermodynamics pipeline; just diff '
                             'the current Biochemistry/ against the baseline')
    parser.add_argument('--display-num', type=int, default=20, metavar='N',
                        help='max rows to show in each diff section '
                             '(default: 20)')
    args = parser.parse_args()

    if args.display_num < 0:
        parser.error('--display-num must be non-negative')

    if args.refresh_baseline or not baseline_present():
        bootstrap_baseline(args.baseline_ref)
    else:
        print('Baseline already present at {0} (use --refresh-baseline to '
              'rebuild)'.format(BASELINE_DIR))

    if not args.no_run:
        run_pipeline()
    else:
        print('Skipping pipeline run (--no-run)')

    baseline = load_reactions(BASELINE_DIR)
    current = load_reactions(BIOCHEM_DIR)
    base_cpds = load_compounds(BASELINE_DIR)
    cur_cpds = load_compounds(BIOCHEM_DIR)
    ok = compare(current, baseline, cur_cpds, base_cpds,
                 max_show=args.display_num)
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
