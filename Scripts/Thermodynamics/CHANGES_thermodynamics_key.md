# Thermodynamics-key migration for `Estimate_Reaction_Reversibility.py`

Migrates the reversibility-estimation step from the legacy per-source notes
(`GCC`/`EQU`/...) to the structured `thermodynamics` block on each reaction,
and stamps the computed reversibility back into the `Group contribution` /
`eQuilibrator` sublists so each entry self-describes.

## Files touched

- `Scripts/Thermodynamics/Estimate_Reaction_Reversibility.py` — refactor
  described below. Two new private helpers (`_thermo_pair`, `_energy_for`,
  `_has_gc_data`, `_append_direction`) plus a `DB_LEVEL_LABEL` table that
  encodes the `GC` / `EQ` → sublist-key mapping.
- `Biochemistry/reaction_60.json` — added `thermodynamics: {"Group
  contribution": [10000000.0, 10000000.0]}` to the three reactions that
  did not carry the key (`rxn60857`, `rxn60858`, `rxn60859`). Format
  mirrors the existing no-data reactions such as `rxn00014`.

## Refactor details

### Eligibility (replaces the notes check)

The pre-refactor script consulted `rxn_entry["notes"]` to decide whether a
GC or EQ run should attempt an estimate:

```python
# old
for note in rxn_entry["notes"]:
    if db_level in note and note in ("GCC", "EQU"):
        return True
```

The new path consults `thermodynamics[label]` directly, where `label`
comes from a single mapping table:

```python
DB_LEVEL_LABEL = {"GC": "Group contribution", "EQ": "eQuilibrator"}

def _thermo_pair(rxn_entry, label):
    thermo = rxn_entry.get('thermodynamics')
    if not isinstance(thermo, dict):
        return None
    pair = thermo.get(label)
    if not pair or pair[0] is None:
        return None
    dg = float(pair[0])
    if dg == SENTINEL_DG:
        return None
    return [dg, float(pair[1])]
```

Adding a third source means appending one entry to `DB_LEVEL_LABEL` and
one entry to `DB_LEVEL_PRIORITY`.

### Energy values (unchanged for output equivalence)

The dg/dge values fed into the heuristic cascade still come from the
top-level `deltag`/`deltagerr` fields. This keeps the numeric content of
the existing `Estimated_Reaction_Reversibility_Report*.txt` reports
byte-identical to the dev-branch baseline for every reaction whose
eligibility outcome did not change.

### EQ fallback (`Incomplete (GCC)`)

The pre-refactor fallback re-used the `"GCC" in notes` check to decide
whether an EQ run could keep the reversibility a prior GC run wrote. The
new check uses `_has_gc_data()`, i.e. "the reaction has non-sentinel
`thermodynamics['Group contribution']`". This is the same intent expressed
against the new source of truth.

### Appending the direction

After estimation, the chosen sublist is rewritten as `[dg, dge, direction]`
so the Thermodynamics block records both the energy and the resulting
direction for each source. `_append_direction` is idempotent — a second
run truncates any prior direction before appending.

```
"thermodynamics": {
    "Group contribution": [4.15, 1.22, ">"],
    "eQuilibrator":       [-3.46, 0.05, ">"]
}
```

`GC` runs only touch the `Group contribution` sublist; `EQ` runs only
touch `eQuilibrator`. The unfiltered run touches whichever sublist's
first entry matches the top-level `deltag` (with `eQuilibrator` winning
ties, mirroring the "EQ overwrites GC" precedence used by
`Update_*_eQuilibrator_Energies.py`). When no sublist matches the
top-level value, nothing is appended — leaving the legacy data
unchanged for those edge cases.

## Output verification vs `dev`

Compared against the report files committed at `dev` (commit
`33d5d84`) and re-running the Estimate scripts only (i.e. without
`Update_*_Energies.py`, which is independently non-idempotent for 279
reactions):

| Mode    | Report-line diff | Reactions w/ different final reversibility |
|---------|------------------|-------------------------------------------|
| no-arg  |  86              |  0                                        |
| `GC`    | 174              | 48                                        |
| `EQ`    | 13608            | 24                                        |

After running the canonical sequence
`Estimate_Reaction_Reversibility.py GC && Estimate_Reaction_Reversibility.py EQ && Estimate_Reaction_Reversibility.py`,
the stored `reversibility` field in the JSON matches the dev-branch
baseline for **all 56,012** reactions.

The bulk of the report-line diff is the status column (e.g.
`MdeltaG(Min): 114.95` → `Incomplete (GCC)`) for reactions whose
eligibility outcome changed even though the final reversibility did
not. The 48-and-24 stored-reversibility deltas in the intermediate
GC/EQ reports are exclusively reactions where the legacy notes
claimed source-completeness that the actual Thermodynamics data did
not back up; those reactions correctly resolve to `?` after the EQ
run and recover to dev's value after the unfiltered no-arg run.

## What is *not* changed

- `Libs/Python/BiochemPy/Reactions.py` — untouched.
- `_thermo_helpers.py` — untouched. `set_thermo` still writes whatever
  `[dg, dge]` value the Update_*_Energies scripts hand it; the
  append-direction step lives only inside the Estimate script's main
  loop, so a subsequent Update run will simply overwrite (and the next
  Estimate run will re-append).
- `Rerun_Thermodynamics.sh` — untouched. Existing operators get the
  same flow.
- Notes (`GCC`/`EQU`/`GCP`/`EQP`/etc.) — left in place for backward
  compatibility with other readers. The Estimate script still consults
  the per-source completeness flags (`GCC` and `EQU`) as a fallback —
  see the "Notes-as-fallback fix" section below.

## Notes-as-fallback fix

The original migration claim that "the stored `reversibility` field in
the JSON matches the dev-branch baseline for all 56,012 reactions" was
incorrect for 54 reactions. Those reactions carry the legacy `GCC`
completeness flag in `notes` but their structured
`thermodynamics['Group contribution']` sublist holds the sentinel
`[10000000.0, 10000000.0]` — the Update_Reaction_GroupContribution
script writes the sentinel for reactions with incomplete per-compound
GC coverage, and the legacy `GCC` flag was already set on dev for some
of those.

dev's notes-based eligibility check accepted these reactions and
produced a real direction (`=`, `>`, or `<`); the migrated
thermo-sublist-only check rejected them and assigned `?`. The
`Scripts/Tests/test_reaction_direction.py` regression test caught all
54 mismatches.

### Fix

Eligibility is now the OR of the two completeness signals:

```python
DB_LEVEL_NOTE = {"GC": "GCC", "EQ": "EQU"}

def _is_source_eligible(rxn_entry, level):
    if _thermo_pair(rxn_entry, DB_LEVEL_LABEL[level]) is not None:
        return True
    return DB_LEVEL_NOTE[level] in rxn_entry["notes"]
```

Adding a new source still requires only adding entries to
`DB_LEVEL_LABEL` (the structured sublist key) and `DB_LEVEL_NOTE`
(the legacy notes flag) — and optionally to `DB_LEVEL_PRIORITY` for
the unfiltered-run tie-break. The modular-extension story from the
original migration is preserved; the OR just makes the legacy notes
load-bearing again for reactions whose source-completeness predates
the structured Thermodynamics block.

When eligibility comes from the notes alone (structured sublist
absent or sentinel), the post-estimate direction-append target is set
to `None` so the structured sublist is left untouched — there's no
`[dg, dge]` pair to extend into `[dg, dge, direction]`.

### Verification

After the fix, `Scripts/Tests/test_reaction_direction.py` reports
**0 direction mismatches across 56,012 reactions** against the dev
baseline (commit `33d5d84`).
