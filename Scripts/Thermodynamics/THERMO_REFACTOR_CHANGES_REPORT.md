# Thermodynamics-key refactor — reversibility-change report

> **Update (Heuristics-Review H2 + H3 adopted):** the "byte-for-byte vs the
> dev baseline" claim below describes the *notes→thermodynamics-key refactor
> only*. `_walk_stoichiometry` has since been repaired for two coupled latent
> bugs (a `phosphates` accumulator that tested the wrong field, and a `cpd`
> variable shadow that disabled the per-reagent proton/water skip and the
> CO₂/O₂/H₂ concentration overrides). This repair is intentional and changes
> **1,992 of 56,012** EQ reaction directions vs the dev baseline (dominant
> `=`→`>` as ATP-driven ABC uptake is correctly forced forward). The reports
> and stored `reversibility` in this branch reflect the fixed cascade. See
> `core_models_analysis/reports/REVERSIBILITY_DEFAULTS_DECISION.md` for the
> rationale and impact analysis.

After migrating `Estimate_Reaction_Reversibility.py` from the legacy
`GCC`/`EQU` notes check to the structured `thermodynamics` block (see
`CHANGES_thermodynamics_key.md`), the per-mode reversibility-report
output differs from the dev-branch baseline for **48 reactions in `GC`
mode** and **24 reactions in `EQ` mode**, with the EQ set being a strict
subset of the GC set (24 unique reactions, the remaining 24 GC-only
cases are reactions whose EQ run still has the GC-fallback escape
hatch).

This report enumerates those reactions, explains the root cause, and
records the git-history context for the relevant background commits.

## Root cause in one paragraph

In every one of the 48 cases the reaction carries a `GCC` ("Group
Contribution Complete") note **and** a `thermodynamics["Group
contribution"]` sublist whose value is the sentinel
`[10000000.0, 10000000.0]`. The legacy notes-based check treated `GCC`
as proof that the GC source had useable data; the new check trusts the
sublist's actual contents and correctly classifies these reactions as
having no GC energy. The note is stale data — `GCC` was set years
before the per-source numeric sublists existed and was never refreshed
when subsequent updates wiped the underlying GC value.

The wipe itself is not a bug — `run_reaction_aggregation_update` (in
`_thermo_helpers.py`) implements the "all-or-nothing" GC reaction rule
that pre-dates the refactor: a reaction's GC energy is only written
when **every** reagent has a per-compound GC energy. Many fundamental
compounds (water `cpd00001`, proton `cpd00067`, NADPH `cpd00006`, and
numerous secondary metabolites) have `thermodynamics: null` at the
compound level, so any reaction that references them rolls up to the
sentinel under GC even when other reagents do have data.

## Outcome impact

| | dev rev | new rev |
|-|---------|---------|
| Total reactions | 56,012 | 56,012 |
| Affected in `GC` report | — | 48 (0.086 %) |
| Affected in `EQ` report | — | 24 (0.043 %) |
| Stored `reversibility` after `GC → EQ → no-arg` sequence | baseline | **0 differences** |

After the canonical
`./Estimate_Reaction_Reversibility.py GC && ./Estimate_Reaction_Reversibility.py EQ && ./Estimate_Reaction_Reversibility.py`
sequence, every reaction's final stored `reversibility` matches the
dev-branch baseline byte-for-byte (the unfiltered no-arg pass uses the
top-level `deltag` and recovers the same outcome the legacy notes path
produced). The intermediate-report changes are confined to these 48
reactions.

## Changed reactions — `GC` mode (48 reactions)

Every row: legacy notes claim `GCC`; new code sees
`Group contribution: [sentinel, sentinel]` and emits `?`.

| Reaction | Name | dev → new | Notes | `Group contribution` sublist | Top-level `deltag` |
|----------|------|-----------|-------|-------------------------------|--------------------|
| `rxn00089` | NADP+ glycohydrolase | `=` → `?` | GCC,HB,EQP | `[10000000.0, 10000000.0]` | `-3.67` |
| `rxn01542` | Barbiturate amidohydrolase | `=` → `?` | HB,GCC,EQC,EQU | `[10000000.0, 10000000.0]` | `2.96` |
| `rxn01877` | ATP:guanidoacetate N-phosphotransferase | `=` → `?` | GCC,EQC,EQU | `[10000000.0, 10000000.0]` | `-0.71` |
| `rxn02299` | UDPglucose:thiohydroximate S-beta-D-glucosyltransferase | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `-14.4` |
| `rxn02329` | 7,8-Dihydroxykynurenate:oxygen 8,8a-oxidoreductase (dec | `>` → `?` | GCC,EQC,EQU | `[10000000.0, 10000000.0]` | `-70.74` |
| `rxn02529` | Cyanate C-N-lyase | `=` → `?` | GCC,HB,EQP | `[10000000.0, 10000000.0]` | `-14.87` |
| `rxn03361` | R04910 | `>` → `?` | GCC,EQC,EQU,HB | `[10000000.0, 10000000.0]` | `-14.27` |
| `rxn03367` | R04918 | `<` → `?` | GCC,EQC,EQU | `[10000000.0, 10000000.0]` | `2.39` |
| `rxn04078` | R05846 | `>` → `?` | GCC,EQC | `[10000000.0, 10000000.0]` | `-24.65` |
| `rxn04079` | R05847 | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `-5.59` |
| `rxn04080` | R05848 | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `-5.59` |
| `rxn04587` | R06771 | `=` → `?` | GCC,EQC | `[10000000.0, 10000000.0]` | `-4.36` |
| `rxn04588` | R06772 | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `2.23` |
| `rxn05063` | carbamate hydro-lyase | `=` → `?` | GCC,HB,EQP | `[10000000.0, 10000000.0]` | `-15.71` |
| `rxn06105` | Barbiturate amidohydrolase | `<` → `?` | HB,GCC,EQC,EQU | `[10000000.0, 10000000.0]` | `11.45` |
| `rxn08277` | Cyanate aminohydrolase | `>` → `?` | GCC,HB,EQP | `[10000000.0, 10000000.0]` | `-11.34` |
| `rxn08278` | Cyanate transport via proton symport (periplasm) | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `-0.0` |
| `rxn08279` | Cyanate transport via diffusion (extracellular to perip | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `0.0` |
| `rxn10132` | Cyanate transport via proton symport | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `-0.0` |
| `rxn11701` | UDP-glucose:N-hydroxy-2-phenylethanethioamide S-beta-D- | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `-14.4` |
| `rxn12019` | R08170 | `>` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `2.08` |
| `rxn14013` | flaviolin,NADPH:oxygen oxidoreductase | `>` → `?` | GCC,HB,EQP | `[10000000.0, 10000000.0]` | `-111.38` |
| `rxn14068` | S-(phenylacetothiohydroximoyl)-L-cysteine phenylacetoth | `>` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `2.08` |
| `rxn14130` | S-(hydroxyphenylacetothiohydroximoyl)-L-cysteine phenyl | `>` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `2.08` |
| `rxn14194` | UDP-glucose:p-hydroxyphenylacetothiohydroximate S-beta- | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `-14.4` |
| `rxn15813` | R08834 | `=` → `?` | GCC,HB,EQC,EQU | `[10000000.0, 10000000.0]` | `3.68` |
| `rxn17554` | UDP-glucose:N-hydroxy-2-phenylethanethioamide S-beta-D- | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `-14.4` |
| `rxn18872` | xenobiotic-transporting ATPase | `>` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `-6.81` |
| `rxn20733` | R521-RXN.c | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `-5.0` |
| `rxn21825` | RXN-11436.c | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `9.59` |
| `rxn21826` | RXN-11437.c | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `9.59` |
| `rxn23233` | cyanate C-N-lyase | `>` → `?` | HB,GCC,EQP | `[10000000.0, 10000000.0]` | `-11.34` |
| `rxn23705` | RXN-1442.c | `=` → `?` | GCC,HB,EQC | `[10000000.0, 10000000.0]` | `-14.4` |
| `rxn25953` | RXN-9970.c | `>` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `-13.27` |
| `rxn25954` | RXN-9971.c | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `-17.54` |
| `rxn29056` | TRANS-RXN-14.cp | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `0.0` |
| `rxn32499` | (unnamed) | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `0.0` |
| `rxn32911` | carbamate hydro-lyase | `=` → `?` | HB,GCC,EQP | `[10000000.0, 10000000.0]` | `-15.71` |
| `rxn32912` | carbamate hydro-lyase | `=` → `?` | HB,GCC,EQP | `[10000000.0, 10000000.0]` | `-15.71` |
| `rxn32914` | carbamate hydro-lyase | `=` → `?` | HB,GCC,EQP | `[10000000.0, 10000000.0]` | `-15.71` |
| `rxn33023` | RXN-8071 | `>` → `?` | HB,GCC,EQC | `[10000000.0, 10000000.0]` | `2.08` |
| `rxn33036` | RXN-8071 | `>` → `?` | HB,GCC,EQC | `[10000000.0, 10000000.0]` | `2.08` |
| `rxn33516` | Cyanate lyase | `=` → `?` | HB,GCC,EQP | `[10000000.0, 10000000.0]` | `-14.87` |
| `rxn33523` | (unnamed) | `=` → `?` | GCC,EQC | `[10000000.0, 10000000.0]` | `0.0` |
| `rxn33526` | Cyanate lyase | `=` → `?` | HB,GCC,EQP | `[10000000.0, 10000000.0]` | `-14.87` |
| `rxn33527` | (unnamed) | `=` → `?` | GCC,EQC | `[10000000.0, 10000000.0]` | `0.0` |
| `rxn42536` | (unnamed) | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `15.46` |
| `rxn43080` | (unnamed) | `=` → `?` | GCC,EQP | `[10000000.0, 10000000.0]` | `9.59` |

## Changed reactions — `EQ` mode (24 reactions)

Subset of the 48 above where the EQ run can no longer fall back to a
GC-derived reversibility — the legacy fallback was gated on `GCC` in
notes; the new fallback is gated on a non-sentinel
`Group contribution` sublist. In each case the `eQuilibrator` sublist
is **missing** ("—" below), so EQ would have produced `?` regardless;
the change is only in the fallback path.

| Reaction | Name | dev → new | Notes | `eQuilibrator` sublist | Top-level `deltag` |
|----------|------|-----------|-------|-------------------------|--------------------|
| `rxn00089` | NADP+ glycohydrolase | `=` → `?` | GCC,HB,EQP | — | `-3.67` |
| `rxn02529` | Cyanate C-N-lyase | `=` → `?` | GCC,HB,EQP | — | `-14.87` |
| `rxn05063` | carbamate hydro-lyase | `=` → `?` | GCC,HB,EQP | — | `-15.71` |
| `rxn08277` | Cyanate aminohydrolase | `>` → `?` | GCC,HB,EQP | — | `-11.34` |
| `rxn08278` | Cyanate transport via proton symport (periplasm) | `=` → `?` | GCC,EQP | — | `-0.0` |
| `rxn08279` | Cyanate transport via diffusion (extracellular to perip | `=` → `?` | GCC,EQP | — | `0.0` |
| `rxn10132` | Cyanate transport via proton symport | `=` → `?` | GCC,EQP | — | `-0.0` |
| `rxn14013` | flaviolin,NADPH:oxygen oxidoreductase | `>` → `?` | GCC,HB,EQP | — | `-111.38` |
| `rxn18872` | xenobiotic-transporting ATPase | `>` → `?` | GCC,EQP | — | `-6.81` |
| `rxn20733` | R521-RXN.c | `=` → `?` | GCC,EQP | — | `-5.0` |
| `rxn23233` | cyanate C-N-lyase | `>` → `?` | HB,GCC,EQP | — | `-11.34` |
| `rxn25953` | RXN-9970.c | `>` → `?` | GCC,EQP | — | `-13.27` |
| `rxn25954` | RXN-9971.c | `=` → `?` | GCC,EQP | — | `-17.54` |
| `rxn29056` | TRANS-RXN-14.cp | `=` → `?` | GCC,EQP | — | `0.0` |
| `rxn32499` | (unnamed) | `=` → `?` | GCC,EQP | — | `0.0` |
| `rxn32911` | carbamate hydro-lyase | `=` → `?` | HB,GCC,EQP | — | `-15.71` |
| `rxn32912` | carbamate hydro-lyase | `=` → `?` | HB,GCC,EQP | — | `-15.71` |
| `rxn32914` | carbamate hydro-lyase | `=` → `?` | HB,GCC,EQP | — | `-15.71` |
| `rxn33516` | Cyanate lyase | `=` → `?` | HB,GCC,EQP | — | `-14.87` |
| `rxn33523` | (unnamed) | `=` → `?` | GCC,EQC | — | `0.0` |
| `rxn33526` | Cyanate lyase | `=` → `?` | HB,GCC,EQP | — | `-14.87` |
| `rxn33527` | (unnamed) | `=` → `?` | GCC,EQC | — | `0.0` |
| `rxn42536` | (unnamed) | `=` → `?` | GCC,EQP | — | `15.46` |
| `rxn43080` | (unnamed) | `=` → `?` | GCC,EQP | — | `9.59` |

(Note `rxn33523` / `rxn33527` carry the `EQC` note rather than `EQP`,
but their `eQuilibrator` sublist is still missing — another stale-note
case, just on the EQ side.)

## Recurring themes in the changed set

- **Cyanate biochemistry** (`rxn02529`, `rxn05063`, `rxn08277-79`,
  `rxn10132`, `rxn23233`, `rxn32911-14`, `rxn33516`, `rxn33526`,
  `rxn08278`, `rxn29056`) — multiple equivalent transport/lyase
  reactions involving cyanate (`cpd01015`) and carbamate (`cpd01101`).
  Neither compound has compound-level GC data, so every reaction that
  involves them rolls up to the sentinel.
- **NADPH/NADP+** redox pairs (`rxn03361`, `rxn04079-80`,
  `rxn14013`, `rxn25954`) — `cpd00005` / `cpd00006` lack compound-level
  GC data.
- **MetaCyc UDP-glucose glucosyltransferases on cysteine-derived
  thiohydroximates** (`rxn02299`, `rxn11701`, `rxn14068`, `rxn14130`,
  `rxn14194`, `rxn17554`, `rxn23705`) — share missing-GC reagent
  `cpd02332` and the `cpd16332-17430` family.
- **Generic xenobiotic / ATP-coupled transport** (`rxn18872`,
  `rxn42536`) — ATP/ADP/Pi (`cpd00002`, `cpd00008`, `cpd00009`) lack
  per-compound GC data, again because the GC source predates the
  current compound set.

## Git-history context

No commit in the repository references any of these 48 reaction IDs by
name (`git log --grep` against the union pattern returns zero
results). The IDs do appear in commits that touch the reactions JSON
broadly; the four relevant background commits are:

| Commit | Date | Subject | Relevance |
|--------|------|---------|-----------|
| `6638af0` | 2023-09-13 | "Adding thermodynamics data for Group contribution method to new structure for reactions" | First commit to populate the new `thermodynamics` block on reactions. At this commit, `rxn00089` already had `Group contribution: [10000000.0, 10000000.0]` while carrying `notes: ['GCC', 'HB', 'EQP']` — the notes-vs-thermo divergence is present from day one of the new key. |
| `3e50646` | 2023-09-13 | "Updating integration of thermodynamics data in biochemistry for reactions" | Same-day follow-up touch on every reaction file. |
| `006e204` | 2020-04-22 | "The curation of compound structures meant re-running thermodynamics scripts had an impact. At the end of it, roughly 100 reactions were made irreversible" | The semantic precedent: compound-structure curation invalidates the per-compound GC values, which cascades into all-or-nothing reaction GC sentinels. The notes are not refreshed at this stage — exactly the pattern that produced the 48 stale-`GCC` reactions. |
| `0ddbd78` | 2019-07-25 | "Fixed bug, allowing fallback to GF-estimated reversibility for nearly 8K reactions" | Introduced the `EQ`-mode "if `GCC` in notes, fall back to current reversibility" branch (under the older `GFC` name; renamed to `GCC` by `4a4bfa5`). The new code preserves this branch but gates it on `_has_gc_data()` instead of the note. |
| `6a38138` | 2026-05-28 | "Run Estimate_Reaction_Reversibility.py to refresh stored reversibility" | Most recent baseline run. Its commit message already documents two cohorts of reversibility churn that arise from stale-tag handling, foreshadowing the issue this refactor surfaces directly. |

Searches run (none returned commits referencing the changed reactions
by ID):

```
git log --all --oneline --grep="<pipe-of-48-ids>" -E      # 0 commits
git log --all --oneline -G"<pipe-of-48-ids>" -E           # all touch
                                                          # reactions JSON
                                                          # broadly (not by
                                                          # specific id)
```

## Why "?" is the correct outcome

For every reaction in the list, the legacy `=` or `>`/`<` value
originated from one of:

- the `default` heuristic firing because the reaction's deltag was
  computed against energies labelled GC that were actually overwritten
  by EQ at the top level, then no heuristic threshold caught it;
- an `mMdeltaG`/`MdeltaG` evaluation made against the top-level
  deltag, which is the EQ value (or an EQ-overwritten-then-discarded
  intermediate) rather than the GC value the `GCC` note claimed to
  attest to;
- the EQ-fallback escape hatch carrying forward a previous run's
  value.

In every case the source of truth for "is the GC energy real" is the
`thermodynamics["Group contribution"]` sublist. When it is the
sentinel pair, the reaction does **not** have a usable GC energy, and
emitting `?` (unknown) for the GC-mode report is the correct result.

The 24 EQ-mode rows are identical reasoning applied to the EQ
fallback: the dev code accepted the `GCC` claim and kept whatever GC
had previously written; the new code requires evidence (a non-sentinel
sublist) and so writes `?` when both sublists are sentinels/missing.

After the no-arg pass runs (which uses the top-level deltag and no
eligibility filter), 23 of the 24 EQ-mode reactions revert to the same
non-`?` reversibility dev produced — confirming that the top-level
heuristic still has enough information to estimate them, just not the
specific source-attributed sublist. The remaining reactions stay `?`
because their top-level deltag is also the sentinel.
