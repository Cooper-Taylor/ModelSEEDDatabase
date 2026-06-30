# Reversibility heuristics: analysis and ideal order

The reaction-direction cascade (`reversibility_heuristics.DEFAULT_HEURISTICS`,
consumed by `Estimate_Reaction_Reversibility.py`) is an ordered list of
heuristics. Each is a function `(ctx) -> (status, operator) | None`; the cascade
runs them in order and the **first non-`None` result wins**. Order therefore
determines which rule gets to speak for a reaction that more than one rule could
classify. This note analyzes each heuristic and justifies the order.

## The six heuristics, by specificity

| # | heuristic | inputs | verdict | specificity |
|---|-----------|--------|---------|-------------|
| 1 | `atp_synthase_heuristic` | structure only (transport flag, proton compartments, exact 5-reagent ATPS signature) | `=` | **exact structural match** |
| 2 | `abc_transporter_heuristic` | structure only (transport flag + ATP among reagents) | sign of the ATP coefficient (`>`/`<`) | structural pattern |
| 3 | `stored_bounds_heuristic` | ΔG′, σ, concentration window | `>` if max ΔG over the window < 0; `<` if min > 0 | general thermodynamic (range) |
| 4 | `mmdeltag_band_heuristic` | ΔG′ at 1 mM | `=` if \|mMΔG\| ≤ 2 kcal/mol | general thermodynamic (point) |
| 5 | `low_energy_heuristic` | ΔG′ at 1 mM + low-energy/phosphate "points" score | directional if points·mMΔG > 2 | weak heuristic |
| 6 | `default_heuristic` | none | `=` (reversible) | universal fallback |

**Energy-free vs energy-dependent.** Heuristics 1, 2, 6 need only the
stoichiometry/structure; 3, 4, 5 need a ΔG estimate. So 1/2 can (and should) be
decided for reactions whose direction is structurally determined regardless of
whether a ΔG number exists.

## The ideal order: most-specific → least-specific

```
atp_synthase  →  abc_transporter  →  stored_bounds  →  mmdeltag_band  →  low_energy  →  default
```

The governing principle: **a rule that encodes specific knowledge about *this kind
of reaction* must out-rank a rule that only knows its ΔG number.** Specificity, not
"thermodynamics first," is the right sort key.

1. **`atp_synthase` first.** ATP synthase is physiologically reversible — its net
   direction is set by the proton-motive force, not by ΔG′° of the chemistry. The
   reaction *also* has a ΔG profile that, read naively, looks committed. So the
   exact structural match must fire before any ΔG rule, or the reaction is
   mislabeled directional. This is the single most specific rule (it matches only
   reactions with the exact 5-reagent, multi-compartment proton signature).
2. **`abc_transporter` second.** An ATP-driven transporter's direction follows the
   ATP coupling (hydrolysis drives uptake forward), not the transported
   metabolite's small ΔG. Structural, so it precedes the general ΔG rules; less
   specific than ATPS (any transport reaction with ATP qualifies), so it sits
   second.
3. **`stored_bounds` third.** The strongest *thermodynamic* test: commit a
   direction only if the ΔG range across the whole physiological concentration
   window (1e-5–0.02 M) stays one-signed. Being range-based it is more
   conservative than the point estimate below it, so it runs first among the ΔG
   rules — but it is a general rule with no knowledge of reaction role, so it must
   yield to the two structural rules above.
4. **`mmdeltag_band` fourth.** If the point ΔG at 1 mM is within ±2 kcal/mol the
   reaction is near equilibrium → `=`. This is the "too small to commit" fallback
   once the stricter range test (3) declined to commit.
5. **`low_energy` fifth.** A weaker, score-based directional call (low-energy
   compounds + phosphate spread) for cases the band left undecided. It is the most
   heuristic of the ΔG rules, so it runs last before the default.
6. **`default` last.** Nothing committed → call it reversible (`=`). The universal
   terminal rule, necessarily least specific.

## Why this matters (measured)

The pre-fix order ran `stored_bounds` **first**, so structurally-determined
reactions were captured by the general ΔG rule before their structural rule could
fire. Moving `atp_synthase` + `abc_transporter` ahead of `stored_bounds` changed
**44 of 56,012 reactions** (14 ATP synthases `<`/`>`→`=`; 30 ATP-driven
transporters now following the ATP-coefficient sign) with **no change to any ΔG
value**. In the KEGG core-model panel the only model-relevant case, F(1)-ATPase
(`rxn08173`, `<`→`=`), lifts biomass flux in 91/100 models with zero grow-flips —
a more realistic flux solution, same growth calls.

## Diagnostic: every heuristic on every reaction

Because the cascade short-circuits at the first match, the stored report only shows
the *winning* rule. To see what **every** heuristic says for **every** reaction
(independent of order), use `evaluate_all_heuristics()` in
`reversibility_heuristics.py` and the driver `Dump_Reaction_Heuristic_Outputs.py`,
which writes a per-reaction × per-heuristic matrix to
`Reaction_Heuristic_Outputs.{csv,json}`. That matrix is what justifies/ء audits the
order: it shows, for each reaction, which rules agree, which disagree, and which
rule the cascade ultimately picked.
