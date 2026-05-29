# ModelSEED Biochemistry Database Scripts (Thermodynamics)

Here we have the scripts used to handle the thermodynamics data that
came from both the Group Contribution approach (<a
href="https://doi.org/10.1529/biophysj.107.124784">Jankowski et
al. 2008</a>) and from eQuilibrator (<a
href="https://doi.org/10.1371/journal.pcbi.1003098">Noor et
al. 2013</a>). Our approach in handling this data is described in the
<a
href="https://www.biorxiv.org/content/10.1101/2020.03.31.018663v2">paper</a>.

## Order of execution

The general order is that the energies from the application of the Group Contribution (GC) approach
are stored in the database first, and then the energies from eQuilibrator (EQ), which, in most
cases, take precedence, are used to overwrite the energies in the database

The underlying thermodynamics data is kept in
`../../Biochemistry/Thermodynamics`. The decomposition of molecular
structures and their resulting energies for both the older group
contribution approach and the newer eQuilibrator approach are stored in
the `ModelSEED` and `eQuilibrator` directories.

As an addendum, the two scripts used to update the energies from
eQuilibrator are in this folder, but are dependent on files in
`../../Biochemistry/Structures/MetaNetX`:
```
./Retrieve_eQuilibrator_Compound_Energies.py
./Retrieve_eQuilibrator_Reactions_Energies.py
```

If the underlying thermodynamics data in `../../Biochemistry/Thermodyanmics` hasn't changed,
then running these six commands should not cause any changes to appear in the database.

```
./Update_Compound_GroupContribution_Energies.py
./Update_Reaction_GroupContribution_Energies.py
./Estimate_Reaction_Reversibility.py GC
./Update_Compound_eQuilibrator_Energies.py
./Update_Reaction_eQuilibrator_Energies.py
./Estimate_Reaction_Reversibility.py EQ
```

These easily run together by running:
```
./Rerun_Thermodynamics.sh
```

## Adding a new energy source

The six update scripts above share a single helper module,
`_thermo_helpers.py`, that encapsulates the bits every source repeats:
path resolution, structure picking, sentinel handling, the `[dg, dge]`
formatter, the iterate-and-save scaffolding, and the two
reaction-aggregation patterns (sum-from-compounds vs. direct lookup).
A new source typically reuses these directly — the per-source script
just supplies a parser, a resolver, and a label.

### A compound source (per-structure)

```python
import _thermo_helpers as th
from BiochemPy import Compounds

LABEL = 'My new source'

# 1. Parse the energy table. parse_two_col_energy_table handles
#    'id<TAB>dg<TAB>dge' extracts; parse_gc_compound_table handles
#    MFAToolkit-style 'MolAnalysis.tbl' files; write your own if neither
#    fits.
energies = th.parse_two_col_energy_table(
    th.thermo_path('MyNewSource', 'Compound_Energies.tbl'))

# 2. Per-compound resolver. Called ONLY when the compound has a
#    structure. Return [dg, dge] to write, or None to skip this compound.
def resolve(cpd, stype, structure, aliases):
    return energies.get(structure)  # or whatever lookup makes sense

# 3. Pick on_no_structure='default' to write the sentinel for
#    structure-less compounds (GC behavior), or 'skip' to leave them
#    untouched.
th.run_compound_update(Compounds(), LABEL, resolve, on_no_structure='default')
```

### A reaction source (sum from compound energies)

```python
import _thermo_helpers as th
from BiochemPy import Compounds, Reactions

LABEL = 'My new source'   # must match the label used on the compounds
th.run_reaction_aggregation_update(Reactions(), Compounds(), LABEL)
```

### A reaction source (precomputed reaction-id table)

```python
import _thermo_helpers as th
from BiochemPy import Reactions

LABEL = 'My new source'
rxn_energies = th.parse_two_col_energy_table(
    th.thermo_path('MyNewSource', 'Reaction_Energies.tbl'))
th.run_reaction_lookup_update(Reactions(), LABEL, rxn_energies)
```

That's it — no scaffolding to copy. Add the new script's invocation to
`Rerun_Thermodynamics.sh` in the order it should run (typically the
compound update before the reaction update, both before
`Estimate_Reaction_Reversibility.py`).

### Reversibility heuristics

`Estimate_Reaction_Reversibility.py` runs after the energies are stored.
The heuristic cascade lives in `estimate_one()`; each rule is a small
helper (`_is_atp_synthase`, `_low_energy_points`, etc.) that can be
imported or replaced without touching the others. The per-reaction
stoichiometry walk is a single helper (`_walk_stoichiometry`) that
returns every accumulator the cascade consumes.

> Two latent bugs in the original stoichiometry walk are preserved
> verbatim in `_walk_stoichiometry` (`cpd` is shadowed by the inner
> phosphate loop; `cpd in rgt` checks the row's dict keys rather than
> its compound id). They are flagged with `NB:` comments. Removing
> either changes the pipeline's reversibility output. See the commit
> message for the refactor that established byte-equivalence.
