#!/usr/bin/env python
"""Write Group-Contribution compound energies into the BiochemPy compound JSON.

Source: MFAToolkit mol-analysis tables under
``Biochemistry/Thermodynamics/ModelSEED/*_MolAnalysis.tbl``.

Per-compound resolution:
  - structure picked via ``InChIKey`` then ``SMILE`` preference
  - aliases restricted to those listed for the *curated* structure in
    ``All_ModelSEED_Structures.txt``
  - lowest dg among matching aliases wins; default sentinel if none

Writes a value for every compound (sentinel when no structure / no match)."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Compounds
import _thermo_helpers as th

LABEL = 'Group contribution'

compounds_helper = Compounds()
gc_table = th.parse_gc_compound_table(th.thermo_path())
curated_aliases = th.parse_curated_structure_aliases(
    th.structures_path('All_ModelSEED_Structures.txt'))


def resolve(cpd, stype, structure, aliases):
    # A small number of compounds pick a curated structure from
    # Unique_ModelSEED_Structures.txt that the regenerated
    # All_ModelSEED_Structures.txt no longer lists (the two files drifted out of
    # sync when All was not regenerated alongside Unique upstream). Treat a
    # missing provenance record as "no curated aliases": the alias filter then
    # keeps nothing and lowest_energy_gc_style falls back to the sentinel
    # default, preserving the provenance-trim intent instead of raising
    # KeyError.
    curated = curated_aliases.get(cpd, {}).get(structure, {})
    return th.lowest_energy_gc_style(
        (a for a in aliases if a in curated), gc_table)


th.run_compound_update(compounds_helper, LABEL, resolve, on_no_structure='default')
