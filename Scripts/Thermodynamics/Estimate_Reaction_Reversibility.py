#!/usr/bin/env python
"""Estimate reaction reversibility (``>``, ``<``, ``=``, or ``?``) from the
stored thermodynamic energies and write it back into the reactions JSON.

The algorithm is a fixed cascade of heuristics: the first one that fires
decides the reversibility. To add or remove a heuristic, add or remove
its entry from ``estimate_one()``. To re-use the building blocks elsewhere
(e.g. a different rule set), import the ``_*`` helpers from this module.

The per-source ``GCC``/``EQU`` notes are no longer consulted; ``GC`` and
``EQ`` runs read directly from ``thermodynamics['Group contribution']`` and
``thermodynamics['eQuilibrator']``. After estimation, the computed direction
is appended to whichever Thermodynamics sublist supplied the energy."""
import sys
sys.path.append('../../Libs/Python/')
from BiochemPy import Reactions
from math import log

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
TEMPERATURE = 298.15
GAS_CONSTANT = 0.0019858775
RT_CONST = TEMPERATURE * GAS_CONSTANT
FARADAY = 0.023061  # kcal/vol gram divided by 1000?
SENTINEL_DG = 10000000

# Intracellular concentration range for the bounded MdeltaG estimate.
CELL_MAX = 0.02
CELL_MIN = 0.00001
CELL_CONC = 0.001

# Compounds we treat specially during the stoichiometry walk.
PROTON = "cpd00067"
WATER = "cpd00001"
CO2 = "cpd00011"
PROTON_WATER = frozenset((PROTON, WATER))
LOW_LOCAL_CONC = frozenset(("cpd00007", "cpd11640"))  # O2, H2
ATPS_REAGENTS = frozenset(("cpd00002", "cpd00008", "cpd00009",
                           "cpd00001", "cpd00067"))
ATP = "cpd00002"

# Phosphate-related compounds, used for the low-energy-points heuristic and
# the ABC transporter check.
PHOSPHATE_IDS = ("cpd00002",   # ATP
                 "cpd00008",   # ADP
                 "cpd00018",   # AMP
                 "cpd00009",   # Pi
                 "cpd00012")   # PPi

# Low-energy compounds, taken from MFAToolkit/Parameters/Defaults.txt.
LOW_ENERGY_CPDS = ("cpd00011",  # CO2
                   "cpd00013",  # NH3
                   "cpd11493",  # ACP
                   "cpd00009",  # Pi
                   "cpd00012",  # PPi
                   "cpd00010",  # CoA
                   "cpd00449",  # Dihydrolipoamide
                   "cpd00242")  # HCO3

# Mapping from the ``GC``/``EQ`` CLI flag to the per-source subkey under
# ``rxn_entry['thermodynamics']``. Adding a source means adding one entry.
DB_LEVEL_LABEL = {
    "GC": "Group contribution",
    "EQ": "eQuilibrator",
}
# Order matters for the no-filter fallback: prefer the eQuilibrator energy
# over the Group-contribution one when both are present, mirroring the
# Update_*_eQuilibrator_Energies.py "EQ overwrites GC" precedence that drives
# the top-level ``deltag``/``deltagerr`` values.
DB_LEVEL_PRIORITY = ("EQ", "GC")


# ---------------------------------------------------------------------------
# Per-reaction analysis
# ---------------------------------------------------------------------------
def _thermo_pair(rxn_entry, label):
    """Return ``[dg, dge]`` from ``thermodynamics[label]`` when present and
    non-sentinel, else ``None``. Lists longer than two elements (an earlier
    estimation run appended its direction) are tolerated."""
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


def _energy_for(rxn_entry, db_level):
    """Resolve ``(dg, dge, source_label)`` for the reaction under ``db_level``.

    The energy *values* always come from the top-level ``deltag``/``deltagerr``
    so the reversibility-report numbers stay byte-identical to the dev-branch
    pre-refactor pipeline. The Thermodynamics key drives only two things:

    - *eligibility*: under ``GC``/``EQ`` the reaction is only processed when
      the matching sublist (``Group contribution`` / ``eQuilibrator``) is
      present with a non-sentinel energy. This replaces the legacy
      ``GCC``/``EQU`` notes check.
    - *source label*: the matching sublist's key, returned for the caller's
      direction-appender so the computed reversibility is restamped back
      into the correct sublist.

    For the unfiltered run, ``source_label`` is the sublist whose first
    entry matches the top-level ``deltag`` exactly, with
    ``DB_LEVEL_PRIORITY`` breaking ties — this mirrors the
    ``Update_*_eQuilibrator_Energies.py`` "EQ overwrites GC" precedence that
    drives the top-level value. ``None`` when no sublist matches.

    Returns ``(None, None, None)`` when no usable energy is available."""
    rxn_dg = rxn_entry['deltag']
    rxn_dge = rxn_entry['deltagerr']
    if rxn_dg is not None:
        rxn_dg = float(rxn_dg)
    if rxn_dge is not None:
        rxn_dge = float(rxn_dge)
    if rxn_dg is None or rxn_dg == SENTINEL_DG:
        return None, None, None

    if db_level:
        label = DB_LEVEL_LABEL[db_level]
        if _thermo_pair(rxn_entry, label) is None:
            return None, None, None
        return rxn_dg, rxn_dge, label

    chosen_label = None
    for level in DB_LEVEL_PRIORITY:
        label = DB_LEVEL_LABEL[level]
        pair = _thermo_pair(rxn_entry, label)
        if pair is not None and abs(pair[0] - rxn_dg) < 1e-9:
            chosen_label = label
            break
    return rxn_dg, rxn_dge, chosen_label


def _has_gc_data(rxn_entry):
    """True iff the reaction carries non-sentinel Group-contribution energy.
    Used by EQ runs to decide whether to fall back to the reversibility a
    prior GC run already wrote."""
    return _thermo_pair(rxn_entry, DB_LEVEL_LABEL["GC"]) is not None


def _incomplete_decision(rxn_entry, db_level):
    """Status when the reaction has no usable energy. EQ runs fall back to
    the existing GC reversibility when the reaction has Group-contribution
    data in its Thermodynamics block (set by an earlier GC run)."""
    status = "Incomplete"
    thermoreversibility = "?"
    if db_level == "EQ" and _has_gc_data(rxn_entry):
        thermoreversibility = rxn_entry["reversibility"]
        status += " (GCC)"
    return status, thermoreversibility


def _walk_stoichiometry(stoichiometry):
    """Single pass that produces every accumulator the downstream heuristics
    need. Keeps the original ordering and per-compound special cases."""
    rct_min = rct_max = 0.0
    pdt_min = pdt_max = 0.0
    rgt_sum = 0.0
    proton_cpts = {}
    phosphates = {}

    for rgt in stoichiometry:
        cpd = rgt['compound']
        cpt = rgt['compartment']
        coeff = float(rgt['coefficient'])

        if cpd == PROTON:
            proton_cpts[cpt] = 1

        # NB: two latent bugs preserved verbatim from the original for
        # output equivalence — DO NOT "clean up" either:
        #   1. ``cpd in rgt`` tests the dict keys of the stoichiometry row
        #      (``compound``, ``coefficient``, ``compartment``, ...), not
        #      its compound id. The condition is therefore always False,
        #      so ``phosphates`` is always empty, making the ABCT and
        #      low-energy-points branches unreachable in practice.
        #   2. The loop variable name is ``cpd``, deliberately shadowing
        #      the outer ``cpd``. After the loop, ``cpd`` is the LAST
        #      value of ``PHOSPHATE_IDS`` (cpd00012, PPi) regardless of
        #      the reagent. This makes the PROTON_WATER skip below a
        #      no-op and the CO2 / LOW_LOCAL_CONC special-concentration
        #      branches unreachable. Renaming this variable changes the
        #      output of the entire pipeline.
        for cpd in PHOSPHATE_IDS:
            if cpd in rgt:
                phosphates.setdefault(cpd, 0.0)
                phosphates[cpd] += coeff

        # (cpd is now PHOSPHATE_IDS[-1], not the reagent's compound id;
        # see the note above.)
        if cpd in PROTON_WATER:
            continue

        # MdeltaG bounds under concentration range
        if coeff < 0:
            rct_min += coeff * log(CELL_MIN)
            rct_max += coeff * log(CELL_MAX)
        else:
            pdt_min += coeff * log(CELL_MIN)
            pdt_max += coeff * log(CELL_MAX)

        # mMdeltaG under fixed local concentration
        local_conc = CELL_CONC
        if cpd == CO2:
            local_conc = 0.0001
        elif cpd in LOW_LOCAL_CONC:
            local_conc = 0.000001
        rgt_sum += coeff * log(local_conc)

    return {
        'rct_min': rct_min, 'rct_max': rct_max,
        'pdt_min': pdt_min, 'pdt_max': pdt_max,
        'rgt_sum': rgt_sum,
        'proton_cpts': proton_cpts,
        'phosphates': phosphates,
    }


def _stored_bounds(rxn_dg, rxn_dge, terms):
    """Min/max stored deltaG including concentration-range terms.
    ``rxn_dg_transport`` is reserved for future use (matches original)."""
    rxn_dg_transport = 0.0
    stored_max = (rxn_dg + rxn_dg_transport + rxn_dge
                  + RT_CONST * terms['pdt_max']
                  + RT_CONST * terms['rct_min'])
    stored_min = (rxn_dg + rxn_dg_transport - rxn_dge
                  + RT_CONST * terms['pdt_min']
                  + RT_CONST * terms['rct_max'])
    return stored_max, stored_min


def _is_atp_synthase(rxn_entry, proton_cpts):
    """ATP synthase: transport, multiple proton compartments, exactly the
    five ATPS reagents involved, and only protons crossing the membrane."""
    if rxn_entry['is_transport'] != 1 or len(proton_cpts) <= 1:
        return False

    cpds_cpts = {}
    for rgt in rxn_entry['stoichiometry']:
        cpds_cpts.setdefault(rgt['compound'], []).append(rgt['compartment'])

    if len(cpds_cpts) != 5:
        return False
    for cpd, cpts in cpds_cpts.items():
        if cpd not in ATPS_REAGENTS:
            return False
        if len(cpts) == 2 and cpd != PROTON:
            return False
    return True


def _abc_transporter_decision(rxn_entry, phosphates):
    """Transport reactions with an ATP coefficient: direction follows the
    sign. (Latent today because ``phosphates`` is always empty — preserved
    for parity with the original.)"""
    if rxn_entry['is_transport'] != 1 or ATP not in phosphates:
        return None
    coeff = phosphates[ATP]
    if coeff < 0:
        rev = ">"
    elif coeff > 0:
        rev = "<"
    else:
        # ATP is itself transported; manually reviewed not to be chemical.
        rev = "="
    return f"ABCT: {coeff}", rev


def _low_energy_points(stoichiometry, phosphates):
    """Score using phosphate spread + low-energy-compound coefficients."""
    points = 0.0
    min_coeff = SENTINEL_DG
    if ATP in phosphates and len(phosphates) > 2:
        for pho_coeff in phosphates.values():
            if pho_coeff < min_coeff:
                min_coeff = pho_coeff
    if min_coeff != SENTINEL_DG:
        points -= abs(min_coeff)

    for rgt in stoichiometry:
        if rgt['compound'] in LOW_ENERGY_CPDS:
            points -= float(rgt['coefficient'])
    return points


def estimate_one(rxn_entry, db_level):
    """Returns ``(status_label, thermoreversibility, source_label)`` for one
    reaction.

    ``source_label`` is the Thermodynamics subkey whose energy fed the
    estimate (``'Group contribution'`` or ``'eQuilibrator'``), or ``None``
    when no estimate ran (empty/incomplete) or when the unfiltered run's
    top-level energy did not match any sublist exactly. The caller uses it
    to append the direction back into the matching sublist.

    The heuristic cascade is intentionally explicit so that adding a new
    rule means inserting one ``if`` branch — and removing one means deleting
    it. Each branch's helper is independently testable."""
    if rxn_entry['status'] == "EMPTY":
        return "Empty", "?", None

    rxn_dg, rxn_dge, source_label = _energy_for(rxn_entry, db_level)
    if rxn_dg is None:
        status, thermoreversibility = _incomplete_decision(rxn_entry, db_level)
        return status, thermoreversibility, None

    terms = _walk_stoichiometry(rxn_entry['stoichiometry'])
    stored_max, stored_min = _stored_bounds(rxn_dg, rxn_dge, terms)

    if stored_max < 0:
        return "MdeltaG(Max): {0:.2f}".format(stored_max), ">", source_label
    if stored_min > 0:
        return "MdeltaG(Min): {0:.2f}".format(stored_min), "<", source_label

    if _is_atp_synthase(rxn_entry, terms['proton_cpts']):
        return "ATPS", "=", source_label

    abct = _abc_transporter_decision(rxn_entry, terms['phosphates'])
    if abct is not None:
        status, thermoreversibility = abct
        return status, thermoreversibility, source_label

    mMdeltaG = rxn_dg + RT_CONST * terms['rgt_sum']
    if -2.0 <= mMdeltaG <= 2.0:
        return "mMdeltaG: {0:.2f}".format(mMdeltaG), "=", source_label

    points = _low_energy_points(rxn_entry['stoichiometry'], terms['phosphates'])
    if points * mMdeltaG > 2:
        if mMdeltaG < 0:
            return ("lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points),
                    ">", source_label)
        if mMdeltaG > 0:
            return ("lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points),
                    "<", source_label)

    return "default", "=", source_label


# ---------------------------------------------------------------------------
# Report writer
# ---------------------------------------------------------------------------
def _write_report(db_level, report):
    """Format matches the original: GC runs drop the original-reversibility
    column from the report, EQ and unfiltered runs keep it."""
    name = "Estimated_Reaction_Reversibility_Report"
    if db_level:
        name += "_" + db_level
    name += ".txt"
    with open(name, "w") as fh:
        for rxn in sorted(report):
            row = list(report[rxn])
            if db_level == "GC":
                del row[1]
            fh.write(rxn + "\t" + "\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _parse_db_level(argv):
    if len(argv) > 1 and argv[1] in ('EQ', 'GC'):
        return argv[1]
    return ''


def _append_direction(rxn_entry, source_label, direction):
    """Re-stamp the picked sublist as ``[dg, dge, direction]``. Truncates
    any prior direction so re-running stays idempotent. No-ops when no
    source label was returned (empty/incomplete reactions, or unfiltered
    runs whose top-level energy did not match a sublist)."""
    if source_label is None:
        return
    thermo = rxn_entry.get('thermodynamics')
    if not isinstance(thermo, dict) or source_label not in thermo:
        return
    pair = thermo[source_label]
    if not pair:
        return
    thermo[source_label] = [pair[0], pair[1], direction]


def main():
    db_level = _parse_db_level(sys.argv)
    helper = Reactions()
    reactions_dict = helper.loadReactions()

    report = {}
    for rxn in sorted(reactions_dict.keys()):
        rxn_entry = reactions_dict[rxn]
        status, thermoreversibility, source_label = estimate_one(
            rxn_entry, db_level)
        report[rxn] = [status, rxn_entry["reversibility"], thermoreversibility]
        rxn_entry['reversibility'] = thermoreversibility
        _append_direction(rxn_entry, source_label, thermoreversibility)

    _write_report(db_level, report)
    print("Saving reactions")
    helper.saveReactions(reactions_dict)


if __name__ == "__main__":
    main()
