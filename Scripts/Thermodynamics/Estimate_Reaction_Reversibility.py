#!/usr/bin/env python
"""Estimate reaction reversibility (``>``, ``<``, ``=``, or ``?``) from the
stored thermodynamic energies and write it back into the reactions JSON.

The algorithm is a fixed cascade of heuristics: the first one that fires
decides the reversibility. To add or remove a heuristic, add or remove
its entry from ``estimate_one()``. To re-use the building blocks elsewhere
(e.g. a different rule set), import the ``_*`` helpers from this module.

Output is intentionally byte-identical to the pre-refactor script."""
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


# ---------------------------------------------------------------------------
# Per-reaction analysis
# ---------------------------------------------------------------------------
def _is_db_eligible(rxn_entry, db_level):
    """When a ``GC``/``EQ`` filter is active, the reaction must carry a
    ``GCC`` (Group Contribution Complete) or ``EQU`` note matching it."""
    if not db_level:
        return True
    for note in rxn_entry["notes"]:
        if db_level in note and note in ("GCC", "EQU"):
            return True
    return False


def _incomplete_decision(rxn_entry, db_level):
    """Status when the reaction has no usable energy. EQ runs fall back to
    the existing GC reversibility when ``GCC`` is in the notes."""
    status = "Incomplete"
    thermoreversibility = "?"
    if db_level == "EQ" and "GCC" in rxn_entry["notes"]:
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
    """Returns ``(status_label, thermoreversibility)`` for a single reaction.

    The heuristic cascade is intentionally explicit so that adding a new
    rule means inserting one ``if`` branch — and removing one means deleting
    it. Each branch's helper is independently testable."""
    if rxn_entry['status'] == "EMPTY":
        return "Empty", "?"

    rxn_dg = rxn_entry['deltag']
    rxn_dge = rxn_entry['deltagerr']

    if (rxn_dg == SENTINEL_DG or rxn_dg is None
            or not _is_db_eligible(rxn_entry, db_level)):
        return _incomplete_decision(rxn_entry, db_level)

    terms = _walk_stoichiometry(rxn_entry['stoichiometry'])
    stored_max, stored_min = _stored_bounds(rxn_dg, rxn_dge, terms)

    if stored_max < 0:
        return "MdeltaG(Max): {0:.2f}".format(stored_max), ">"
    if stored_min > 0:
        return "MdeltaG(Min): {0:.2f}".format(stored_min), "<"

    if _is_atp_synthase(rxn_entry, terms['proton_cpts']):
        return "ATPS", "="

    abct = _abc_transporter_decision(rxn_entry, terms['phosphates'])
    if abct is not None:
        return abct

    mMdeltaG = rxn_dg + RT_CONST * terms['rgt_sum']
    if -2.0 <= mMdeltaG <= 2.0:
        return "mMdeltaG: {0:.2f}".format(mMdeltaG), "="

    points = _low_energy_points(rxn_entry['stoichiometry'], terms['phosphates'])
    if points * mMdeltaG > 2:
        if mMdeltaG < 0:
            return ("lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points), ">")
        if mMdeltaG > 0:
            return ("lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points), "<")

    return "default", "="


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


def main():
    db_level = _parse_db_level(sys.argv)
    helper = Reactions()
    reactions_dict = helper.loadReactions()

    report = {}
    for rxn in sorted(reactions_dict.keys()):
        rxn_entry = reactions_dict[rxn]
        status, thermoreversibility = estimate_one(rxn_entry, db_level)
        report[rxn] = [status, rxn_entry["reversibility"], thermoreversibility]
        rxn_entry['reversibility'] = thermoreversibility

    _write_report(db_level, report)
    print("Saving reactions")
    helper.saveReactions(reactions_dict)


if __name__ == "__main__":
    main()
