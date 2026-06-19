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

# Mapping from the ``GC``/``EQ``/``DGP`` CLI flag to the per-source subkey
# under ``rxn_entry['thermodynamics']``. Adding a source means adding one
# entry here (and, when a matching legacy note exists, in ``DB_LEVEL_NOTE``
# below — dGPredictor has no legacy note so it is intentionally absent).
DB_LEVEL_LABEL = {
    "GC": "Group contribution",
    "EQ": "eQuilibrator",
    "DGP": "dGPredictor",
}
# Legacy per-source completeness flag in ``rxn_entry['notes']``. A reaction
# is eligible under a CLI flag when EITHER its ``thermodynamics`` sublist
# carries a non-sentinel energy OR the legacy note is present. The OR
# preserves the upstream behavior for the 54 reactions whose pre-existing
# ``GCC`` note is the only completeness signal (the structured sublist
# holds the sentinel because group-contribution coverage was partial).
# dGPredictor has no legacy note: its only completeness signal is the
# structured sublist, so it is intentionally absent from this map and
# ``_is_source_eligible`` falls through to the sublist-only check.
DB_LEVEL_NOTE = {
    "GC": "GCC",
    "EQ": "EQU",
}
# Order matters for the no-filter fallback: prefer the eQuilibrator energy
# over the Group-contribution one when both are present, mirroring the
# Update_*_eQuilibrator_Energies.py "EQ overwrites GC" precedence that drives
# the top-level ``deltag``/``deltagerr`` values. ``DGP`` is appended last
# (additive-only; never overrides the canonical top-level deltag), so it
# is only picked up when its energy happens to match the top-level value
# exactly — a no-op in practice but it keeps the iteration covering all
# sources.
DB_LEVEL_PRIORITY = ("EQ", "GC", "DGP")


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


def _is_source_eligible(rxn_entry, level):
    """A reaction is eligible under ``level`` (``"GC"`` / ``"EQ"`` /
    ``"DGP"``) when EITHER ``thermodynamics[DB_LEVEL_LABEL[level]]`` carries
    a non-sentinel pair OR the legacy ``DB_LEVEL_NOTE[level]`` flag (if one
    exists for this source) is present in ``notes``. Sources without a
    legacy note (e.g. ``DGP``) fall through to the sublist-only check —
    ``DB_LEVEL_NOTE.get`` returns ``None`` and the second clause short-
    circuits."""
    if _thermo_pair(rxn_entry, DB_LEVEL_LABEL[level]) is not None:
        return True
    note = DB_LEVEL_NOTE.get(level)
    return note is not None and note in rxn_entry["notes"]


def _energy_for(rxn_entry, db_level):
    """Resolve ``(dg, dge, source_label)`` for the reaction under ``db_level``.

    The energy *values* always come from the top-level ``deltag``/``deltagerr``
    so the reversibility-report numbers stay byte-identical to the
    pre-refactor pipeline. The Thermodynamics key + legacy notes drive two
    things:

    - *eligibility*: under ``GC``/``EQ`` the reaction is only processed
      when ``_is_source_eligible`` says so (structured sublist OR legacy
      note — see that helper).
    - *source label*: the matching sublist's key, returned for the caller's
      direction-appender so the computed reversibility is restamped back
      into the correct sublist. ``None`` when the only eligibility signal
      was the legacy note and the structured sublist is absent or sentinel
      — no append target exists in that case.

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
        if not _is_source_eligible(rxn_entry, db_level):
            return None, None, None
        label = DB_LEVEL_LABEL[db_level]
        # Only point the append target at the sublist when it actually
        # carries the energy that drove eligibility. Notes-only eligibility
        # leaves the structured sublist untouched (label=None).
        append_label = label if _thermo_pair(rxn_entry, label) is not None else None
        return rxn_dg, rxn_dge, append_label

    chosen_label = None
    for level in DB_LEVEL_PRIORITY:
        label = DB_LEVEL_LABEL[level]
        pair = _thermo_pair(rxn_entry, label)
        if pair is not None and abs(pair[0] - rxn_dg) < 1e-9:
            chosen_label = label
            break
    return rxn_dg, rxn_dge, chosen_label


def _has_gc_data(rxn_entry):
    """True iff the reaction has Group-contribution coverage — either via
    a non-sentinel ``thermodynamics['Group contribution']`` sublist or the
    legacy ``GCC`` note. Used by EQ runs to decide whether to fall back to
    the reversibility a prior GC run already wrote."""
    return _is_source_eligible(rxn_entry, "GC")


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
    need. Keeps the original ordering and per-compound special cases.

    The phosphate accumulator and the CO2 / O2 / H2 concentration overrides
    were both unreachable in the historical code (see the inline note below);
    Heuristics-Review fixes H2 and H3 restore them."""
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

        # Accumulate phosphate-bearing reagents (ATP/ADP/AMP/Pi/PPi) by
        # compound id; this feeds the ABC-transporter rule and the
        # phosphate-spread term of the low-energy-points heuristic.
        #
        # Heuristics-Review H2 + H3 (repaired here): the original code had two
        # coupled typos in this block. (1) the accumulator looped over
        # PHOSPHATE_IDS and tested ``cpd in rgt`` — the stoichiometry row's
        # dict keys (``compound``/``coefficient``/``compartment``), never a
        # compound id — so the condition was always False and ``phosphates``
        # stayed empty, leaving the ABCT and phosphate-spread branches dead
        # (H3). (2) that loop reused the name ``cpd``, shadowing the reagent's
        # compound id with PHOSPHATE_IDS[-1] (PPi) for the rest of the
        # iteration, which silently disabled the per-reagent PROTON_WATER skip
        # and made the CO2 / LOW_LOCAL_CONC concentration overrides below
        # unreachable (the O2/H2 override is H2). Testing the real compound id
        # against PHOSPHATE_IDS repairs all of them at once.
        if cpd in PHOSPHATE_IDS:
            phosphates.setdefault(cpd, 0.0)
            phosphates[cpd] += coeff

        if cpd in PROTON_WATER:
            continue

        # MdeltaG bounds under concentration range
        if coeff < 0:
            rct_min += coeff * log(CELL_MIN)
            rct_max += coeff * log(CELL_MAX)
        else:
            pdt_min += coeff * log(CELL_MIN)
            pdt_max += coeff * log(CELL_MAX)

        # mMdeltaG under fixed local concentration. CO2 sits near 0.1 mM and
        # dissolved O2/H2 near 1 µM, rather than the 1 mM cytoplasmic default.
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
    sign. ATP consumed (negative coeff) drives uptake forward (``>``); ATP
    produced drives it reverse (``<``). Reachable again now that the
    ``phosphates`` accumulator is repaired (Heuristics-Review H3)."""
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


def _cascade(rxn_entry, rxn_dg, rxn_dge):
    """Run the heuristic cascade against an explicit ``(rxn_dg, rxn_dge)``
    pair and return ``(status_label, operator)``.

    Factored out of :func:`estimate_one` so the ``reversibility_from_energy``
    shim can reuse the same logic with per-source (rather than top-level)
    energy values. The eligibility check and the top-level deltag pick are
    intentionally NOT replayed here — callers supply the energies directly,
    matching the upstream ``_estimate_core`` semantics.

    Adding a heuristic still means inserting one ``if`` branch — and removing
    one means deleting it. Each branch's helper is independently testable."""
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
            return ("lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points),
                    ">")
        if mMdeltaG > 0:
            return ("lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points),
                    "<")

    return "default", "="


def estimate_one(rxn_entry, db_level):
    """Returns ``(status_label, thermoreversibility, source_label)`` for one
    reaction.

    ``source_label`` is the Thermodynamics subkey whose energy fed the
    estimate (``'Group contribution'`` or ``'eQuilibrator'``), or ``None``
    when no estimate ran (empty/incomplete) or when the unfiltered run's
    top-level energy did not match any sublist exactly. The caller uses it
    to append the direction back into the matching sublist."""
    if rxn_entry['status'] == "EMPTY":
        return "Empty", "?", None

    rxn_dg, rxn_dge, source_label = _energy_for(rxn_entry, db_level)
    if rxn_dg is None:
        status, thermoreversibility = _incomplete_decision(rxn_entry, db_level)
        return status, thermoreversibility, None

    status, thermoreversibility = _cascade(rxn_entry, rxn_dg, rxn_dge)
    return status, thermoreversibility, source_label


def reversibility_from_energy(rxn_entry, rxn_dg, rxn_dge):
    """Compute the thermodynamic direction operator for a single per-source
    ``(dg, dge)`` pair without consulting the source-eligibility filter or
    the top-level deltag pick.

    Returns one of ``'>'`` / ``'<'`` / ``'='`` / ``'?'``. Used by the
    per-source updaters (``Update_Reaction_dGPredictor_Energies.py``) and by
    the operator backfill (``Add_Reaction_Thermodynamics_Operators.py``) to
    stamp each sublist's own direction.

    Input coercion mirrors the upstream per-source updater:
      * ``rxn_entry['status'] == 'EMPTY'`` -> ``'?'``
      * ``rxn_dg`` that cannot be ``float()``-coerced (``None``, bools,
        ``'nan'``-strings, etc.) -> ``'?'``
      * ``rxn_dg == SENTINEL_DG`` -> ``'?'``
      * ``rxn_dge`` that cannot be coerced -> treated as ``0.0``
    Otherwise the cascade runs and its operator is returned."""
    if isinstance(rxn_entry, dict) and rxn_entry.get('status') == 'EMPTY':
        return '?'

    # Reject bools explicitly: ``float(True) == 1.0`` would otherwise sneak
    # in. Same defensive treatment for ``None`` / non-numeric strings.
    if isinstance(rxn_dg, bool) or rxn_dg is None:
        return '?'
    try:
        dg = float(rxn_dg)
    except (TypeError, ValueError):
        return '?'
    if dg != dg:  # NaN
        return '?'
    if dg == SENTINEL_DG:
        return '?'

    if isinstance(rxn_dge, bool) or rxn_dge is None:
        dge = 0.0
    else:
        try:
            dge = float(rxn_dge)
        except (TypeError, ValueError):
            dge = 0.0
        if dge != dge:  # NaN
            dge = 0.0

    _status, operator = _cascade(rxn_entry, dg, dge)
    return operator


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
    if len(argv) > 1 and argv[1] in ('EQ', 'GC', 'DGP'):
        return argv[1]
    return ''


def main():
    db_level = _parse_db_level(sys.argv)
    helper = Reactions()
    reactions_dict = helper.loadReactions()

    report = {}
    for rxn in sorted(reactions_dict.keys()):
        rxn_entry = reactions_dict[rxn]
        # The third element returned by ``estimate_one`` (the source label
        # of the cascade winner) is intentionally ignored: per-source
        # operators are written at energy-table time by ``_thermo_helpers``
        # (and backfilled by ``Add_Reaction_Thermodynamics_Operators`` for
        # any legacy 2-element entries), each using THAT source's own dG.
        # This step only updates the canonical top-level reversibility.
        status, thermoreversibility, _ = estimate_one(rxn_entry, db_level)
        report[rxn] = [status, rxn_entry["reversibility"], thermoreversibility]
        rxn_entry['reversibility'] = thermoreversibility

    _write_report(db_level, report)
    print("Saving reactions")
    helper.saveReactions(reactions_dict)


if __name__ == "__main__":
    main()
