"""Composable reaction-reversibility cascade: pluggable heuristics + energy sources.

This is the reusable core behind ``Estimate_Reaction_Reversibility.py``. The cascade
is no longer a hardcoded if/elif chain — it is an ordered list of **heuristics** run
against a reaction with a chosen **energy source**:

    heuristic   : (ctx: Context) -> (status_label, operator) | None      # first non-None wins
    energy_source: (rxn_entry)    -> (dg, dge, source_label) | (None, None, label)

Build any rule set by composing a list of heuristics and picking an energy source::

    src   = per_source_energy("eQuilibrator")          # or top_level_energy("EQ"), explicit_energy(dg, dge)
    rules = DEFAULT_HEURISTICS + [make_ln_reversibility_index_heuristic(ln_ri_map)]
    status, op, label = run_reversibility(rxn_entry, src, rules)

``DEFAULT_HEURISTICS`` reproduces the historical fixed cascade byte-for-byte (same order,
same status strings), so ``Estimate_Reaction_Reversibility.estimate_one`` /
``reversibility_from_energy`` are unchanged in behaviour. Add a heuristic by writing one
``(ctx) -> tuple|None`` function and inserting it in the list; add an energy source by
writing one ``(rxn_entry) -> (dg, dge, label)`` callable.
"""
from dataclasses import dataclass
from math import log

# ---------------------------------------------------------------------------
# Constants (verbatim from the pre-refactor module)
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

# CLI flag -> per-source subkey under ``rxn_entry['thermodynamics']``.
DB_LEVEL_LABEL = {
    "GC": "Group contribution",
    "EQ": "eQuilibrator",
    "DGP": "dGPredictor",
}
# Legacy per-source completeness flag in ``rxn_entry['notes']`` (dGPredictor has none).
DB_LEVEL_NOTE = {
    "GC": "GCC",
    "EQ": "EQU",
}
# Precedence for the no-filter fallback: EQ over GC, DGP last.
DB_LEVEL_PRIORITY = ("EQ", "GC", "DGP")

# Default reversibility-index cutoff (ln(1000) == |log10 gamma| of 3, Noor 2012).
LN_RI_THRESHOLD = 6.9077552789821


# ---------------------------------------------------------------------------
# Energy / eligibility helpers (verbatim)
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
    """Eligible under ``level`` when the structured sublist carries a
    non-sentinel pair OR the legacy ``DB_LEVEL_NOTE[level]`` flag is present."""
    if _thermo_pair(rxn_entry, DB_LEVEL_LABEL[level]) is not None:
        return True
    note = DB_LEVEL_NOTE.get(level)
    return note is not None and note in rxn_entry["notes"]


def _energy_for(rxn_entry, db_level):
    """Resolve ``(dg, dge, source_label)`` for the reaction under ``db_level``.

    Energy *values* always come from the top-level ``deltag``/``deltagerr`` so the
    reversibility-report numbers stay byte-identical to the pre-refactor pipeline.
    The Thermodynamics key + legacy notes drive eligibility and the append target.
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
    """True iff the reaction has Group-contribution coverage."""
    return _is_source_eligible(rxn_entry, "GC")


def _incomplete_decision(rxn_entry, db_level):
    """Status when the reaction has no usable energy. EQ runs fall back to the
    existing GC reversibility when the reaction has Group-contribution data."""
    status = "Incomplete"
    thermoreversibility = "?"
    if db_level == "EQ" and _has_gc_data(rxn_entry):
        thermoreversibility = rxn_entry["reversibility"]
        status += " (GCC)"
    return status, thermoreversibility


# ---------------------------------------------------------------------------
# Heuristic building blocks (verbatim)
# ---------------------------------------------------------------------------
def _walk_stoichiometry(stoichiometry):
    """Single pass producing every accumulator the heuristics need (H2+H3 fixed)."""
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

        if cpd in PHOSPHATE_IDS:
            phosphates.setdefault(cpd, 0.0)
            phosphates[cpd] += coeff

        if cpd in PROTON_WATER:
            continue

        if coeff < 0:
            rct_min += coeff * log(CELL_MIN)
            rct_max += coeff * log(CELL_MAX)
        else:
            pdt_min += coeff * log(CELL_MIN)
            pdt_max += coeff * log(CELL_MAX)

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
    """Min/max stored deltaG including concentration-range terms."""
    rxn_dg_transport = 0.0
    stored_max = (rxn_dg + rxn_dg_transport + rxn_dge
                  + RT_CONST * terms['pdt_max']
                  + RT_CONST * terms['rct_min'])
    stored_min = (rxn_dg + rxn_dg_transport - rxn_dge
                  + RT_CONST * terms['pdt_min']
                  + RT_CONST * terms['rct_max'])
    return stored_max, stored_min


def _is_atp_synthase(rxn_entry, proton_cpts):
    """ATP synthase: transport, multiple proton compartments, exactly the five
    ATPS reagents involved, and only protons crossing the membrane."""
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
    """Transport reactions with an ATP coefficient: direction follows the sign."""
    if rxn_entry['is_transport'] != 1 or ATP not in phosphates:
        return None
    coeff = phosphates[ATP]
    if coeff < 0:
        rev = ">"
    elif coeff > 0:
        rev = "<"
    else:
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


# ---------------------------------------------------------------------------
# Context + composable heuristics
# ---------------------------------------------------------------------------
@dataclass
class Context:
    """Everything a heuristic may need for one reaction, computed once.

    ``terms`` (the stoichiometry walk) and ``mMdeltaG`` are cached lazily so
    composing many heuristics costs only one walk per reaction."""
    rxn_entry: dict
    dg: float
    dge: float
    rt: float = RT_CONST
    _terms: dict = None

    @property
    def terms(self):
        if self._terms is None:
            self._terms = _walk_stoichiometry(self.rxn_entry['stoichiometry'])
        return self._terms

    @property
    def stored_bounds(self):
        return _stored_bounds(self.dg, self.dge, self.terms)

    @property
    def mMdeltaG(self):
        return self.dg + self.rt * self.terms['rgt_sum']


def stored_bounds_heuristic(ctx):
    """MdeltaG bounds over the concentration range: forward if the max stays <0,
    reverse if the min stays >0."""
    stored_max, stored_min = ctx.stored_bounds
    if stored_max < 0:
        return "MdeltaG(Max): {0:.2f}".format(stored_max), ">"
    if stored_min > 0:
        return "MdeltaG(Min): {0:.2f}".format(stored_min), "<"
    return None


def atp_synthase_heuristic(ctx):
    if _is_atp_synthase(ctx.rxn_entry, ctx.terms['proton_cpts']):
        return "ATPS", "="
    return None


def abc_transporter_heuristic(ctx):
    return _abc_transporter_decision(ctx.rxn_entry, ctx.terms['phosphates'])


def mmdeltag_band_heuristic(ctx):
    mMdeltaG = ctx.mMdeltaG
    if -2.0 <= mMdeltaG <= 2.0:
        return "mMdeltaG: {0:.2f}".format(mMdeltaG), "="
    return None


def low_energy_heuristic(ctx):
    mMdeltaG = ctx.mMdeltaG
    points = _low_energy_points(ctx.rxn_entry['stoichiometry'], ctx.terms['phosphates'])
    if points * mMdeltaG > 2:
        if mMdeltaG < 0:
            return "lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points), ">"
        if mMdeltaG > 0:
            return "lowE: {0:.2f}".format(mMdeltaG) + ":" + str(points), "<"
    return None


def default_heuristic(ctx):
    """Terminal heuristic — always fires, mirroring the cascade's final fallback."""
    return "default", "="


# The historical fixed cascade, as a composable list. Order + status strings are
# load-bearing (the reversibility report is byte-compared in the regression test).
DEFAULT_HEURISTICS = [
    atp_synthase_heuristic,
    abc_transporter_heuristic,
    stored_bounds_heuristic,
    mmdeltag_band_heuristic,
    low_energy_heuristic,
    default_heuristic,
]

# Heuristics that need only stoichiometry/structure (no ΔG). These can still be
# evaluated for reactions that carry no stored energy; the rest require a ΔG and
# are reported as "no-energy" by ``evaluate_all_heuristics`` in that case.
ENERGY_FREE_HEURISTICS = (
    atp_synthase_heuristic,
    abc_transporter_heuristic,
    default_heuristic,
)


def make_ln_reversibility_index_heuristic(ln_ri_by_rxn, threshold=LN_RI_THRESHOLD):
    """OPTIONAL heuristic (not in DEFAULT_HEURISTICS): eQuilibrator's reversibility
    index. ``ln_ri_by_rxn`` maps rxn id -> ln(gamma); |ln gamma| > threshold gives a
    directional call (ln gamma < 0 -> ``>``). Demonstrates composing a new,
    energy-derived heuristic into a custom cascade."""
    def heuristic(ctx):
        ln_ri = ln_ri_by_rxn.get(ctx.rxn_entry['id'])
        if ln_ri is not None and abs(ln_ri) > threshold:
            return "lnRI: {0:.2f}".format(ln_ri), (">" if ln_ri < 0 else "<")
        return None
    return heuristic


# ---------------------------------------------------------------------------
# Pluggable energy sources: (rxn_entry) -> (dg, dge, source_label)
# ---------------------------------------------------------------------------
def top_level_energy(db_level):
    """Default source: top-level ``deltag``/``deltagerr`` gated by ``db_level``
    eligibility (the historical behaviour). ``db_level`` is ``''`` / ``'GC'`` /
    ``'EQ'`` / ``'DGP'``."""
    def resolve(rxn_entry):
        return _energy_for(rxn_entry, db_level)
    return resolve


def per_source_energy(label):
    """Source that returns the reaction's OWN ``thermodynamics[label]`` dG (not the
    canonical top-level deltag). ``label`` is a Thermodynamics subkey, e.g.
    ``'Group contribution'`` / ``'eQuilibrator'`` / ``'dGPredictor'``."""
    def resolve(rxn_entry):
        pair = _thermo_pair(rxn_entry, label)
        if pair is None:
            return None, None, label
        return pair[0], pair[1], label
    return resolve


def explicit_energy(dg, dge):
    """Source wrapping an explicit ``(dg, dge)`` pair (for reversibility_from_energy)."""
    def resolve(rxn_entry):
        return dg, dge, None
    return resolve


# ---------------------------------------------------------------------------
# Cascade runner
# ---------------------------------------------------------------------------
def run_reversibility(rxn_entry, energy_source, heuristics=DEFAULT_HEURISTICS):
    """Resolve energy via ``energy_source`` then run ``heuristics`` (first non-None
    wins). Returns ``(status, operator, source_label)``, or ``(None, None, label)``
    when the energy source yields no usable energy — the caller decides how to handle
    that (EMPTY / incomplete fallbacks live in the orchestrator)."""
    dg, dge, source_label = energy_source(rxn_entry)
    if dg is None:
        return None, None, source_label
    ctx = Context(rxn_entry, float(dg), float(dge))
    for heuristic in heuristics:
        result = heuristic(ctx)
        if result is not None:
            return result[0], result[1], source_label
    return "default", "=", source_label


def evaluate_all_heuristics(rxn_entry, energy_source, heuristics=DEFAULT_HEURISTICS):
    """Run EVERY heuristic on one reaction WITHOUT the cascade short-circuit.

    Unlike :func:`run_reversibility` (first non-``None`` wins), this records what
    *every* heuristic would say for the reaction -- for diagnostics and for
    auditing the cascade order (which rules agree, which disagree, which the
    cascade picks).

    Returns ``(source_label, has_energy, results)`` where ``results`` maps each
    heuristic's ``__name__`` to one of:
      * ``(status, operator)`` -- the heuristic fired,
      * ``None``               -- it ran but abstained,
      * ``"no-energy"``        -- it needs a ΔG and the reaction has none,
      * ``"empty"``            -- the reaction is EMPTY (no chemistry).
    """
    if rxn_entry.get("status") == "EMPTY":
        return None, False, {h.__name__: "empty" for h in heuristics}

    dg, dge, source_label = energy_source(rxn_entry)
    has_energy = dg is not None
    results = {}
    if has_energy:
        ctx = Context(rxn_entry, float(dg), float(dge))
        for h in heuristics:
            results[h.__name__] = h(ctx)
    else:
        # No ΔG available: only the structure-only heuristics can be evaluated;
        # the energy-dependent ones are reported as "no-energy" rather than run
        # against a bogus ΔG of 0.
        ctx = Context(rxn_entry, 0.0, 0.0)
        for h in heuristics:
            results[h.__name__] = h(ctx) if h in ENERGY_FREE_HEURISTICS else "no-energy"
    return source_label, has_energy, results
