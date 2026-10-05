"""
run_ecoli_reconstruction_regimes.py — Comparative Centler-style inflow-regime
analysis across successive E. coli genome-scale reconstruction generations.

Runs the SAME environmental questions across four reconstructions of
Escherichia coli str. K-12 MG1655, of increasing curation depth:

    e_coli_core (2000, "textbook" core, 95 rxns, 7 native food species)
      -> iAF1260 (2007, 2382 rxns, 22 native food species)
      -> iJO1366 (2011, 2583 rxns, 25 native food species)
      -> iML1515 (2017, 2712 rxns, 24 native food species)

No prior work applying chemical organization theory comparatively across
this reconstruction lineage was found in a literature search this session
-- this is a novel angle, not a replication of an existing study.

Two kinds of scenario, chosen so the comparison means two different things:

  1. CORE regimes (aerobic_core / anaerobic_core / carbon_starvation_core)
     use only the 7 species common to ALL FOUR models (co2_e, glc__D_e,
     h2o_e, h_e, nh4_e, o2_e, pi_e) -- these answer "does 15 years of
     curation change the organizational structure for the SAME
     environmental question", holding the environment fixed and varying
     only reconstruction depth.

  2. FULL_NATIVE regime (the larger models only) uses each model's own
     complete native inflow set (~22-25 species, including trace metals
     e_coli_core doesn't represent at all: cobalamin, molybdate, several
     transition metals) -- this answers "what does the added biological
     detail itself unlock", within one reconstruction.

Fundamental organizations only (include_spurious=False throughout, per
organizations.py's module docstring) -- required at this scale; the
spurious/free-species tail is exponential and was never the point.

Usage (from repo root):
    python projects/COT_Fundamental_Generators_Exploration/scripts/run_ecoli_reconstruction_regimes.py [model_name ...]

    With no arguments, runs all three larger reconstructions
    (e_coli_core was already run separately as the size/speed baseline).
"""
import os
import sys
import time

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))

sys.path.insert(0, _here)
from inflow_regime_analysis import analyze_inflow_regimes

CORE_FOOD = ['co2_e', 'glc__D_e', 'h2o_e', 'h_e', 'nh4_e', 'o2_e', 'pi_e']

NATIVE_FOOD = {
    'iAF1260': ['ca2_e', 'cbl1_e', 'cl_e', 'co2_e', 'cobalt2_e', 'cu2_e', 'fe2_e',
                'fe3_e', 'glc__D_e', 'h2o_e', 'h_e', 'k_e', 'mg2_e', 'mn2_e',
                'mobd_e', 'na1_e', 'nh4_e', 'o2_e', 'pi_e', 'so4_e', 'tungs_e', 'zn2_e'],
    'iJO1366': ['ca2_e', 'cbl1_e', 'cl_e', 'co2_e', 'cobalt2_e', 'cu2_e', 'fe2_e',
                'fe3_e', 'glc__D_e', 'h2o_e', 'h_e', 'k_e', 'mg2_e', 'mn2_e',
                'mobd_e', 'na1_e', 'nh4_e', 'ni2_e', 'o2_e', 'pi_e', 'sel_e',
                'slnt_e', 'so4_e', 'tungs_e', 'zn2_e'],
    'iML1515': ['ca2_e', 'cl_e', 'co2_e', 'cobalt2_e', 'cu2_e', 'fe2_e', 'fe3_e',
                'glc__D_e', 'h2o_e', 'h_e', 'k_e', 'mg2_e', 'mn2_e', 'mobd_e',
                'na1_e', 'nh4_e', 'ni2_e', 'o2_e', 'pi_e', 'sel_e', 'slnt_e',
                'so4_e', 'tungs_e', 'zn2_e'],
}

# Conservative per-model SO-hierarchy order cap -- start low, these networks
# are 15-35x e_coli_core's ERC count (per this session's timing: e_coli_core's
# 49 ERCs -> 0.7s; iAF1260's 1471 ERCs -> ~66s for ERC/hierarchy/relations/
# elementary SOs alone, before the higher-order search/LP). Raise only after
# confirming a given model finishes in reasonable time at the current cap.
MAX_SO_ORDER = {
    'iAF1260': 1,
    'iJO1366': 1,
    'iML1515': 1,
}


def run_one(model: str):
    scenarios = [
        ('aerobic_core', CORE_FOOD),
        ('anaerobic_core', [s for s in CORE_FOOD if s != 'o2_e']),
        ('carbon_starvation_core', [s for s in CORE_FOOD if s != 'glc__D_e']),
        ('full_native', NATIVE_FOOD[model]),
    ]
    path = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG',
                         f'bigg_{model}.txt')
    t0 = time.perf_counter()
    out_dir = analyze_inflow_regimes(
        path, model, scenarios,
        run_organizations=True, run_hasse=True,
        max_so_order=MAX_SO_ORDER[model],
        # The isolated (no-inflow) scenario's fundamental-relations graph
        # can be far denser than any fed scenario -- confirmed on iAF1260:
        # 91,322 fundamental synergies isolated vs. 4,387 fed (~20x), which
        # made the organization search run for hours and consume 17GB+
        # with zero output. Phase A/B (cheap ERC/hierarchy/relations stats,
        # incl. that 91k/4.4k comparison itself) still runs; only the
        # organization search is skipped for this scenario.
        skip_organizations_for=('isolated',),
    )
    print(f"[run_ecoli_reconstruction_regimes] {model} done in "
          f"{time.perf_counter() - t0:.1f}s -> {out_dir}")


if __name__ == '__main__':
    targets = sys.argv[1:] or ['iAF1260', 'iJO1366', 'iML1515']
    for m in targets:
        run_one(m)
