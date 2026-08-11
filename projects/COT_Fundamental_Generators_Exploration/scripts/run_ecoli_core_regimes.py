"""
run_ecoli_core_regimes.py — Centler-style inflow-regime comparison on a
real genome-scale-derived E. coli model (BiGG e_coli_core), extending the
reconstruction of Centler et al. 2006 beyond their own hand-built model.

e_coli_core's 7 native inflow (food) species (confirmed this session):
    co2_e, glc__D_e, h2o_e, h_e, nh4_e, o2_e, pi_e

Scenarios chosen to mirror Centler's own "different environments" logic
(starvation / aerobic / anaerobic), the smallest meaningful set that
still varies something biologically real about this specific model:

  aerobic_full        : all 7 native inflows (the model's own default)
  anaerobic           : all 7 minus o2_e (fermentative growth)
  carbon_starvation   : all 7 minus glc__D_e (no carbon/energy source)

Usage (from repo root):
    python projects/COT_Fundamental_Generators_Exploration/scripts/run_ecoli_core_regimes.py
"""
import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))

sys.path.insert(0, _here)
from inflow_regime_analysis import analyze_inflow_regimes

NATIVE_FOOD = ['co2_e', 'glc__D_e', 'h2o_e', 'h_e', 'nh4_e', 'o2_e', 'pi_e']

SCENARIOS = [
    ('aerobic_full', NATIVE_FOOD),
    ('anaerobic', [s for s in NATIVE_FOOD if s != 'o2_e']),
    ('carbon_starvation', [s for s in NATIVE_FOOD if s != 'glc__D_e']),
]

if __name__ == '__main__':
    path = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG',
                         'bigg_e_coli_core.txt')
    analyze_inflow_regimes(
        path, 'e_coli_core',
        SCENARIOS,
        run_organizations=True,
        max_espm_order=2,
    )
