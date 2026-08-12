"""
run_ecoli_epm_comparison.py — EPM-only comparative analysis across E. coli
genome-scale reconstruction generations.

Deliberately stops at EPMs (order-0 elementary persistent modules): no
ESPM search, no LP self-maintenance verification, no decomposition, no
Hasse diagram. This is a scope cut made explicit by the user after
diagnosing where time actually goes in the full pipeline on a genome-scale
network (iAF1260, 1471 ERCs): ERC->hierarchy->fundamental relations->EPM
took under 2 seconds; the LP verification step (Stage 5 of
compute_organizations), which rebuilds a full pyCOT sub-network object
per candidate rather than reusing a matrix, is the actual bottleneck and
is out of scope here entirely -- not being computed, not being worked
around.

For each of four successive reconstructions of Escherichia coli str. K-12
MG1655 (e_coli_core -> iAF1260 -> iJO1366 -> iML1515), computes EPMs under
several inflow regimes and reports count + size-distribution statistics,
directly comparable across models because the CORE regimes use only the
7 species common to all four (co2_e, glc__D_e, h2o_e, h_e, nh4_e, o2_e,
pi_e).

The "isolated" (no-inflow) scenario is DELIBERATELY EXCLUDED for the three
larger reconstructions: confirmed on iAF1260 to have a ~20x denser
fundamental-relations graph than any fed scenario (91,322 vs 4,387
fundamental synergies), and the EPM DFS itself did not complete within
several minutes there (unclear whether this is a single pathological
DFS state or genuine combinatorial blow-up -- not diagnosed further,
out of scope for this pass). It remains included for e_coli_core, where
it is fast.

Usage (from repo root):
    python projects/COT_Fundamental_Generators_Exploration/scripts/run_ecoli_epm_comparison.py
"""
from __future__ import annotations

import os
import sys
import time
import csv

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations import (
    build_rndata, compute_ercs, build_hierarchy,
    compute_synergies_basis_first, compute_complementarities,
)
from pyCOT.analysis.organizations.epm import compute_epms

sys.path.insert(0, _here)
from inflow_regime_analysis import build_scenario_txt

OUT_DIR = os.path.join(_repo_root, 'projects', 'COT_Fundamental_Generators_Exploration',
                        'outputs', 'ecoli_epm_comparison')
os.makedirs(OUT_DIR, exist_ok=True)

CORE_FOOD = ['co2_e', 'glc__D_e', 'h2o_e', 'h_e', 'nh4_e', 'o2_e', 'pi_e']

NATIVE_FOOD = {
    'e_coli_core': CORE_FOOD,
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

MODELS = ['e_coli_core', 'iAF1260', 'iJO1366', 'iML1515']

RESULTS_CSV = os.path.join(OUT_DIR, 'results.csv')
REPORT_MD = os.path.join(OUT_DIR, 'report.md')

FIELDS = ['model', 'scenario', 'food_species', 'n_species_total', 'n_reactions',
          'e0_size', 'n_ercs', 'n_persistent_ercs', 'n_fundamental_synergies',
          'n_fundamental_complementarities', 'n_epms', 'epm_size_min',
          'epm_size_mean', 'epm_size_max', 'time_s', 'status']


def _scenarios_for(model: str, include_isolated: bool) -> list[tuple[str, list[str]]]:
    s = []
    if include_isolated:
        s.append(('isolated', []))
    s += [
        ('aerobic_core', CORE_FOOD),
        ('anaerobic_core', [x for x in CORE_FOOD if x != 'o2_e']),
        ('carbon_starvation_core', [x for x in CORE_FOOD if x != 'glc__D_e']),
    ]
    if model != 'e_coli_core':
        s.append(('full_native', NATIVE_FOOD[model]))
    return s


def run_one(model: str, label: str, food_tokens: list[str], writer, f_csv):
    path = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG',
                         f'bigg_{model}.txt')
    with open(path, 'r', encoding='utf-8') as fh:
        base_txt = fh.read()
    txt = build_scenario_txt(base_txt, food_tokens)

    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False,
                                      encoding='utf-8') as fh:
        fh.write(txt)
        tmp_path = fh.name

    row = {'model': model, 'scenario': label,
           'food_species': ','.join(food_tokens) or '(none)'}
    t0 = time.perf_counter()
    try:
        rn = read_txt(tmp_path, exact_names=True)
        rn_data = build_rndata(rn, network_id=f'{model}__{label}')
        ercs = compute_ercs(rn_data, verify=False)
        hier = build_hierarchy(ercs)
        syn = compute_synergies_basis_first(ercs, hier)
        comp = compute_complementarities(ercs, hier, syn)
        epm_res = compute_epms(rn_data, ercs, hier, syn, comp, verbose=False)

        sizes = []
        for m in epm_res.all_epm_masks:
            sizes.append(bin(m | rn_data.E0_mask).count('1'))

        row.update({
            'n_species_total': rn_data.n_species,
            'n_reactions': rn_data.n_reactions,
            'e0_size': bin(rn_data.E0_mask).count('1'),
            'n_ercs': len(ercs),
            'n_persistent_ercs': sum(1 for e in ercs if e.is_persistent()),
            'n_fundamental_synergies': len(syn.fundamental),
            'n_fundamental_complementarities': len(comp.fundamental),
            'n_epms': len(epm_res.all_epm_masks),
            'epm_size_min': min(sizes) if sizes else 0,
            'epm_size_mean': round(sum(sizes) / len(sizes), 1) if sizes else 0,
            'epm_size_max': max(sizes) if sizes else 0,
            'time_s': round(time.perf_counter() - t0, 2),
            'status': 'ok',
        })
        print(f"  {model:12s} / {label:24s}: {row['n_ercs']:5d} ERCs, "
              f"{row['n_epms']:5d} EPMs (sizes {row['epm_size_min']}-{row['epm_size_max']}), "
              f"{row['time_s']:.2f}s", flush=True)
    except Exception as e:
        row.update({'status': f'ERROR: {e}', 'time_s': round(time.perf_counter() - t0, 2)})
        print(f"  {model:12s} / {label:24s}: ERROR: {e}", flush=True)
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass

    writer.writerow(row)
    f_csv.flush()
    return row


def main():
    print(f"[run_ecoli_epm_comparison] writing incrementally to {RESULTS_CSV}", flush=True)
    with open(RESULTS_CSV, 'w', newline='', encoding='utf-8') as f_csv:
        writer = csv.DictWriter(f_csv, fieldnames=FIELDS)
        writer.writeheader()
        for model in MODELS:
            include_isolated = (model == 'e_coli_core')
            for label, food in _scenarios_for(model, include_isolated):
                run_one(model, label, food, writer, f_csv)
    print(f"[run_ecoli_epm_comparison] done -> {RESULTS_CSV}", flush=True)


if __name__ == '__main__':
    main()
