"""
run_ecoli_eso_eo_comparison.py -- Elementary Semi-Organization (ESO) and
Elementary Organization (EO) comparative analysis across E. coli
genome-scale reconstruction generations.

Supersedes the retired run_ecoli_epm_comparison.py: same scenarios, same
models, same underlying ERC->hierarchy->fundamental-relations->ESO pipeline
(unchanged; what that script called "EPMs" are ESOs, which are also what
the reworked companion paper calls SO0 -- see eso_eo_analysis.py's module
docstring for the terminology note and the code-level confirmation that
the search is ERC-based, not species-based, throughout). This script adds
exactly one further stage per (model, scenario): for every ESO found,
verify LP self-maintenance (self_maintenance.check_self_maintenance) to
determine whether it is also an Elementary Organization (EO).

That stage was previously scoped OUT of this pipeline because LP
verification was diagnosed as the genome-scale bottleneck of the FULL
compute_organizations() search (which explores a much larger candidate
space: free-species extension and latent-join completion for spurious,
non-elementary organizations). Checking only the already-computed ESO list
is a completely different, much smaller cost: ~0.03s/ESO even on iAF1260
(measured directly before writing this script), i.e. a few seconds to at
most ~1 minute of LP time added per (model, scenario) row -- negligible
next to the ERC/hierarchy stage's own ~1-2 minutes at genome scale.

Usage (from repo root):
    python projects/COT_Fundamental_Generators_Exploration/scripts/run_ecoli_eso_eo_comparison.py
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
from pyCOT.analysis.organizations.so_search import compute_elementary_sos

sys.path.insert(0, _here)
from inflow_regime_analysis import build_scenario_txt
from eso_eo_analysis import compute_eso_eo_summary

OUT_DIR = os.path.join(_repo_root, 'projects', 'COT_Fundamental_Generators_Exploration',
                        'outputs', 'ecoli_eso_eo_comparison')
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

FIELDS = ['model', 'scenario', 'food_species', 'n_species_total', 'n_reactions',
          'e0_size', 'n_ercs', 'n_persistent_ercs', 'n_fundamental_synergies',
          'n_fundamental_complementarities',
          'n_eso', 'eso_size_min', 'eso_size_mean', 'eso_size_max',
          'n_eo', 'eo_fraction', 'eo_size_min', 'eo_size_mean', 'eo_size_max',
          'time_s', 'lp_time_s', 'status']


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
        eso_res = compute_elementary_sos(rn_data, ercs, hier, syn, comp, verbose=False)

        eso_eo = compute_eso_eo_summary(rn, rn_data, eso_res)

        row.update({
            'n_species_total': rn_data.n_species,
            'n_reactions': rn_data.n_reactions,
            'e0_size': bin(rn_data.E0_mask).count('1'),
            'n_ercs': len(ercs),
            'n_persistent_ercs': sum(1 for e in ercs if e.is_persistent()),
            'n_fundamental_synergies': len(syn.fundamental),
            'n_fundamental_complementarities': len(comp.fundamental),
            'n_eso': eso_eo['n_eso'],
            'eso_size_min': min(eso_eo['eso_sizes']) if eso_eo['eso_sizes'] else 0,
            'eso_size_mean': round(sum(eso_eo['eso_sizes']) / len(eso_eo['eso_sizes']), 1) if eso_eo['eso_sizes'] else 0,
            'eso_size_max': max(eso_eo['eso_sizes']) if eso_eo['eso_sizes'] else 0,
            'n_eo': eso_eo['n_eo'],
            'eo_fraction': round(eso_eo['eo_fraction'], 3),
            'eo_size_min': min(eso_eo['eo_sizes']) if eso_eo['eo_sizes'] else 0,
            'eo_size_mean': round(sum(eso_eo['eo_sizes']) / len(eso_eo['eo_sizes']), 1) if eso_eo['eo_sizes'] else 0,
            'eo_size_max': max(eso_eo['eo_sizes']) if eso_eo['eo_sizes'] else 0,
            'time_s': round(time.perf_counter() - t0, 2),
            'lp_time_s': round(eso_eo['lp_time_s'], 2),
            'status': 'ok',
        })
        print(f"  {model:12s} / {label:24s}: {row['n_ercs']:5d} ERCs, "
              f"{row['n_eso']:5d} ESO (sizes {row['eso_size_min']}-{row['eso_size_max']}) -> "
              f"{row['n_eo']:5d} EO ({row['eo_fraction']*100:.0f}%), "
              f"{row['time_s']:.2f}s (+{row['lp_time_s']:.2f}s LP)", flush=True)
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
    print(f"[run_ecoli_eso_eo_comparison] writing incrementally to {RESULTS_CSV}", flush=True)
    with open(RESULTS_CSV, 'w', newline='', encoding='utf-8') as f_csv:
        writer = csv.DictWriter(f_csv, fieldnames=FIELDS)
        writer.writeheader()
        for model in MODELS:
            include_isolated = (model == 'e_coli_core')
            for label, food in _scenarios_for(model, include_isolated):
                run_one(model, label, food, writer, f_csv)
    print(f"[run_ecoli_eso_eo_comparison] done -> {RESULTS_CSV}", flush=True)


if __name__ == '__main__':
    main()
