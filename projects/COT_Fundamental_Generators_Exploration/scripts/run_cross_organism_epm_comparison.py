"""
run_cross_organism_epm_comparison.py -- EPM-only comparison across FOUR
different organisms' BiGG genome-scale reconstructions, extending the
within-E.-coli reconstruction-lineage comparison (run_ecoli_epm_comparison.py)
to a between-organism axis: does the "fundamental generative structure is
size-invariant" finding hold only for E. coli, or more broadly?

Organisms (all real BiGG genome-scale reconstructions, roughly comparable
network size to iAF1260/iJO1366/iML1515):
  - iJN678  Synechocystis sp. PCC 6803 -- cyanobacterium, photoautotroph.
            Native medium uses no3_e (nitrate), not nh4_e, as N source --
            a real, auto-detected physiological difference from E. coli.
  - iYO844  Bacillus subtilis -- Gram-positive (Firmicutes), spore-former.
  - iND750  Saccharomyces cerevisiae -- eukaryote (yeast).
  - iAF987  Geobacter metallireducens -- obligate anaerobic metal-reducer.
            Has NO glc__D_e or o2_e exchange at all (strict anaerobe, grows
            on acetate + Fe(III)/metal reduction) -- cannot run the shared
            aerobic/anaerobic/carbon-starved core scenarios used for the
            other three, so it is run under its own full_native medium only.

Scenario design mirrors run_ecoli_epm_comparison.py exactly: the three
CORE regimes use the same 7 species as the E. coli comparison
(co2_e, glc__D_e, h2o_e, h_e, nh4_e, o2_e, pi_e) so EPM counts/sizes are
directly comparable across the organism axis, not just within E. coli.
full_native uses each organism's own default BiGG exchange set, auto-
extracted from empty-LHS "=> species_e" reactions in the raw file (the
same convention e_coli_core's hand-curated 7-species list already matches
exactly, confirmed before writing this script).

Usage (from repo root):
    python projects/COT_Fundamental_Generators_Exploration/scripts/run_cross_organism_epm_comparison.py
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
                        'outputs', 'cross_organism_epm_comparison')
os.makedirs(OUT_DIR, exist_ok=True)

CORE_FOOD = ['co2_e', 'glc__D_e', 'h2o_e', 'h_e', 'nh4_e', 'o2_e', 'pi_e']

ORGANISMS = {
    'iJN678':  'Synechocystis sp. PCC 6803',
    'iYO844':  'Bacillus subtilis',
    'iND750':  'Saccharomyces cerevisiae',
    'iAF987':  'Geobacter metallireducens',
}

# Organisms that can run the shared 7-species core (need glc__D_e + o2_e
# exchanges present in the model at all). Geobacter has neither.
CORE_CAPABLE = {'iJN678', 'iYO844', 'iND750'}

RESULTS_CSV = os.path.join(OUT_DIR, 'results.csv')

FIELDS = ['model', 'organism', 'scenario', 'food_species', 'n_species_total', 'n_reactions',
          'e0_size', 'n_ercs', 'n_persistent_ercs', 'n_fundamental_synergies',
          'n_fundamental_complementarities', 'n_epms', 'epm_size_min',
          'epm_size_mean', 'epm_size_max', 'time_s', 'status']


def _native_food(path: str) -> list[str]:
    """Auto-extract each model's own default medium: species with a pure
    empty-LHS '=> species_e' inflow reaction in the raw BiGG file. Verified
    to reproduce e_coli_core's hand-curated 7-species list exactly."""
    food = []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            body = line.split(';', 1)[0]
            if ':' not in body or '=>' not in body:
                continue
            lhs = body.split(':', 1)[1].split('=>', 1)[0].strip()
            rhs = body.split('=>', 1)[1].strip()
            if lhs == '' and rhs != '':
                parts = rhs.split()
                sp = parts[-1] if parts else rhs
                if sp.endswith('_e'):
                    food.append(sp)
    return food


def _scenarios_for(model: str, native_food: list[str]) -> list[tuple[str, list[str]]]:
    s = []
    if model in CORE_CAPABLE:
        s += [
            ('aerobic_core', CORE_FOOD),
            ('anaerobic_core', [x for x in CORE_FOOD if x != 'o2_e']),
            ('carbon_starvation_core', [x for x in CORE_FOOD if x != 'glc__D_e']),
        ]
    s.append(('full_native', native_food))
    return s


def run_one(model: str, organism: str, label: str, food_tokens: list[str], writer, f_csv):
    path = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG', f'bigg_{model}.txt')
    with open(path, 'r', encoding='utf-8') as fh:
        base_txt = fh.read()
    txt = build_scenario_txt(base_txt, food_tokens)

    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False,
                                      encoding='utf-8') as fh:
        fh.write(txt)
        tmp_path = fh.name

    row = {'model': model, 'organism': organism, 'scenario': label,
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
        print(f"  {model:10s} ({organism:28s}) / {label:24s}: {row['n_ercs']:5d} ERCs, "
              f"{row['n_epms']:5d} EPMs (sizes {row['epm_size_min']}-{row['epm_size_max']}), "
              f"{row['time_s']:.2f}s", flush=True)
    except Exception as e:
        row.update({'status': f'ERROR: {e}', 'time_s': round(time.perf_counter() - t0, 2)})
        print(f"  {model:10s} ({organism:28s}) / {label:24s}: ERROR: {e}", flush=True)
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass

    writer.writerow(row)
    f_csv.flush()
    return row


def main():
    print(f"[run_cross_organism_epm_comparison] writing incrementally to {RESULTS_CSV}", flush=True)
    with open(RESULTS_CSV, 'w', newline='', encoding='utf-8') as f_csv:
        writer = csv.DictWriter(f_csv, fieldnames=FIELDS)
        writer.writeheader()
        for model, organism in ORGANISMS.items():
            path = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG', f'bigg_{model}.txt')
            native = _native_food(path)
            for label, food in _scenarios_for(model, native):
                run_one(model, organism, label, food, writer, f_csv)
    print(f"[run_cross_organism_epm_comparison] done -> {RESULTS_CSV}", flush=True)


if __name__ == '__main__':
    main()
