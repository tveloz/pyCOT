"""
yeast_compartment_eso_analysis.py -- tests a specific biochemist hypothesis
about S. cerevisiae (iND750): since eukaryotes have a dual evolutionary
origin (archaeal cytoplasmic/nuclear lineage + bacterial mitochondrial
lineage, itself the product of an ancient endosymbiosis), do ESOs of a
real eukaryotic reconstruction separate into compartment-pure cytoplasmic
vs. compartment-pure mitochondrial modules, or do they "fuse" across the
mitochondrial membrane the way trace-metal cofactors fuse otherwise
separate ESOs in E. coli (paper_ecoli_comparative Sec. 3.3)?

This reuses BiGG's own compartment tags (species suffix _c/_n/_r/_g/_v/_x
= various non-mitochondrial internal compartments, _m = mitochondrial,
_e = extracellular/food, excluded from the classification itself) --
no new modelling assumption, just reading metadata iND750 already carries.

Classification per ESO (species with an _e suffix, i.e. food/exchange
species, are ignored for this classification -- they say nothing about
which internal compartment a module belongs to):
  - "cytoplasmic-only" : every non-food species carries a non-_m suffix
  - "mitochondrial-only": every non-food species carries the _m suffix
  - "compartment-hybrid": draws from both

Also builds a supplemented medium: iND750's auto-extracted native food set
(6 species: nh4_e, o2_e, pi_e, so4_e, glc__D_e, h2o_e -- missing co2_e and
h_e relative to the shared 7-species core, and entirely missing exchange
reactions for iron, zinc, manganese, copper, magnesium, calcium, cobalt,
molybdate, folate and niacin -- this reconstruction simply does not
declare boundary reactions for most trace elements a real yeast culture
needs) is extended with every vitamin/mineral exchange iND750 DOES
declare (biotin, riboflavin, thiamin, pantothenate, myo-inositol,
potassium, sodium, sulfate) on top of the 7-species core, to test whether
mitochondrial engagement changes once the model's own richest available
medium is supplied -- while being explicit that full trace-metal
supplementation is not representable in this particular reconstruction.
"""
from __future__ import annotations

import os
import sys
import time

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations import (
    build_rndata, compute_ercs, build_hierarchy,
    compute_synergies_basis_first, compute_complementarities,
)
from pyCOT.analysis.organizations.so_search import compute_elementary_sos
from pyCOT.analysis.organizations.self_maintenance import check_self_maintenance

sys.path.insert(0, _here)
from inflow_regime_analysis import build_scenario_txt

MODEL_PATH = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG', 'bigg_iND750.txt')

CORE_FOOD = ['co2_e', 'glc__D_e', 'h2o_e', 'h_e', 'nh4_e', 'o2_e', 'pi_e']
NATIVE_FOOD = ['nh4_e', 'o2_e', 'pi_e', 'so4_e', 'glc__D_e', 'h2o_e']
# every vitamin/mineral exchange iND750 actually declares, beyond the core
SUPPLEMENT_FOOD = CORE_FOOD + ['so4_e', 'k_e', 'na1_e', 'btn_e', 'ribflv_e', 'thm_e',
                                'pnto__R_e', 'inost_e']

SCENARIOS = {
    'aerobic_core': CORE_FOOD,
    'native (auto-extracted)': NATIVE_FOOD,
    'supplemented (core + available vitamins/minerals)': SUPPLEMENT_FOOD,
}


def compartment_of(sp: str) -> str:
    for suf in ('_c', '_n', '_r', '_g', '_v', '_x'):
        if sp.endswith(suf):
            return 'cyto'
    if sp.endswith('_m'):
        return 'mito'
    if sp.endswith('_e'):
        return 'food'
    return 'other'


def classify(sp_set: set) -> str:
    comps = {compartment_of(sp) for sp in sp_set if compartment_of(sp) != 'food'}
    if comps == {'cyto'}:
        return 'cytoplasmic-only'
    if comps == {'mito'}:
        return 'mitochondrial-only'
    if 'cyto' in comps and 'mito' in comps:
        return 'compartment-hybrid'
    return 'other/empty'


def run_scenario(label: str, food_tokens: list[str]):
    with open(MODEL_PATH, 'r', encoding='utf-8') as fh:
        base_txt = fh.read()
    txt = build_scenario_txt(base_txt, food_tokens)

    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, encoding='utf-8') as fh:
        fh.write(txt)
        tmp_path = fh.name

    t0 = time.perf_counter()
    rn = read_txt(tmp_path, exact_names=True)
    rn_data = build_rndata(rn, network_id=f'iND750__{label}')
    ercs = compute_ercs(rn_data, verify=False)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)
    eso_res = compute_elementary_sos(rn_data, ercs, hier, syn, comp, verbose=False)
    os.remove(tmp_path)

    names = rn_data.species_names
    species_objs = {s.name: s for s in rn.species()}
    counts = {'cytoplasmic-only': 0, 'mitochondrial-only': 0, 'compartment-hybrid': 0, 'other/empty': 0}
    eo_counts = {'cytoplasmic-only': 0, 'mitochondrial-only': 0, 'compartment-hybrid': 0, 'other/empty': 0}
    sizes_by_cat = {'cytoplasmic-only': [], 'mitochondrial-only': [], 'compartment-hybrid': []}

    for mask in eso_res.all_elementary_masks:
        full_mask = mask | rn_data.E0_mask
        sp_set = {names[j] for j in range(rn_data.n_species) if (full_mask >> j) & 1}
        cat = classify(sp_set)
        counts[cat] += 1
        if cat in sizes_by_cat:
            sizes_by_cat[cat].append(len(sp_set))
        sp_list = [species_objs[n] for n in sp_set]
        is_eo, _flux, _prod = check_self_maintenance(sp_list, rn)
        if is_eo:
            eo_counts[cat] += 1

    elapsed = time.perf_counter() - t0
    n_eso = len(eso_res.all_elementary_masks)
    print(f"\n=== iND750 / {label}  (food: {food_tokens}) ===")
    print(f"  {n_eso} ESOs total, {len(ercs)} ERCs, {elapsed:.1f}s")
    print(f"  {'category':<22} {'count':>6} {'%':>6} {'EO':>6} {'EO%':>6} {'mean size':>10}")
    for cat in ['cytoplasmic-only', 'mitochondrial-only', 'compartment-hybrid', 'other/empty']:
        c, e = counts[cat], eo_counts[cat]
        pct = (c / n_eso * 100) if n_eso else 0.0
        epct = (e / c * 100) if c else 0.0
        mean_sz = (sum(sizes_by_cat[cat]) / len(sizes_by_cat[cat])) if sizes_by_cat.get(cat) else 0.0
        print(f"  {cat:<22} {c:>6} {pct:>5.0f}% {e:>6} {epct:>5.0f}% {mean_sz:>10.1f}")

    return {'label': label, 'n_eso': n_eso, 'counts': counts, 'eo_counts': eo_counts,
            'sizes_by_cat': sizes_by_cat}


if __name__ == '__main__':
    results = {}
    for label, food in SCENARIOS.items():
        results[label] = run_scenario(label, food)
