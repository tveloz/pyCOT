"""
espm_feasibility_test.py -- quick timing probe: is ESPM search (order >= 1)
tractable at genome scale WITHOUT LP verification?

Background: earlier profiling found the genome-scale bottleneck was LP
self-maintenance verification (Stage 5 of compute_organizations, which
rebuilds a full sub-network object per candidate), NOT the EPM/ESPM
combinatorial search itself. compute_espm() is a pure BFS/DFS extension
over the fundamental-relations graph -- no LP inside it. This script
checks whether that recollection still holds on the current code, on
iAF1260 (the largest-but-one reconstruction), under its aerobic_core
regime (the regime used throughout the EPM comparison).

Usage (from repo root):
    python projects/COT_Fundamental_Generators_Exploration/scripts/espm_feasibility_test.py
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
from pyCOT.analysis.organizations.epm import compute_epms, compute_espm

sys.path.insert(0, _here)
from inflow_regime_analysis import build_scenario_txt

CORE_FOOD = ['co2_e', 'glc__D_e', 'h2o_e', 'h_e', 'nh4_e', 'o2_e', 'pi_e']


def run(model: str, max_order: int):
    path = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG', f'bigg_{model}.txt')
    with open(path, 'r', encoding='utf-8') as fh:
        base_txt = fh.read()
    txt = build_scenario_txt(base_txt, CORE_FOOD)

    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, encoding='utf-8') as fh:
        fh.write(txt)
        tmp_path = fh.name

    t0 = time.perf_counter()
    rn = read_txt(tmp_path, exact_names=True)
    rn_data = build_rndata(rn, network_id=f'{model}__espm_test')
    ercs = compute_ercs(rn_data, verify=False)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)
    t_pre = time.perf_counter()
    print(f"[{model}] pre-EPM stages (ERC/hier/syn/comp): {t_pre - t0:.2f}s "
          f"({len(ercs)} ERCs, {len(syn.fundamental)} fund. synergies, "
          f"{len(comp.fundamental)} fund. complementarities)", flush=True)

    epm_res = compute_epms(rn_data, ercs, hier, syn, comp, verbose=False)
    t_epm = time.perf_counter()
    print(f"[{model}] EPM stage: {t_epm - t_pre:.2f}s -> {len(epm_res.all_epm_masks)} EPMs", flush=True)

    espm_res = compute_espm(rn_data, ercs, hier, syn, comp, epm_res,
                             max_order=max_order, verbose=True)
    t_espm = time.perf_counter()
    print(f"[{model}] ESPM stage (max_order={max_order}): {t_espm - t_epm:.2f}s -> "
          f"{espm_res.total_espm()} ESPMs (order>=1), max order reached = {espm_res.max_order()}, "
          f"{len(espm_res.all_so_masks)} total SOs (EPM+ESPM)", flush=True)
    print(f"[{model}] TOTAL: {t_espm - t0:.2f}s", flush=True)

    try:
        os.remove(tmp_path)
    except OSError:
        pass


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('model', nargs='?', default='iAF1260')
    ap.add_argument('--max-order', type=int, default=2)
    args = ap.parse_args()
    run(args.model, args.max_order)
