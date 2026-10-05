"""
fusigenic_inflow_search.py -- generic per-species inflow sweep (THEORY.md
Section 2.1, "Question 1").

Given a network and a baseline food set, tests each candidate additional
inflow species (added on top of the baseline, one at a time) and reports how
it changes elementary-SO count/size relative to baseline:

    Delta_count = n_elementary(F+c) - n_elementary(F)
    Delta_mean_size, Delta_max_size
    fusion_score = Delta_mean_size * (1 + max(0, -Delta_count))   if Delta_count <= 0
                                                                      and Delta_mean_size > 0
                 = 0                                              otherwise (not fusion by
                                                                      this definition)

    The (1 + max(0, -Delta_count)) factor rewards genuine count REDUCTION
    (classic fusion, as in the E. coli trace-metal case) with a bonus on top
    of the size-growth term, without zeroing the score to 0 when a candidate
    grows an elementary SO without changing count at all (Delta_count == 0) -- an
    earlier version multiplied by Delta_count directly, which incorrectly
    scored every Delta_count==0 candidate as 0 regardless of how much it
    grew mean size; caught by hand-checking the toy host network, where
    fescluster_h (Delta_count=0, Delta_mean=+2.0) was scored 0 by that bug.

Candidates default to every species with a nonzero req_mask across the
network's persistent ERCs (i.e. species that at least one persistent-once-
extended ERC actually needs) union every species that appears as a reactant
somewhere -- in practice, for tractable networks, we just sweep every
species not already in the food set; for large genome-scale networks pass
`max_candidates` to sample instead of exhaustively sweeping (still cheap
since a single elementary-SO computation is seconds, but n_species per-species sweeps
add up).

Usage
-----
    python fusigenic_inflow_search.py <network.txt> <food_sp1,food_sp2,...>
        [--max-candidates N] [--top K]
"""
from __future__ import annotations

import os
import sys
import time
import random

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations import (
    build_rndata, compute_ercs, build_hierarchy,
    compute_synergies_basis_first, compute_complementarities,
)
from pyCOT.analysis.organizations.so_search import compute_elementary_sos

sys.path.insert(0, os.path.join(_repo_root, 'projects',
                                 'COT_Fundamental_Generators_Exploration', 'scripts'))
from inflow_regime_analysis import build_scenario_txt


def _elem_sizes(rn_data, elem_res):
    return [bin(m | rn_data.E0_mask).count('1') for m in elem_res.all_elementary_masks]


def _run_once(base_txt: str, food_tokens: list[str], network_id: str):
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False,
                                      encoding='utf-8') as fh:
        fh.write(build_scenario_txt(base_txt, food_tokens))
        tmp_path = fh.name
    try:
        rn = read_txt(tmp_path, exact_names=True)
        rn_data = build_rndata(rn, network_id=network_id)
        ercs = compute_ercs(rn_data, verify=False)
        hier = build_hierarchy(ercs)
        syn = compute_synergies_basis_first(ercs, hier)
        comp = compute_complementarities(ercs, hier, syn)
        elem_res = compute_elementary_sos(rn_data, ercs, hier, syn, comp, verbose=False)
        sizes = _elem_sizes(rn_data, elem_res)
        return {
            'n_species': rn_data.n_species, 'n_reactions': rn_data.n_reactions,
            'n_ercs': len(ercs), 'n_elementary': len(sizes),
            'mean_size': (sum(sizes) / len(sizes)) if sizes else 0.0,
            'max_size': max(sizes) if sizes else 0,
            'min_size': min(sizes) if sizes else 0,
        }
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass


def search_fusigenic_inflows(
    network_path: str,
    base_food: list[str],
    *,
    candidates: list[str] | None = None,
    max_candidates: int | None = None,
    seed: int = 0,
    verbose: bool = True,
):
    with open(network_path, 'r', encoding='utf-8') as f:
        base_txt = f.read()

    baseline = _run_once(base_txt, base_food, "baseline")
    if verbose:
        print(f"[baseline] food={base_food}  n_elementary={baseline['n_elementary']}  "
              f"mean_size={baseline['mean_size']:.1f}  max_size={baseline['max_size']}  "
              f"({baseline['n_species']} species, {baseline['n_reactions']} reactions, "
              f"{baseline['n_ercs']} ERCs)", flush=True)

    if candidates is None:
        # Every species token appearing anywhere in the raw file, minus the
        # base food set -- the simple, exhaustive-by-default candidate pool.
        import re
        toks = set(re.findall(r'[A-Za-z_][A-Za-z0-9_]*', base_txt))
        # Drop things that are clearly reaction-name-ish or keywords by
        # cross-checking against species actually present after one parse.
        rn0 = read_txt(network_path, exact_names=True)
        rn_data0 = build_rndata(rn0, network_id="probe")
        all_species = set(rn_data0.species_names)
        candidates = sorted(all_species - set(base_food))

    if max_candidates is not None and len(candidates) > max_candidates:
        rng = random.Random(seed)
        candidates = rng.sample(candidates, max_candidates)

    results = []
    t0 = time.perf_counter()
    for i, c in enumerate(candidates):
        r = _run_once(base_txt, base_food + [c], f"cand_{c}")
        d_count = r['n_elementary'] - baseline['n_elementary']
        d_mean = r['mean_size'] - baseline['mean_size']
        d_max = r['max_size'] - baseline['max_size']
        fusion_score = (d_mean * (1 + max(0, -d_count))) if (d_count <= 0 and d_mean > 0) else 0.0
        results.append({
            'species': c, 'n_elementary': r['n_elementary'], 'mean_size': r['mean_size'],
            'max_size': r['max_size'], 'd_count': d_count, 'd_mean': d_mean,
            'd_max': d_max, 'fusion_score': fusion_score,
        })
        if verbose and (i + 1) % max(1, len(candidates) // 10) == 0:
            elapsed = time.perf_counter() - t0
            print(f"  ... {i+1}/{len(candidates)} candidates done ({elapsed:.1f}s)", flush=True)

    results.sort(key=lambda r: r['fusion_score'], reverse=True)
    return baseline, results


def print_report(baseline, results, top: int = 15):
    print(f"\n=== Top {top} fusigenic candidates (by fusion_score = "
          f"Delta_mean_size * (1 + max(0,-Delta_count)), fusion only) ===")
    print(f"{'species':<20} {'n_elementary':>12} {'d_count':>8} {'mean_size':>10} "
          f"{'d_mean':>8} {'max_size':>9} {'d_max':>7} {'fusion_score':>13}")
    for r in results[:top]:
        print(f"{r['species']:<20} {r['n_elementary']:>12} {r['d_count']:>8} "
              f"{r['mean_size']:>10.1f} {r['d_mean']:>+8.1f} {r['max_size']:>9} "
              f"{r['d_max']:>+7} {r['fusion_score']:>13.2f}")

    proliferative = sorted([r for r in results if r['d_count'] > 0],
                            key=lambda r: -r['d_count'])[:5]
    if proliferative:
        print(f"\n=== Most proliferative (Delta_count > 0, NOT fusion) ===")
        for r in proliferative:
            print(f"{r['species']:<20} d_count={r['d_count']:+d}  d_mean={r['d_mean']:+.1f}")

    inert = [r for r in results if r['d_count'] == 0 and abs(r['d_mean']) < 1e-9]
    print(f"\n{len(inert)}/{len(results)} candidates had NO effect on elementary-SO structure at all.")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("network_path")
    ap.add_argument("base_food", help="comma-separated species tokens")
    ap.add_argument("--max-candidates", type=int, default=None)
    ap.add_argument("--candidates", type=str, default=None,
                     help="comma-separated explicit candidate species list "
                          "(use for genome-scale networks -- exhaustive sweep "
                          "is NOT tractable there, see PROGRESS.md)")
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    baseline, results = search_fusigenic_inflows(
        args.network_path, args.base_food.split(','),
        candidates=(args.candidates.split(',') if args.candidates else None),
        max_candidates=args.max_candidates,
    )
    print_report(baseline, results, top=args.top)
