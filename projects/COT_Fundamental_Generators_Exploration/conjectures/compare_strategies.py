"""
compare_strategies.py — run both §6.4 conjectures on the benchmark suite and
record a controlled, apples-to-apples comparison.

For every network in conjectures.benchmark_suite.full_suite():
  1. Compute the shared Stage 0-4 pipeline ONCE (ERCs, hierarchy, fundamental
     synergies, fundamental complementarities) — identical inputs for both
     strategies, so any difference in what they find or how much work they
     do is attributable to the one conceptual difference (vertical lift).
  2. Run run_contained_bfs (Conjecture 1) and run_independent_seed
     (Conjecture 2).
  3. For small networks (|ERCs| <= ORACLE_ERC_LIMIT), also run the
     brute-force oracle (oracles/so_oracle.py) and record whether each
     strategy's result matches it exactly.
  4. Write one CSV row per (network, strategy).

HOW TO RUN
----------
    python projects/COT_Fundamental_Generators_Exploration/conjectures/compare_strategies.py
"""
from __future__ import annotations

import os
import sys
import csv
import time

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from itertools import combinations

from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity import compute_complementarities
from pyCOT.analysis.organizations.so_search import latent_join

from conjectures.benchmark_suite import full_suite, BenchmarkNetwork
from conjectures.strategies import run_contained_bfs, run_independent_seed, StrategyRun
from oracles.so_oracle import so_hierarchy_oracle


def _matches_oracle_modulo_latent_joins(rn_data, found: set[int], oracle_all: set[int]) -> bool:
    """
    so_search.py deliberately never searches for pure "latent joins" --
    unions of already-persistent modules that share no fundamental synergy
    or complementarity edge (see so_search.py's latent_join() docstring).
    So `found` is expected to be a SOUND subset of the oracle's full
    enumeration, with the gap exactly explained by repeatedly applying
    latent_join() to what was found. This mirrors
    tests/test_so_search.py's _assert_matches_oracle_modulo_latent_joins,
    returning a bool instead of asserting.
    """
    if not (found <= oracle_all):
        return False   # unsound: found something the oracle says isn't a real SO

    missing = oracle_all - found
    if not missing:
        return True

    pool = set(found)
    changed = True
    while changed:
        changed = False
        for a, b in combinations(sorted(pool), 2):
            if (a & b) == 0:
                continue
            j = latent_join(rn_data, a, b)
            if j is not None and j not in pool and j in missing:
                pool.add(j)
                missing.discard(j)
                changed = True
    return not missing

# ── Configuration ────────────────────────────────────────────────────────────
MAX_REACTIONS = 100
SAMPLE_CAP = 50            # real BioModels/BiGG networks sampled across the size range
ORACLE_ERC_LIMIT = 20      # brute force is 2^|ERCs| subsets -- keep this small
MAX_SO_ORDER = 15
OUT_DIR = os.path.join(_proj, "outputs", "conjecture_comparison")
OUT_CSV = os.path.join(OUT_DIR, "results.csv")

FIELDNAMES = [
    "network", "category", "n_reactions", "n_species", "n_ercs",
    "n_fundamental_synergies", "n_fundamental_complementarities",
    "strategy",
    "n_so_total", "n_elementary", "max_order_reached",
    "states_explored", "comp_extensions", "syn_extensions",
    "lift_candidates", "canonical_pruned",
    "wall_time_s",
    "matches_oracle", "matches_other_strategy",
    "error",
]


def _run_one_network(bn: BenchmarkNetwork) -> list[dict]:
    rows: list[dict] = []
    try:
        ercs = compute_ercs(bn.rn_data, verify=False)
        hier = build_hierarchy(ercs)
        syn = compute_synergies_basis_first(ercs, hier)
        comp = compute_complementarities(ercs, hier, syn)
    except Exception as exc:
        return [{
            "network": bn.name, "category": bn.category,
            "n_reactions": bn.n_reactions, "n_species": bn.n_species,
            "strategy": "n/a", "error": f"pipeline_error: {exc}",
        }]

    base = {
        "network": bn.name, "category": bn.category,
        "n_reactions": bn.n_reactions, "n_species": bn.n_species,
        "n_ercs": len(ercs),
        "n_fundamental_synergies": len(syn.fundamental),
        "n_fundamental_complementarities": len(comp.fundamental),
    }

    runs: dict[str, StrategyRun] = {}
    for label, fn in (("contained_bfs", run_contained_bfs), ("independent_seed", run_independent_seed)):
        row = dict(base)
        row["strategy"] = label
        try:
            run = fn(bn.rn_data, ercs, hier, syn, comp, max_order=MAX_SO_ORDER)
            runs[label] = run
            row.update({
                "n_so_total": run.n_so_total(),
                "n_elementary": run.n_elementary,
                "max_order_reached": run.max_order_reached,
                "states_explored": run.states_explored,
                "comp_extensions": run.comp_extensions,
                "syn_extensions": run.syn_extensions,
                "lift_candidates": run.lift_candidates,
                "canonical_pruned": run.canonical_pruned,
                "wall_time_s": round(run.wall_time_s, 4),
                "error": "",
            })
        except Exception as exc:
            row["error"] = f"strategy_error: {exc}"
        rows.append(row)

    # ── Cross-strategy agreement ────────────────────────────────────────────
    if "contained_bfs" in runs and "independent_seed" in runs:
        match = set(runs["contained_bfs"].all_so_masks) == set(runs["independent_seed"].all_so_masks)
        for row in rows:
            if row.get("strategy") in ("contained_bfs", "independent_seed"):
                row["matches_other_strategy"] = match

    # ── Oracle cross-validation (small networks only) ───────────────────────
    # compute_elementary_sos additionally reports E0 itself as an extra
    # elementary SO whenever E0_mask != 0 (erc.py's module docstring) -- a
    # deliberate departure from the paper's Def 30 "E_emptyset = emptyset"
    # convention, which the brute-force oracle still encodes literally.
    # Stripped here exactly as tests/test_so_search.py does, rather than
    # taught to the oracle.
    if len(ercs) <= ORACLE_ERC_LIMIT:
        try:
            e0 = bn.rn_data.E0_mask
            oracle_all = {m for m in so_hierarchy_oracle(bn.rn_data, ercs)["all_sos"] if m != e0}
            for row in rows:
                label = row.get("strategy")
                if label in runs:
                    found = {m for m in runs[label].all_so_masks if m != e0}
                    row["matches_oracle"] = _matches_oracle_modulo_latent_joins(
                        bn.rn_data, found, set(oracle_all))
        except Exception as exc:
            for row in rows:
                row["error"] = (row.get("error") or "") + f"  oracle_error: {exc}"

    return rows


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"[compare_strategies] assembling benchmark suite "
          f"(<= {MAX_REACTIONS} reactions, sample_cap={SAMPLE_CAP})...")
    suite = full_suite(max_reactions=MAX_REACTIONS, sample_cap=SAMPLE_CAP)
    print(f"[compare_strategies] {len(suite)} networks "
          f"({sum(1 for b in suite if b.category in ('gold', 'worked_example'))} small/oracle-checkable, "
          f"{sum(1 for b in suite if b.category not in ('gold', 'worked_example'))} real BioModels/BiGG)")

    all_rows: list[dict] = []
    t0 = time.perf_counter()
    for i, bn in enumerate(suite, 1):
        print(f"  [{i:>3}/{len(suite)}] {bn.name:<32} "
              f"rxn={bn.n_reactions:>5} sp={bn.n_species:>5} ...", end="", flush=True)
        rows = _run_one_network(bn)
        all_rows.extend(rows)
        errs = [r["error"] for r in rows if r.get("error")]
        if errs:
            print(f"  ERROR: {errs[0][:80]}")
        else:
            n_ercs = rows[0].get("n_ercs", "?")
            n1 = next((r["n_so_total"] for r in rows if r["strategy"] == "contained_bfs"), "?")
            n2 = next((r["n_so_total"] for r in rows if r["strategy"] == "independent_seed"), "?")
            agree = rows[0].get("matches_other_strategy")
            print(f"  ERCs={n_ercs}  SOs: contained_bfs={n1} independent_seed={n2}  "
                  f"agree={agree if agree is not None else '?'}")

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    total_s = time.perf_counter() - t0
    print(f"\n[compare_strategies] done in {total_s:.1f}s -- wrote {len(all_rows)} rows to {OUT_CSV}")

    # ── Quick summary ────────────────────────────────────────────────────────
    divergences = [
        r for r in all_rows
        if r.get("strategy") == "independent_seed" and r.get("matches_other_strategy") is False
    ]
    oracle_checked = [r for r in all_rows if r.get("matches_oracle") is not None]
    oracle_fails = [r for r in oracle_checked if not r.get("matches_oracle")]
    print(f"[compare_strategies] independent_seed diverged from contained_bfs on "
          f"{len(divergences)} network(s): {[r['network'] for r in divergences]}")
    print(f"[compare_strategies] oracle-checked {len(oracle_checked)} (network, strategy) rows, "
          f"{len(oracle_fails)} mismatch(es)")


if __name__ == "__main__":
    main()
