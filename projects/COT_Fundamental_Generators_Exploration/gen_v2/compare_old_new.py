"""
compare_old_new.py -- old (production) vs new (gen_v2) engine comparison
on large real reaction networks, where no brute-force oracle is feasible.

No oracle check is possible at this scale, so old-vs-new agreement is the
only available cross-check: if they ever disagree, that is evidence of a
real bug in one of them (gen_v2 is the untrusted one by default, but a
disagreement does not by itself say which side is wrong -- it says
"go look"). Agreement is not proof of correctness either; it just means
nothing has surfaced yet at this scale.

Network selection
------------------
Scans the full data/biochemical_databases/ corpus (not just the curated
category folders) for .txt reaction networks with MIN_REACTIONS <=
n_reactions <= MAX_REACTIONS, deduplicated by basename. Reports however
many it actually finds (and runs all of them) rather than padding the
count up to TARGET_COUNT if the corpus has fewer.

Safety: each (network, engine) run happens in its own subprocess with a
hard wall-clock timeout (TIMEOUT_S) -- these are real, large networks and
neither engine has been stress-tested at this size before. A run that
times out is recorded as such, not silently skipped, and does not block
the rest of the sweep. Results are written to CSV incrementally (one row
as soon as each run finishes) so a long sweep's progress is never lost.

Checkpoint / resume for the new engine (added 2026-10-04): the genome-
scale networks (iAB_RBC_283, iIT341, iIS312_Amastigote -- 600+ reactions)
confirmed they hit TIMEOUT_S not during ERC/hierarchy construction (both
finish in a few seconds even here) but during gen_v2's own exploration --
a search space wide enough that a single 300s run cannot finish it. The
new engine is given its own internal time budget (TIMEOUT_S minus a
safety margin, so it can save a checkpoint and return cleanly before the
subprocess's hard kill) and a checkpoint file under CKPT_DIR, named by
network_id. A network whose new-engine run did not finish in time is
reported with status "partial" (not "timeout") and its progress is saved
to that checkpoint -- simply re-running this script (or a script that
calls gen_v2.engine.explore with the same checkpoint_path) continues
from exactly where it left off, accumulating more of the search on each
run, rather than restarting from nothing every time. This is NOT applied
to the old (production) engine, which stays exactly as it has always
behaved -- out of scope for this exploration project.

Staleness warning: a checkpoint is keyed only by network_id, not by any
hash of gen_v2's own code. If engine.py's search logic changes, delete
CKPT_DIR (or the specific network's .pkl) before re-running -- otherwise
an already-complete checkpoint from before the change is returned as-is,
silently skipping the new logic entirely.

Progress visibility (added 2026-10-04, per request): both engines now run
verbose. The old (production) engine uses its own existing verbose mode
(percentage/rate progress lines, already built into so_search.py). The
new (gen_v2) engine prints a stage marker at each pipeline phase (ERC
generation, hierarchy+fundamental relations, elementary-SO search, then
higher-order exploration ONE ORDER AT A TIME), an immediate, uncapped
line the moment ANY new SO is found, and a rolling window of the last
N_LOG individual ERC-addition steps, refreshed every PRINT_EVERY steps --
so a run that is still going can be read for "where is it right now"
instead of only reporting at the end (or at a timeout). This is a log,
not a terminal redraw: each refresh prints a fresh block, it does not
overwrite previous output in place.

Run:
    python -m gen_v2.compare_old_new
"""
from __future__ import annotations

import os
import sys
import csv
import time
import glob
import multiprocessing as mp

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

_DATA_ROOT = os.path.join(_repo, "data", "biochemical_databases")

# ── Configuration ────────────────────────────────────────────────────────
MIN_REACTIONS = 80
MAX_REACTIONS = 800
TARGET_COUNT  = 20          # requested; actual available corpus may be smaller
TIMEOUT_S     = 300         # per (network, engine) run
MAX_ORDER     = 15
N_LOG         = 30          # rolling trace window size (construction steps)
PRINT_EVERY   = 500         # how often (in construction steps) to flush the trace window
OUT_DIR = os.path.join(_proj, "outputs", "gen_v2_comparison")
OUT_CSV = os.path.join(OUT_DIR, f"old_vs_new_{MIN_REACTIONS}_{MAX_REACTIONS}.csv")
CKPT_DIR = os.path.join(_proj, "outputs", "gen_v2_checkpoints")
CKPT_MARGIN_S = 20   # reserved out of TIMEOUT_S so explore() can save+return before the hard kill

FIELDNAMES = [
    "network", "n_reactions", "n_species", "n_ercs",
    "n_fundamental_synergies", "n_fundamental_complementarities",
    "engine", "status", "wall_time_s",
    "n_elementary", "n_so_total", "max_order_reached",
    "states_explored", "comp_extensions", "syn_extensions", "lift_candidates",
    "error",
]


# ---------------------------------------------------------------------------
# Network discovery
# ---------------------------------------------------------------------------

def _cheap_reaction_count(path: str) -> int:
    n = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            if "=>" in line:
                n += 1
    return n


def discover_networks(min_reactions: int, max_reactions: int) -> list[tuple[str, str, int]]:
    """Full-corpus scan (every folder under data/biochemical_databases/,
    not just the curated category subset), deduplicated by basename."""
    seen: dict[str, tuple[str, int]] = {}
    for root, _dirs, files in os.walk(_DATA_ROOT):
        for fname in sorted(files):
            if not fname.endswith(".txt"):
                continue
            path = os.path.join(root, fname)
            name = os.path.splitext(fname)[0]
            if name.endswith("_manyOrgs"):
                continue
            if name.startswith("bigg_"):
                name = name[len("bigg_"):]
            if name in seen:
                continue
            try:
                n = _cheap_reaction_count(path)
            except Exception:
                continue
            if min_reactions <= n <= max_reactions:
                seen[name] = (path, n)
    return sorted(((n, p, r) for n, (p, r) in seen.items()), key=lambda t: t[2])


# ---------------------------------------------------------------------------
# Subprocess worker -- runs ONE (network, engine) combination, self-contained
# (reloads and rebuilds the shared pipeline itself -- cheap relative to the
# search, avoids pickling FundamentalGraph-derived objects across the
# process boundary).
# ---------------------------------------------------------------------------

def _worker(path: str, network_id: str, engine: str, q: "mp.Queue"):
    # The old engine's existing verbose mode prints Unicode (checkmarks,
    # arrows) that crashes with UnicodeEncodeError under Windows' default
    # console codepage when stdout is captured/redirected rather than
    # attached to a real console. Force UTF-8 with replacement so a
    # display-only issue can never masquerade as an algorithm failure.
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    try:
        t_start = time.perf_counter()

        def _stage(msg: str) -> None:
            print(f"[STAGE t={time.perf_counter() - t_start:7.1f}s] [{engine}] {msg}", flush=True)

        from pyCOT.io.functions import read_txt
        from pyCOT.analysis.organizations.io_pyCOT import build_rndata
        from pyCOT.analysis.organizations.erc import compute_ercs
        from pyCOT.analysis.organizations.hierarchy import build_hierarchy
        from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
        from pyCOT.analysis.organizations.complementarity import compute_complementarities

        rn_pycot = read_txt(path)
        rn = build_rndata(rn_pycot, network_id=network_id)

        _stage(f"ERC generation: starting ({rn.n_reactions} reactions, {rn.n_species} species)")
        ercs = compute_ercs(rn, verify=False)
        _stage(f"ERC generation: done -- {len(ercs)} ERCs")

        hier = build_hierarchy(ercs)
        syn = compute_synergies_basis_first(ercs, hier)
        comp = compute_complementarities(ercs, hier, syn)
        _stage(f"Hierarchy with fundamental relations: done -- "
               f"{len(syn.fundamental)} fundamental synergies, "
               f"{len(comp.fundamental)} fundamental complementarities")

        base = {
            "n_reactions": rn.n_reactions, "n_species": rn.n_species,
            "n_ercs": len(ercs),
            "n_fundamental_synergies": len(syn.fundamental),
            "n_fundamental_complementarities": len(comp.fundamental),
        }

        status = "ok"
        error = ""
        if engine == "old":
            from pyCOT.analysis.organizations.so_search import (
                compute_elementary_sos, compute_so_hierarchy,
            )
            elem = compute_elementary_sos(rn, ercs, hier, syn, comp, verbose=True)
            hr = compute_so_hierarchy(rn, ercs, hier, syn, comp, elem,
                                       use_vertical_lift=True, max_order=MAX_ORDER, verbose=True)
            all_masks = sorted(hr.all_so_masks)
            n_elem = len(elem.all_elementary_masks)
            max_ord = hr.max_order()
            agg = {"states_explored": 0, "comp_extensions": 0, "syn_extensions": 0, "lift_candidates": 0}
            for k in agg:
                agg[k] += elem.stats.get(k, 0)
                for round_stats in hr.stats_by_order.values():
                    agg[k] += round_stats.get(k, 0)
        else:
            from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph
            from gen_v2.engine import explore
            g = FundamentalGraph(ercs, hier, syn, comp)
            os.makedirs(CKPT_DIR, exist_ok=True)
            ckpt_path = os.path.join(CKPT_DIR, f"{network_id}.pkl")
            result = explore(g, verbose=True, n_log=N_LOG, print_every=PRINT_EVERY,
                              checkpoint_path=ckpt_path,
                              time_budget_s=max(TIMEOUT_S - CKPT_MARGIN_S, 1),
                              checkpoint_every_s=CKPT_MARGIN_S,
                              network_label=f"{network_id} ({rn.n_reactions} rxn, {len(ercs)} ERCs)")
            all_masks = sorted(result.all_so_masks)
            n_elem = len(result.elementary_masks)
            max_ord = max(result.so_by_order.keys(), default=0)
            agg = {
                "states_explored": result.stats.get("states_resolved", 0),
                "comp_extensions": 0, "syn_extensions": 0, "lift_candidates": 0,
            }
            if not result.complete:
                status = "partial"
                error = (f"time budget reached -- checkpoint saved at {ckpt_path} "
                         f"({agg['states_explored']} states resolved so far); "
                         f"re-run to continue from there")

        wall = time.perf_counter() - t_start
        q.put({
            "status": status, "error": error, "wall_time_s": round(wall, 3),
            "n_elementary": n_elem, "n_so_total": len(all_masks), "max_order_reached": max_ord,
            "all_masks": all_masks, **base, **agg,
        })
    except Exception as exc:
        import traceback
        q.put({"status": "error", "error": f"{exc}\n{traceback.format_exc()[-500:]}"})


def run_one(path: str, network_id: str, engine: str, timeout_s: int) -> dict:
    q: mp.Queue = mp.Queue()
    p = mp.Process(target=_worker, args=(path, network_id, engine, q))
    p.start()
    p.join(timeout_s)
    if p.is_alive():
        p.terminate()
        p.join()
        return {"status": "timeout", "error": f"exceeded {timeout_s}s"}
    if not q.empty():
        return q.get()
    return {"status": "error", "error": "process exited without a result (crash?)"}


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    candidates = discover_networks(MIN_REACTIONS, MAX_REACTIONS)
    print(f"[compare_old_new] {len(candidates)} networks found with "
          f"{MIN_REACTIONS}-{MAX_REACTIONS} reactions (requested {TARGET_COUNT}).")
    if len(candidates) <= TARGET_COUNT:
        if len(candidates) < TARGET_COUNT:
            print(f"[compare_old_new] NOTE: only {len(candidates)} available in the corpus "
                  f"-- running all of them, not padding to {TARGET_COUNT}.")
        selected = candidates
    else:
        # Evenly-spread sample across the size range (not just the smallest
        # TARGET_COUNT), so the sample still covers the full range requested.
        step = len(candidates) / TARGET_COUNT
        selected = [candidates[int(i * step)] for i in range(TARGET_COUNT)]

    rows_written = 0
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()

        for i, (name, path, n_est) in enumerate(selected, 1):
            print(f"\n[{i:>2}/{len(selected)}] {name}  (~{n_est} reactions)")
            results: dict[str, dict] = {}
            for engine in ("old", "new"):
                print(f"    running {engine}...", end="", flush=True)
                t0 = time.perf_counter()
                res = run_one(path, name, engine, TIMEOUT_S)
                dt = time.perf_counter() - t0
                results[engine] = res
                if res["status"] == "ok":
                    print(f" done in {dt:.1f}s -- {res['n_so_total']} SOs, "
                          f"max_order={res['max_order_reached']}")
                else:
                    print(f" {res['status'].upper()} after {dt:.1f}s ({res.get('error', '')[:80]})")

                row = {"network": name, "engine": engine, "status": res["status"],
                       "error": res.get("error", "")}
                for k in ("n_reactions", "n_species", "n_ercs",
                          "n_fundamental_synergies", "n_fundamental_complementarities",
                          "wall_time_s", "n_elementary", "n_so_total", "max_order_reached",
                          "states_explored", "comp_extensions", "syn_extensions", "lift_candidates"):
                    row[k] = res.get(k, "")
                writer.writerow(row)
                f.flush()
                rows_written += 1

            if results["old"]["status"] == "ok" and results["new"]["status"] == "ok":
                old_masks = set(results["old"]["all_masks"])
                new_masks = set(results["new"]["all_masks"])
                if old_masks == new_masks:
                    print(f"    AGREE: both found {len(old_masks)} SOs, identical sets.")
                else:
                    print(f"    *** DISAGREE *** old-only={len(old_masks - new_masks)}  "
                          f"new-only={len(new_masks - old_masks)}")

    print(f"\n[compare_old_new] wrote {rows_written} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
