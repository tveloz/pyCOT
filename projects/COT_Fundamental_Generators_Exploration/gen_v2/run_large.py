"""
run_large.py -- targeted, checkpointed gen_v2 runs for the networks that
don't finish inside compare_old_new.py's per-run time budget.

Why this exists, separate from compare_old_new.py: that script's job is
old-vs-new AGREEMENT across a broad sample, re-running everything on
every sweep. Once a network is known to need multiple sessions' worth of
wall-clock time (the genome-scale ones -- iAB_RBC_283, iIT341,
iIS312_Amastigote all confirmed: ERC generation + hierarchy finish in a
few seconds, the exploration itself is what doesn't fit in 300s), there
is no reason to keep re-running the other 17 networks or the old engine
just to push one big one further. This script runs ONLY the new (gen_v2)
engine, against the SAME checkpoint files compare_old_new.py uses
(outputs/gen_v2_checkpoints/<name>.pkl) -- so progress made here counts
the next time compare_old_new.py runs, and vice versa.

No subprocess wrapper, no hard kill: this is a direct, attended run.
Ctrl+C at any point is safe -- the last periodic checkpoint (every
CHECKPOINT_EVERY_S seconds by default) is already on disk, nothing
earlier is lost, and some recent work is simply redone on the next
invocation (cheaply -- see engine.py's module docstring, "Checkpoint /
resume", for why).

HOW TO RUN: edit the CONFIG block below, then just run this file (the
IDE's Run button, or `python -m gen_v2.run_large`) -- no command-line
arguments needed. CLI flags are still accepted if you want them (handy
for one-off overrides without editing the file); anything left
unspecified on the command line falls back to the CONFIG values here.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

_DATA_ROOT = os.path.join(_repo, "data", "biochemical_databases")
CKPT_DIR = os.path.join(_proj, "outputs", "gen_v2_checkpoints")

# ── CONFIG -- edit these, then just run the file ────────────────────────────
# Explicit network names (.txt basename, "bigg_" prefix optional) take
# priority over the size window below. Leave as [] to use the window instead.
NETWORK_NAMES: list[str] = []

# Reaction-count window (inclusive) used when NETWORK_NAMES is empty. The
# default (600-1500) picks up right where compare_old_new.py's 80-800 sweep
# left off: it starts just below 800 so it re-includes the three networks
# that timed out there (iAB_RBC_283=645, iIT341=737, iIS312_Amastigote=784)
# and reaches a bit further into the genome-scale corpus.
MIN_REACTIONS = 600
MAX_REACTIONS = 1500

# Safety cap on how many networks a single "press play" picks up (smallest
# first within the window) -- set to None for no limit (all matches).
MAX_NETWORKS = 30

# Seconds to spend on EACH network before checkpointing and moving to the
# next; None = no limit, run each to completion (or until you Ctrl+C).
TIME_BUDGET_S = 2400

# How often (seconds) to save a checkpoint WHILE a network is still running,
# independent of TIME_BUDGET_S -- bounds how much is ever at risk from an
# unplanned interruption (e.g. closing the IDE).
CHECKPOINT_EVERY_S = 600.0

# True = ignore and discard any existing checkpoint, start every listed
# network from nothing. Leave False to resume prior progress as normal.
FRESH = False

# Seconds to pause after reporting a network that was ALREADY complete when
# this run reached it (so there's time to actually read the numbers before
# the next network's output starts scrolling by).
PAUSE_ON_CACHED_S = 10
# ─────────────────────────────────────────────────────────────────────────


def find_network(name: str) -> str:
    """Same basename-matching rule compare_old_new.py's discover_networks
    uses (bigg_ prefix stripped, first match by that name wins), but with
    no reaction-count range filter -- for NETWORK_NAMES, which may name
    something outside the configured window on purpose."""
    candidates = []
    for root, _dirs, files in os.walk(_DATA_ROOT):
        for fname in files:
            if not fname.endswith(".txt"):
                continue
            stem = fname[:-4]
            if stem.startswith("bigg_"):
                stem = stem[len("bigg_"):]
            if stem == name:
                candidates.append(os.path.join(root, fname))
    if not candidates:
        raise SystemExit(f"No network file found matching '{name}' under {_DATA_ROOT}")
    if len(candidates) > 1:
        print(f"  NOTE: {len(candidates)} files match '{name}', using the first: {candidates[0]}")
    return candidates[0]


def resolve_networks(names: list[str], min_reactions: int | None, max_reactions: int | None,
                      max_networks: int | None) -> list[tuple[str, str]]:
    if names:
        return [(n, find_network(n)) for n in names]

    from gen_v2.compare_old_new import discover_networks
    found = discover_networks(min_reactions if min_reactions is not None else 0,
                               max_reactions if max_reactions is not None else 10 ** 9)
    print(f"{len(found)} networks found with "
          f"{min_reactions if min_reactions is not None else 0}-"
          f"{max_reactions if max_reactions is not None else 'inf'} reactions.")
    if max_networks is not None and len(found) > max_networks:
        print(f"  MAX_NETWORKS={max_networks}: using the {max_networks} smallest "
              f"(raise MAX_NETWORKS, or set it to None, to run more).")
        found = found[:max_networks]
    for n, _p, r in found:
        print(f"    {n:30s} {r} reactions")
    return [(n, p) for n, p, _r in found]


def run_one(name: str, path: str, time_budget_s: float | None, fresh: bool) -> None:
    from gen_v2.engine import explore, load_checkpoint, build_result

    os.makedirs(CKPT_DIR, exist_ok=True)
    ckpt_path = os.path.join(CKPT_DIR, f"{name}.pkl")
    if fresh and os.path.exists(ckpt_path):
        os.remove(ckpt_path)
        print(f"  fresh=True: removed existing checkpoint {ckpt_path}")

    # Peek BEFORE touching the ERC/hierarchy/synergy/complementarity
    # pipeline at all -- a network that's already fully resolved should
    # cost nothing to report on, not just skip the search itself.
    if not fresh:
        peek = load_checkpoint(ckpt_path)
        if peek is not None and peek.complete:
            # getattr, not .network_label: a checkpoint saved before this
            # field existed won't have the attribute at all after unpickling.
            label = getattr(peek, "network_label", "") or name
            result = build_result(peek.discovered_sos, peek.stats, peek.parent_of, complete=True)
            print(f"\n=== {label}  -- ALREADY COMPLETE (skipped rebuild entirely) ===")
            print(f"  {len(result.all_so_masks)} semi-organizations "
                  f"({len(result.elementary_masks)} elementary, "
                  f"max order {max(result.so_by_order.keys(), default=0)})")
            print(f"  {result.stats.get('states_resolved', 0)} closed sets resolved")
            print(f"  checkpoint: {ckpt_path}")
            if PAUSE_ON_CACHED_S > 0:
                print(f"  (pausing {PAUSE_ON_CACHED_S}s before moving on -- Ctrl+C to stop here)")
                time.sleep(PAUSE_ON_CACHED_S)
            return

    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations.io_pyCOT import build_rndata
    from pyCOT.analysis.organizations.erc import compute_ercs
    from pyCOT.analysis.organizations.hierarchy import build_hierarchy
    from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
    from pyCOT.analysis.organizations.complementarity import compute_complementarities
    from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph

    print(f"\n=== {name}  ({path}) ===")
    t0 = time.perf_counter()
    rn = build_rndata(read_txt(path), network_id=name)
    print(f"  {rn.n_reactions} reactions, {rn.n_species} species")
    ercs = compute_ercs(rn, verify=False)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)
    print(f"  {len(ercs)} ERCs, {len(syn.fundamental)} fundamental synergies, "
          f"{len(comp.fundamental)} fundamental complementarities  ({time.perf_counter() - t0:.1f}s)")

    network_label = f"{name} ({rn.n_reactions} rxn, {len(ercs)} ERCs)"
    g = FundamentalGraph(ercs, hier, syn, comp)
    result = explore(g, verbose=True, checkpoint_path=ckpt_path,
                      time_budget_s=time_budget_s, checkpoint_every_s=CHECKPOINT_EVERY_S,
                      network_label=network_label)

    status = "COMPLETE" if result.complete else "PARTIAL (time budget reached -- re-run to continue)"
    print(f"\n  {name}: {status}")
    print(f"  {len(result.all_so_masks)} semi-organizations found so far "
          f"({len(result.elementary_masks)} elementary, "
          f"max order {max(result.so_by_order.keys(), default=0)})")
    print(f"  {result.stats.get('states_resolved', 0)} closed sets resolved")
    print(f"  checkpoint: {ckpt_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("networks", nargs="*", default=None,
                     help="explicit network name(s); overrides NETWORK_NAMES/the size window")
    ap.add_argument("--min-reactions", type=int, default=None)
    ap.add_argument("--max-reactions", type=int, default=None)
    ap.add_argument("--max-networks", type=int, default=None)
    ap.add_argument("--time-budget", type=float, default=None,
                     help="seconds per network before checkpointing and moving on")
    ap.add_argument("--fresh", action="store_true",
                     help="discard existing checkpoints and start over")
    args = ap.parse_args()

    names = args.networks if args.networks else NETWORK_NAMES
    min_r = args.min_reactions if args.min_reactions is not None else MIN_REACTIONS
    max_r = args.max_reactions if args.max_reactions is not None else MAX_REACTIONS
    max_n = args.max_networks if args.max_networks is not None else MAX_NETWORKS
    time_budget = args.time_budget if args.time_budget is not None else TIME_BUDGET_S
    fresh = args.fresh or FRESH

    targets = resolve_networks(names, min_r, max_r, max_n)
    if not targets:
        raise SystemExit("No networks to run -- widen MIN_REACTIONS/MAX_REACTIONS or set NETWORK_NAMES.")

    for name, path in targets:
        run_one(name, path, time_budget, fresh)


if __name__ == "__main__":
    main()
