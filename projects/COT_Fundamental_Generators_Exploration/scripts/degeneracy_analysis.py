"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  degeneracy_analysis.py — How much do generative pathways converge?         ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT IT MEASURES
----------------
At every branching state visited during the EPM search (Mode-1 DFS), the
algorithm proposes N candidate extensions (complementarity producers +
synergy partners that survive the existing canonical-ordering / coverage
gates). Each candidate is a DIFFERENT ERC to add, but after `extend_state`'s
synergy closure runs, several different candidates can land on the exact
same resulting species closure (same `.sp`).

This script measures that convergence -- "generative degeneracy":
  - Per-state: N (candidates tried) vs D (distinct resulting closures).
  - Globally: build a Counter of resulting-closure -> incoming-edge-count
    across the WHOLE traversal, and check whether that distribution is
    Pareto/power-law-like (a small fraction of distinct closures absorbing
    most of the raw candidate edges) -- rank-frequency log-log slope, Gini
    coefficient, and "top X% of closures account for Y% of edges".

Hypothesis under test (user's, paper-relevant): reaction networks that hold
real persistent modules should show HIGH degeneracy (many pathways converge
on few self-productive structures) BECAUSE self-productive structures are
attractors of the generative process almost by definition -- whereas a
network with no real persistent structure would show expansions fanning out
to genuinely distinct closures with little convergence.

WHAT IT DOES NOT DO
--------------------
Does not modify cot_gen/epm.py or fundamental_graph.py. This is read-only
instrumentation built on top of the real FundamentalGraph, wrapping (not
replacing) the real extend_state/erc_syn_close logic, so numbers reflect
exactly what the shipped search actually does.

HOW TO RUN
----------
  1. Edit the CONFIGURATION section below (pick a network, EPM vs ESPM).
  2. Press ▶ (play) in VS Code, or run:
       python projects/COT_fundamental_Generators/scripts/degeneracy_analysis.py
  3. Console prints summary stats; a PNG + CSV are written to
       projects/COT_fundamental_Generators/outputs/degeneracy/<network>/
"""

# ── Path setup (do not edit) ──────────────────────────────────────────────────
from __future__ import annotations
import os, sys, csv, math

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — edit this block, then press Play                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

NETWORK = "e_coli_core"     # short catalogue name or direct .txt path
ANALYZE_ESPM_TOO = True          # also run the (slower, secondary) Mode-2 pass
ESPM_MAX_ORDER = 4               # cap Mode-2 rounds analyzed (can be slow/huge)

# ╔══════════════════════════════════════════════════════════════════════════════╝

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations.io_pyCOT import build_rndata
from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity import compute_complementarities
from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph
from pyCOT.analysis.organizations.epm import _bits, _single_erc_epms


# ---------------------------------------------------------------------------
# Network loading (mirrors run_network.py's resolution logic)
# ---------------------------------------------------------------------------
def _load_network(name: str):
    if os.path.isfile(name):
        path = name
    else:
        candidates = []
        for root, _, files in os.walk(os.path.join(_repo, "data", "biomodels")):
            for f in files:
                if f == f"{name}.txt" or f == f"bigg_{name}.txt":
                    candidates.append(os.path.join(root, f))
        if not candidates:
            raise FileNotFoundError(f"Could not locate a network file for '{name}'")
        path = candidates[0]
    rn_pycot = read_txt(path)
    return build_rndata(rn_pycot, network_id=name), path


# ---------------------------------------------------------------------------
# Instrumented Mode-1 DFS: identical logic to cot_gen/epm.py's _mode1_dfs,
# but additionally records, for every branching state, the full raw list of
# (candidate_erc_idx, resulting_sp) pairs BEFORE distinctness dedup.
# ---------------------------------------------------------------------------
def instrumented_mode1_dfs(seed_states, g, visited_sp, *, state_records, global_edge_counter):
    """
    state_records: list to append (n_candidates, n_distinct, top_share) per
                    branching state (states with >=2 candidates only --
                    "no choice" states aren't informative about degeneracy).
    global_edge_counter: Counter[int] mapping resulting sp -> total raw
                    candidate-edge count across the WHOLE traversal.
    """
    ssm_states, leaf_states = [], []
    stack = [s for s in seed_states if s.sp not in visited_sp]
    while stack:
        state = stack.pop()
        if state.sp in visited_sp:
            continue
        visited_sp.add(state.sp)
        if state.is_ssm:
            ssm_states.append(state)
            continue

        made_progress = False
        min_seed = min(state.erc_set)
        raw_targets: list[int] = []   # resulting .sp per candidate tried

        # Option A: complementarity
        for s_bit in _bits(state.req):
            for prod_idx in g.comp_by_species.get(s_bit, []):
                if prod_idx < state.min_ext:
                    continue
                if (g.species_mask[prod_idx] & state.sp) == g.species_mask[prod_idx]:
                    continue
                new_state = g.extend_state(state, prod_idx)
                newly_added = new_state.erc_set - state.erc_set
                if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                    continue
                raw_targets.append(new_state.sp)
                global_edge_counter[new_state.sp] = global_edge_counter.get(new_state.sp, 0) + 1
                if new_state.sp not in visited_sp:
                    stack.append(new_state)
                    made_progress = True

        # Option B: synergy
        for erc_i in state.erc_set:
            for (j, k) in g.syn_from.get(erc_i, []):
                if j in state.erc_set:
                    continue
                if j < state.min_ext:
                    continue
                if not (g.prod_mask[k] & state.req):
                    continue
                new_state = g.extend_state(state, j)
                newly_added = new_state.erc_set - state.erc_set
                if any(i < min_seed and not g.is_persistent[i] for i in newly_added):
                    continue
                raw_targets.append(new_state.sp)
                global_edge_counter[new_state.sp] = global_edge_counter.get(new_state.sp, 0) + 1
                if new_state.sp not in visited_sp:
                    stack.append(new_state)
                    made_progress = True

        if len(raw_targets) >= 2:
            n = len(raw_targets)
            distinct = len(set(raw_targets))
            top_count = max(raw_targets.count(t) for t in set(raw_targets))
            state_records.append((n, distinct, top_count / n))

        if not made_progress:
            leaf_states.append(state)

    return ssm_states, leaf_states


def run_epm_degeneracy(ercs, hier, syn, comp):
    g = FundamentalGraph(ercs, hier, syn, comp)
    single_idx, single_masks = _single_erc_epms(ercs, hier)
    visited_sp: set[int] = set(single_masks)
    non_p_seeds = [g.make_seed_state(i) for i in range(len(ercs)) if not ercs[i].is_persistent()]

    state_records: list[tuple[int, int, float]] = []
    global_edge_counter: dict[int, int] = {}

    ssm_states, leaf_states = instrumented_mode1_dfs(
        non_p_seeds, g, visited_sp,
        state_records=state_records, global_edge_counter=global_edge_counter,
    )
    return {
        "g": g, "visited_sp": visited_sp, "ssm_states": ssm_states,
        "state_records": state_records, "global_edge_counter": global_edge_counter,
        "n_ssm": len(ssm_states), "n_leaf": len(leaf_states), "n_states": len(visited_sp),
    }


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------
def gini_coefficient(counts: list[int]) -> float:
    if not counts:
        return 0.0
    values = sorted(counts)
    n = len(values)
    cum = 0
    total = sum(values)
    if total == 0:
        return 0.0
    weighted_sum = sum((i + 1) * v for i, v in enumerate(values))
    return (2 * weighted_sum) / (n * total) - (n + 1) / n


def top_k_share(counts: list[int], frac: float) -> float:
    """What fraction of TOTAL edges do the top `frac` of distinct targets absorb?"""
    values = sorted(counts, reverse=True)
    total = sum(values)
    if total == 0:
        return 0.0
    k = max(1, int(round(len(values) * frac)))
    return sum(values[:k]) / total


def power_law_fit_slope(counts: list[int]) -> float | None:
    """Rough log-log rank-frequency slope (classic power-law diagnostic)."""
    values = sorted(counts, reverse=True)
    xs, ys = [], []
    for rank, c in enumerate(values, start=1):
        if c > 0:
            xs.append(math.log(rank))
            ys.append(math.log(c))
    if len(xs) < 5:
        return None
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    den = sum((x - mean_x) ** 2 for x in xs)
    if den == 0:
        return None
    return num / den


def summarize(label: str, state_records, global_edge_counter, out_dir):
    print(f"\n{'='*72}\n{label}\n{'='*72}")
    if not state_records:
        print("  (no branching states with >=2 candidates -- nothing to analyze)")
        return

    ratios = [d / n for (n, d, _) in state_records]
    top_shares = [ts for (_, _, ts) in state_records]
    ns = [n for (n, _, _) in state_records]

    print(f"  branching states analyzed: {len(state_records)}")
    print(f"  candidates per state (N):  min={min(ns)}  median={sorted(ns)[len(ns)//2]}  max={max(ns)}")
    print(f"  distinct/N ratio (D/N):    mean={sum(ratios)/len(ratios):.3f}  "
          f"median={sorted(ratios)[len(ratios)//2]:.3f}  min={min(ratios):.3f}  max={max(ratios):.3f}")
    print(f"  top-target share of N:     mean={sum(top_shares)/len(top_shares):.3f}  "
          f"median={sorted(top_shares)[len(top_shares)//2]:.3f}")

    n_full_degenerate = sum(1 for r in ratios if r < 0.5)
    print(f"  states where <50% of candidates are distinct: {n_full_degenerate}/{len(ratios)} "
          f"({100*n_full_degenerate/len(ratios):.1f}%)")

    counts = list(global_edge_counter.values())
    print(f"\n  GLOBAL target-closure edge distribution ({len(counts)} distinct closures "
          f"received {sum(counts)} total candidate edges):")
    for frac in (0.1, 0.2, 0.5):
        print(f"    top {int(frac*100):>2d}% of distinct closures absorb "
              f"{100*top_k_share(counts, frac):.1f}% of all edges")
    gini = gini_coefficient(counts)
    print(f"    Gini coefficient of edge distribution: {gini:.3f}  (0=uniform, 1=maximally concentrated)")
    slope = power_law_fit_slope(counts)
    if slope is not None:
        print(f"    rank-frequency log-log slope: {slope:.3f}  "
              f"(power-law-like if roughly linear / slope typically in [-3,-0.5])")

    # ── write CSV ──
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"{label.replace(' ', '_')}_state_records.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n_candidates", "n_distinct", "top_target_share"])
        for row in state_records:
            w.writerow(row)
    print(f"\n  wrote {csv_path}")

    # ── plot ──
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

        axes[0].hist(ratios, bins=30, color="#3776ab", edgecolor="white")
        axes[0].set_xlabel("distinct / N  (1.0 = no convergence)")
        axes[0].set_ylabel("number of branching states")
        axes[0].set_title(f"{label}\nPer-state convergence ratio")

        sorted_counts = sorted(counts, reverse=True)
        ranks = list(range(1, len(sorted_counts) + 1))
        axes[1].loglog(ranks, sorted_counts, marker="o", markersize=3, linestyle="none", color="#c0392b")
        axes[1].set_xlabel("rank of target closure (log)")
        axes[1].set_ylabel("incoming candidate edges (log)")
        axes[1].set_title("Rank-frequency (power-law check)")

        cum = 0
        xs_l, ys_l = [0.0], [0.0]
        total = sum(sorted_counts)
        for i, c in enumerate(sorted_counts, start=1):
            cum += c
            xs_l.append(i / len(sorted_counts))
            ys_l.append(cum / total if total else 0)
        axes[2].plot(xs_l, ys_l, color="#27ae60", label="observed")
        axes[2].plot([0, 1], [0, 1], "--", color="gray", label="perfect equality")
        axes[2].set_xlabel("fraction of distinct target closures")
        axes[2].set_ylabel("cumulative fraction of edges")
        axes[2].set_title(f"Lorenz curve (Gini={gini:.3f})")
        axes[2].legend()

        fig.tight_layout()
        png_path = os.path.join(out_dir, f"{label.replace(' ', '_')}.png")
        fig.savefig(png_path, dpi=130)
        plt.close(fig)
        print(f"  wrote {png_path}")
    except ImportError:
        print("  (matplotlib not available -- skipped plot, CSV still written)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 72)
    print(f"Generative Degeneracy Analysis — {NETWORK}")
    print("=" * 72)

    rn, path = _load_network(NETWORK)
    print(f"loaded: {path}")
    ercs = compute_ercs(rn)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn_result=syn)
    print(f"n_ercs={len(ercs)}  n_fundamental_syn={len(syn.fundamental)}  "
          f"n_fundamental_comp={len(comp.fundamental)}")

    out_dir = os.path.join(_proj, "outputs", "degeneracy", NETWORK)

    # ── EPM pass ──
    epm_out = run_epm_degeneracy(ercs, hier, syn, comp)
    print(f"\nEPM search: {epm_out['n_states']} states explored, "
          f"{epm_out['n_ssm']} SSMs, {epm_out['n_leaf']} dead-ends")
    summarize("EPM (Mode-1 DFS)", epm_out["state_records"], epm_out["global_edge_counter"], out_dir)

    # ── ESPM pass (secondary, optional) ──
    if ANALYZE_ESPM_TOO:
        g = epm_out["g"]
        visited_sp = epm_out["visited_sp"]  # continue from where the EPM pass left off

        current_layer = epm_out["ssm_states"]
        espm_state_records: list[tuple[int, int, float]] = []
        espm_edge_counter: dict[int, int] = {}

        for order in range(1, ESPM_MAX_ORDER + 1):
            if not current_layer:
                break
            mode2_seeds = []
            for so_state in current_layer:
                raw_targets = []
                for i in so_state.erc_set:
                    for (j, _) in g.syn_from.get(i, []):
                        if j in so_state.erc_set:
                            continue
                        ext_state = g.extend_state(so_state, j)
                        raw_targets.append(ext_state.sp)
                        espm_edge_counter[ext_state.sp] = espm_edge_counter.get(ext_state.sp, 0) + 1
                        if ext_state.sp not in visited_sp:
                            mode2_seeds.append(ext_state)
                    for a in g.parents[i]:
                        if a in so_state.erc_set:
                            continue
                        if (g.species_mask[a] & so_state.sp) == g.species_mask[a]:
                            continue
                        ext_state = g.extend_state(so_state, a)
                        raw_targets.append(ext_state.sp)
                        espm_edge_counter[ext_state.sp] = espm_edge_counter.get(ext_state.sp, 0) + 1
                        if ext_state.sp not in visited_sp:
                            mode2_seeds.append(ext_state)
                for s_bit in _bits(so_state.prod):
                    for cons_idx in g.comp_consumers_by_species.get(s_bit, []):
                        if cons_idx in so_state.erc_set:
                            continue
                        if (g.species_mask[cons_idx] & so_state.sp) == g.species_mask[cons_idx]:
                            continue
                        ext_state = g.extend_state(so_state, cons_idx)
                        raw_targets.append(ext_state.sp)
                        espm_edge_counter[ext_state.sp] = espm_edge_counter.get(ext_state.sp, 0) + 1
                        if ext_state.sp not in visited_sp:
                            mode2_seeds.append(ext_state)
                if len(raw_targets) >= 2:
                    n = len(raw_targets)
                    distinct = len(set(raw_targets))
                    top_count = max(raw_targets.count(t) for t in set(raw_targets))
                    espm_state_records.append((n, distinct, top_count / n))

            if not mode2_seeds:
                break
            new_ssm_states, _ = instrumented_mode1_dfs(
                mode2_seeds, g, visited_sp,
                state_records=espm_state_records, global_edge_counter=espm_edge_counter,
            )
            print(f"  [ESPM round {order}] {len(current_layer)} SOs -> {len(mode2_seeds)} seeds "
                  f"-> {len(new_ssm_states)} new SOs")
            current_layer = new_ssm_states
            if not new_ssm_states:
                break

        summarize("ESPM (Mode-2 rounds)", espm_state_records, espm_edge_counter, out_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
