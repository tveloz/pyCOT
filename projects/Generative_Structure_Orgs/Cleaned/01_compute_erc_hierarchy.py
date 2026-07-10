#!/usr/bin/env python3
"""
01_compute_erc_hierarchy.py
============================
Batch-compute ERC hierarchy statistics for every reaction network in
SCAN_FOLDERS.  Two CSV files are produced:

  outputs/forest_stats.csv    — one row per network, forest-level metrics
  outputs/per_tree_stats.csv  — one row per ERC hierarchy tree (connected
                                component of the ERC DAG with >= MIN_TREE_NODES)

Both are incremental: already-computed networks are skipped on re-runs.
Run this script before the plot scripts 03 and 04.

forest_stats.csv columns
-------------------------
  file, dataset, n_species, n_reactions, n_ercs, erc_reaction_ratio,
  n_trees, max_depth, mean_depth, mean_tree_depth,
  n_roots, n_leaves, mean_out_degree, max_out_degree, frac_isolated,
  time_ercs, time_forest

per_tree_stats.csv columns
---------------------------
  file, dataset, group, n_nodes, n_ercs, depth, branching
"""

import os
import sys
import time
import statistics

import pandas as pd
import networkx as nx

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, '..', '..', '..', 'src'))
sys.path.insert(0, _SCRIPT_DIR)

from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy
from utils_ercs import load_ercs
from config import (
    SCAN_FOLDERS, MAX_REACTIONS, MAX_ERCS, MIN_ERCS, MIN_TREE_NODES,
    OUT_DIR, VIZ_DIR, assign_group,
)

os.makedirs(OUT_DIR, exist_ok=True)

FOREST_CSV   = os.path.join(OUT_DIR, 'forest_stats.csv')
PER_TREE_CSV = os.path.join(OUT_DIR, 'per_tree_stats.csv')

# ── Helpers ───────────────────────────────────────────────────────────────────

def _node_depths(G):
    """Longest path from any root (in-degree 0) to each node (root = depth 0)."""
    depth = {n: 0 for n in G.nodes()}
    for node in nx.topological_sort(G):
        for child in G.successors(node):
            if depth[node] + 1 > depth[child]:
                depth[child] = depth[node] + 1
    return depth


def forest_metrics(ercs, G):
    """Return dict of per-network forest statistics."""
    n = len(ercs)
    if n == 0:
        return {
            'n_ercs': 0, 'n_trees': 0,
            'max_depth': 0, 'mean_depth': float('nan'), 'mean_tree_depth': float('nan'),
            'n_roots': 0, 'n_leaves': 0,
            'mean_out_degree': float('nan'), 'max_out_degree': 0,
            'frac_isolated': float('nan'),
        }

    components   = list(nx.connected_components(G.to_undirected()))
    depths       = _node_depths(G)
    depth_vals   = list(depths.values())
    max_depth    = max(depth_vals)
    mean_depth   = statistics.mean(depth_vals)

    tree_max = []
    for comp in components:
        sub = G.subgraph(comp)
        sd  = _node_depths(sub)
        tree_max.append(max(sd.values()) if sd else 0)
    mean_tree_depth = statistics.mean(tree_max)

    n_roots  = sum(1 for v in G.nodes() if G.in_degree(v) == 0)
    n_leaves = sum(1 for v in G.nodes() if G.out_degree(v) == 0)
    nonleaf_out = [G.out_degree(v) for v in G.nodes() if G.out_degree(v) > 0]
    mean_out = statistics.mean(nonleaf_out) if nonleaf_out else float('nan')
    max_out  = max(nonleaf_out) if nonleaf_out else 0
    n_iso    = sum(1 for v in G.nodes() if G.degree(v) == 0)

    return {
        'n_ercs':          n,
        'n_trees':         len(components),
        'max_depth':       max_depth,
        'mean_depth':      round(mean_depth, 4),
        'mean_tree_depth': round(mean_tree_depth, 4),
        'n_roots':         n_roots,
        'n_leaves':        n_leaves,
        'mean_out_degree': round(mean_out, 4) if mean_out == mean_out else float('nan'),
        'max_out_degree':  max_out,
        'frac_isolated':   round(n_iso / n, 4),
    }


def per_tree_stats(G_sub):
    """Return (depth, branching, n_nodes) for one connected component."""
    n_nodes    = G_sub.number_of_nodes()
    depths     = _node_depths(G_sub)
    depth      = max(depths.values()) if depths else 0
    out_degs   = [G_sub.out_degree(v) for v in G_sub.nodes() if G_sub.out_degree(v) > 0]
    branching  = float(sum(out_degs) / len(out_degs)) if out_degs else 0.0
    return depth, branching, n_nodes


# ── Load existing CSVs so we can skip already-done networks ──────────────────

def _load_done(csv_path, id_col='file'):
    if os.path.exists(csv_path):
        try:
            df = pd.read_csv(csv_path)
            return set(df[id_col].astype(str))
        except Exception:
            pass
    return set()


forest_done   = _load_done(FOREST_CSV)
per_tree_done = _load_done(PER_TREE_CSV)
already_done  = forest_done & per_tree_done   # skip only if both are complete

forest_existing   = pd.read_csv(FOREST_CSV)   if os.path.exists(FOREST_CSV)   else pd.DataFrame()
per_tree_existing = pd.read_csv(PER_TREE_CSV) if os.path.exists(PER_TREE_CSV) else pd.DataFrame()

print(f"forest_stats.csv:   {len(forest_done)} networks cached")
print(f"per_tree_stats.csv: {len(per_tree_done)} networks cached")
print(f"Will skip {len(already_done)} networks already in both CSVs.\n")

# ── Main loop ─────────────────────────────────────────────────────────────────

forest_rows   = []
per_tree_rows = []
skipped       = []

for folder, dataset_label in SCAN_FOLDERS.items():
    if not os.path.isdir(folder):
        continue
    txt_files = sorted(f for f in os.listdir(folder) if f.endswith('.txt'))
    print(f"\n=== {dataset_label}: {len(txt_files)} files ===")

    for fname in txt_files:
        if fname in already_done:
            print(f"  CACHED  {fname}")
            continue

        txt_path = os.path.join(folder, fname)
        try:
            RN   = read_txt(txt_path)
            n_sp = len(RN.species())
            n_rx = len(RN.reactions())

            if n_rx > MAX_REACTIONS:
                print(f"  SKIP {fname}: {n_rx} reactions > {MAX_REACTIONS}")
                skipped.append((fname, f'too many reactions: {n_rx}'))
                continue

            t0 = time.time()
            ercs, _cached = load_ercs(txt_path, RN, ERC)
            ercs = [e for e in ercs if len(e.get_closure_names(RN)) > 0]
            t_ercs = time.time() - t0

            n_ercs = len(ercs)
            cache_tag = ' [cache]' if _cached else ''
            print(f"  {fname}: {n_ercs} ERCs ({t_ercs:.1f}s){cache_tag}")

            if n_ercs > MAX_ERCS:
                print(f"    SKIP: too many ERCs ({n_ercs})")
                skipped.append((fname, f'too many ERCs: {n_ercs}'))
                continue
            if n_ercs < MIN_ERCS:
                skipped.append((fname, f'too few ERCs: {n_ercs}'))
                continue

            t1 = time.time()
            hierarchy  = ERC_Hierarchy(RN, ercs)
            G          = hierarchy.graph
            metrics    = forest_metrics(ercs, G)
            t_forest   = time.time() - t1

            forest_rows.append({
                'file':               fname,
                'dataset':            dataset_label,
                'n_species':          n_sp,
                'n_reactions':        n_rx,
                'erc_reaction_ratio': round(n_ercs / n_rx, 4) if n_rx > 0 else float('nan'),
                'time_ercs':          round(t_ercs, 3),
                'time_forest':        round(t_forest, 3),
                **metrics,
            })

            # Per-tree stats from connected components
            G_undir = G.to_undirected()
            n_trees_added = 0
            for comp in nx.connected_components(G_undir):
                if len(comp) < MIN_TREE_NODES:
                    continue
                sub = G.subgraph(comp)
                depth, branching, n_nodes = per_tree_stats(sub)
                if depth == 0:
                    continue
                per_tree_rows.append({
                    'file':      fname,
                    'dataset':   dataset_label,
                    'group':     assign_group(dataset_label),
                    'n_nodes':   n_nodes,
                    'n_ercs':    n_ercs,
                    'depth':     depth,
                    'branching': round(branching, 4),
                })
                n_trees_added += 1

            print(f"    trees={metrics['n_trees']}  max_depth={metrics['max_depth']}"
                  f"  qualifying_trees={n_trees_added}")

            # Flush CSVs after each network so progress is not lost on crash
            _fs = pd.concat([forest_existing, pd.DataFrame(forest_rows)], ignore_index=True) \
                  if not forest_existing.empty else pd.DataFrame(forest_rows)
            _pt = pd.concat([per_tree_existing, pd.DataFrame(per_tree_rows)], ignore_index=True) \
                  if not per_tree_existing.empty else pd.DataFrame(per_tree_rows)
            _fs.to_csv(FOREST_CSV,   index=False)
            _pt.to_csv(PER_TREE_CSV, index=False)

        except Exception as exc:
            import traceback
            print(f"  ERR {fname}: {exc}")
            traceback.print_exc()
            skipped.append((fname, str(exc)))

# ── Final save ────────────────────────────────────────────────────────────────

new_forest   = pd.DataFrame(forest_rows)
new_per_tree = pd.DataFrame(per_tree_rows)

df_forest = pd.concat([forest_existing, new_forest], ignore_index=True) \
            if not forest_existing.empty and not new_forest.empty \
            else (new_forest if not new_forest.empty else forest_existing)

df_per_tree = pd.concat([per_tree_existing, new_per_tree], ignore_index=True) \
              if not per_tree_existing.empty and not new_per_tree.empty \
              else (new_per_tree if not new_per_tree.empty else per_tree_existing)

df_forest.to_csv(FOREST_CSV,   index=False)
df_per_tree.to_csv(PER_TREE_CSV, index=False)

print(f"\n=== Done ===")
print(f"forest_stats.csv:   {len(df_forest)} networks total  ({len(new_forest)} new)")
print(f"per_tree_stats.csv: {len(df_per_tree)} trees total  ({len(new_per_tree)} new)")
if skipped:
    print(f"Skipped {len(skipped)} networks:")
    for fn, reason in skipped:
        print(f"  {fn}: {reason}")
