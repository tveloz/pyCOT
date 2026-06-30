#!/usr/bin/env python3
"""
script_erc_forest_stats.py
==========================
Batch-compute ERC hierarchy forest metrics for every .txt reaction network
found in SCAN_FOLDERS.  Results are saved to a CSV and appended on re-runs
(already-processed networks are skipped unless FORCE_RECOMPUTE=True).

Metrics per network
-------------------
  n_species, n_reactions
  n_ercs                  — total number of ERCs
  erc_reaction_ratio      — n_ercs / n_reactions  (or NaN if n_reactions==0)
  n_trees                 — number of connected components in the *undirected*
                            version of the hierarchy DAG  (= isolated ERCs count
                            as their own tree)
  max_depth               — maximum chain length from any root to a leaf
                            (root = in-degree 0 in the DAG, i.e. maximal ERC)
  mean_depth              — mean depth across all nodes (root = depth 0)
  n_roots                 — ERCs with in-degree 0 (top of hierarchy, maximal)
  n_leaves                — ERCs with out-degree 0 (bottom, minimal)
  mean_out_degree         — mean out-degree among non-leaf nodes (branching)
  max_out_degree          — maximum out-degree  (widest branch point)
  frac_isolated           — fraction of ERCs that are isolated (no edges)
  time_ercs               — wall-clock seconds for ERC computation
  time_forest             — wall-clock seconds for hierarchy + metrics

Scanning
--------
SCAN_FOLDERS is a dict  folder_path -> dataset_label  so each row in the
output CSV carries a 'dataset' tag (e.g. 'BioMD_metabolic', 'BiGG').

The script uses pkl caches for ERC computation exactly as the synergy/
complementarity scripts do: if a .pkl exists next to the .txt it is loaded
instead of recomputed.
"""

import os
import sys
import time
import pickle
import traceback
import statistics

import pandas as pd
import networkx as nx

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))
sys.path.insert(0, _SCRIPT_DIR)

from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from utils_ercs import load_ercs

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
_BIOMD = os.path.join(_PYCOT_ROOT, 'data', 'biomodels')

SCAN_FOLDERS = {
    os.path.join(_BIOMD, 'BioMD_metabolic'):       'BioMD_metabolic',
    os.path.join(_BIOMD, 'BioMD_cell_cycle'):       'BioMD_cell_cycle',
    os.path.join(_BIOMD, 'BioMD_circadian'):        'BioMD_circadian',
    os.path.join(_BIOMD, 'BioMD_signaling'):        'BioMD_signaling',
    os.path.join(_BIOMD, 'BioMD_gene_regulation'):  'BioMD_gene_regulation',
    os.path.join(_BIOMD, 'BioMD_apoptosis'):        'BioMD_apoptosis',
    os.path.join(_BIOMD, 'BioMD_immune'):           'BioMD_immune',
    os.path.join(_BIOMD, 'BioMD_other'):            'BioMD_other',
    os.path.join(_BIOMD, 'BiGG'):                   'BiGG',
    os.path.join(_BIOMD, 'Other'):                  'Other',
}

OUT_DIR  = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', 'outputs', 'forest_stats'))
CSV_FILE = os.path.join(OUT_DIR, 'forest_stats.csv')

MAX_REACTIONS = 1000   # skip networks exceeding this (pre-ERC filter)
MAX_ERCS      = 900   # skip networks exceeding this (post-ERC filter)
MIN_ERCS      = 4     # skip networks with fewer than this many ERCs (post-ERC filter)
MAX_TIME_ERCS = 600   # seconds; abort ERC computation if exceeded

FORCE_RECOMPUTE = False   # set True to ignore the CSV cache

os.makedirs(OUT_DIR, exist_ok=True)


# ---------------------------------------------------------------------------
# ERC loading (with pkl cache)
# ---------------------------------------------------------------------------

def _load_ercs(txt_path, RN):
    """Load ERCs from pkl cache if available, else compute and save."""
    pkl_path = txt_path.replace('.txt', '.pkl')
    if os.path.exists(pkl_path):
        with open(pkl_path, 'rb') as f:
            cache = pickle.load(f)
        raw = cache.get('ERCs', [])
        # raw ERCs are stored as lists [min_generators, closure_species, reactions, label]
        ercs = ERC.from_cache(raw, RN) if hasattr(ERC, 'from_cache') else None
        if ercs is not None:
            return ercs, True  # (ercs, from_cache)
    ercs = ERC.ERCs(RN)
    # Save to pkl
    raw = [[e.min_generators, list(e.get_closure_names(RN)), [], e.label]
           for e in ercs]
    with open(pkl_path, 'wb') as f:
        pickle.dump({'ERCs': raw}, f)
    return ercs, False


def _compute_ercs_timed(txt_path, RN):
    """Returns (ercs, t_ercs, from_cache).  Raises on timeout (approximate)."""
    t0 = time.time()
    pkl_path = txt_path.replace('.txt', '.pkl')
    from_cache = False
    if os.path.exists(pkl_path):
        with open(pkl_path, 'rb') as f:
            cache = pickle.load(f)
        raw = cache.get('ERCs', [])
        if raw:
            ercs = ERC.ERCs(RN)   # still need proper ERC objects for hierarchy
            from_cache = True
            return ercs, time.time() - t0, True
    ercs = ERC.ERCs(RN)
    t_ercs = time.time() - t0
    # cache
    try:
        with open(pkl_path, 'wb') as f:
            pickle.dump({'ERCs': [[e.min_generators,
                                   list(e.get_closure_names(RN)),
                                   [], e.label] for e in ercs]}, f)
    except Exception:
        pass
    return ercs, t_ercs, False


# ---------------------------------------------------------------------------
# Forest metric computation
# ---------------------------------------------------------------------------

def _node_depths(G):
    """
    Return a dict {node: depth} where depth = length of the longest path
    from any root (in-degree 0) down to the node.
    Roots have depth 0.
    Uses topological sort for O(V+E).
    """
    depth = {n: 0 for n in G.nodes()}
    for node in nx.topological_sort(G):
        for child in G.successors(node):
            if depth[node] + 1 > depth[child]:
                depth[child] = depth[node] + 1
    return depth


def forest_metrics(ercs, hierarchy_graph):
    """
    Compute all forest statistics from the ERC list and its hierarchy DAG.

    The ERC hierarchy for a single RN is a FOREST (multiple disconnected
    components / trees).  Per-RN metrics therefore come in two flavours:

      - Global (across the whole forest): max_depth, mean_depth, n_roots, n_leaves
      - Per-tree averages: mean_tree_depth = mean of each component's max depth

    max_depth is dominated by the single deepest component; mean_tree_depth
    gives the typical depth of an individual tree in the forest.

    Parameters
    ----------
    ercs : list of ERC objects
    hierarchy_graph : nx.DiGraph  (edge from parent -> child, larger -> smaller)

    Returns
    -------
    dict with all forest metrics
    """
    G = hierarchy_graph
    n_ercs = len(ercs)

    if n_ercs == 0:
        return {
            'n_ercs': 0, 'erc_reaction_ratio': float('nan'),
            'n_trees': 0, 'max_depth': 0, 'mean_depth': float('nan'),
            'mean_tree_depth': float('nan'),
            'n_roots': 0, 'n_leaves': 0,
            'mean_out_degree': float('nan'), 'max_out_degree': 0,
            'frac_isolated': float('nan'),
        }

    # Connected components (undirected) = number of "trees" in the forest
    G_undir = G.to_undirected()
    components = list(nx.connected_components(G_undir))
    n_trees = len(components)

    # Global depth across the whole forest
    depths = _node_depths(G)
    depth_vals = list(depths.values())
    max_depth  = max(depth_vals)
    mean_depth = statistics.mean(depth_vals) if depth_vals else float('nan')

    # Per-component (per-tree) max depth, then averaged across trees
    # This gives the TYPICAL tree depth, not just the deepest outlier
    tree_max_depths = []
    for component_nodes in components:
        sub = G.subgraph(component_nodes)
        sub_depths = _node_depths(sub)
        tree_max_depths.append(max(sub_depths.values()) if sub_depths else 0)
    mean_tree_depth = statistics.mean(tree_max_depths) if tree_max_depths else float('nan')

    # Roots (in-degree 0) and leaves (out-degree 0)
    n_roots  = sum(1 for n in G.nodes() if G.in_degree(n) == 0)
    n_leaves = sum(1 for n in G.nodes() if G.out_degree(n) == 0)

    # Branching: out-degree of non-leaf nodes (isolated nodes have out-degree 0
    # and are correctly excluded here)
    non_leaf_out = [G.out_degree(n) for n in G.nodes() if G.out_degree(n) > 0]
    mean_out_degree = statistics.mean(non_leaf_out) if non_leaf_out else float('nan')
    max_out_degree  = max(non_leaf_out) if non_leaf_out else 0

    # Isolated nodes (no edges at all)
    n_isolated  = sum(1 for n in G.nodes() if G.degree(n) == 0)
    frac_isolated = n_isolated / n_ercs

    return {
        'n_ercs':           n_ercs,
        'n_trees':          n_trees,
        'max_depth':        max_depth,
        'mean_depth':       round(mean_depth, 4) if mean_depth == mean_depth else float('nan'),
        'mean_tree_depth':  round(mean_tree_depth, 4) if mean_tree_depth == mean_tree_depth else float('nan'),
        'n_roots':          n_roots,
        'n_leaves':         n_leaves,
        'mean_out_degree':  round(mean_out_degree, 4) if mean_out_degree == mean_out_degree else float('nan'),
        'max_out_degree':   max_out_degree,
        'frac_isolated':    round(frac_isolated, 4),
    }


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def main():
    # Load existing CSV to skip already-computed networks.
    # If new metric columns are missing (schema evolved), force full recompute —
    # pkl caches make this fast since ERC computation is skipped.
    _REQUIRED_COLS = {'mean_tree_depth'}
    done = set()
    rows_existing = []
    force = FORCE_RECOMPUTE
    if os.path.exists(CSV_FILE) and not force:
        df_ex = pd.read_csv(CSV_FILE)
        missing_cols = _REQUIRED_COLS - set(df_ex.columns)
        if missing_cols:
            print(f"[INFO] New metric columns detected ({missing_cols}); recomputing all rows.")
            force = True
        else:
            df_ex = df_ex[df_ex['n_ercs'] >= MIN_ERCS]
            done = set(df_ex['file'].tolist())
            rows_existing = df_ex.to_dict('records')
            print(f"Loaded {len(done)} previously computed networks from cache.")

    rows_new = []
    skipped  = 0
    errors   = 0

    for folder, dataset_label in SCAN_FOLDERS.items():
        if not os.path.isdir(folder):
            print(f"[WARN] folder not found: {folder}")
            continue

        txt_files = sorted(f for f in os.listdir(folder) if f.endswith('.txt'))
        print(f"\n=== {dataset_label}: {len(txt_files)} networks ===")

        for fname in txt_files:
            if fname in done:
                continue

            txt_path = os.path.join(folder, fname)

            try:
                RN = read_txt(txt_path)
                n_sp = len(RN.species())
                n_rx = len(RN.reactions())

                if n_rx > MAX_REACTIONS:
                    print(f"  SKIP {fname}: {n_rx} reactions > {MAX_REACTIONS}")
                    skipped += 1
                    continue

                # Compute ERCs (exclude E_∅: the empty-closure ERC present in all)
                t0 = time.time()
                ercs, _cached = load_ercs(txt_path, RN, ERC)
                ercs = [e for e in ercs if len(e.get_closure_names(RN)) > 0]
                t_ercs = time.time() - t0

                if len(ercs) > MAX_ERCS:
                    print(f"  SKIP {fname}: {len(ercs)} ERCs > {MAX_ERCS}")
                    skipped += 1
                    continue
                if len(ercs) < MIN_ERCS:
                    print(f"  SKIP {fname}: only {len(ercs)} ERCs < {MIN_ERCS}")
                    skipped += 1
                    continue

                # Build hierarchy
                t1 = time.time()
                hierarchy = ERC_Hierarchy(RN, ercs)
                G = hierarchy.graph
                metrics = forest_metrics(ercs, G)
                t_forest = time.time() - t1

                row = {
                    'file':             fname,
                    'dataset':          dataset_label,
                    'n_species':        n_sp,
                    'n_reactions':      n_rx,
                    'erc_reaction_ratio': round(len(ercs) / n_rx, 4) if n_rx > 0 else float('nan'),
                    'time_ercs':        round(t_ercs, 3),
                    'time_forest':      round(t_forest, 3),
                }
                row.update(metrics)
                rows_new.append(row)
                print(f"  OK  {fname:45s} ercs={len(ercs):3d}  trees={metrics['n_trees']:3d}"
                      f"  depth={metrics['max_depth']:2d}  ratio={row['erc_reaction_ratio']:.3f}")

            except Exception as e:
                print(f"  ERR {fname}: {e}")
                errors += 1

    # Merge and save
    all_rows = rows_existing + rows_new
    df = pd.DataFrame(all_rows)
    df.to_csv(CSV_FILE, index=False)

    print(f"\nDone.  New: {len(rows_new)}  Skipped: {skipped}  Errors: {errors}")
    print(f"Total rows in CSV: {len(df)}")
    print(f"Saved -> {CSV_FILE}")

    if not df.empty:
        print("\nOverall forest stats summary:")
        cols = ['n_ercs', 'erc_reaction_ratio', 'n_trees', 'max_depth',
                'mean_depth', 'n_roots', 'n_leaves', 'mean_out_degree',
                'frac_isolated']
        print(df[cols].describe().round(3).to_string())


if __name__ == '__main__':
    main()
