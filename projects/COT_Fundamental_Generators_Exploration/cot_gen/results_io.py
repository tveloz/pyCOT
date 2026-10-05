"""
results_io.py — Persistent CSV results store for the COT pipeline.

One row per network run.  Dynamic column groups grow automatically when a
new run introduces a larger length / order than any previously stored row:

  elem_len{N} — elementary SOs (SO0) whose erc_set has exactly N members
  so_o{K}     — SOs of order K (order 0 = elementary, order 1 = first extension layer, …)
  leaves_len{N} — DFS dead-ends (req≠0, no extension) with erc_set size N

All other columns are fixed.

Public API
----------
make_row(network_id, rn, ercs_list, syn_fund, comp, elem, so_hier,
         timing_ms, status="ok") -> dict

update_results_csv(csv_path, row) -> None
    Read the existing CSV (if any), replace or append the row (matched by
    "network"), expand dynamic columns if needed, rewrite atomically.
"""
from __future__ import annotations
import csv
import os

# ── Fixed column groups ───────────────────────────────────────────────────────

_FRONT = [
    "network", "species", "reactions", "ercs", "p_ercs",
    "fund_syn", "fund_comp",
]

_TOTALS = [
    "n_elem",          # total elementary SOs (SO0)
    "n_so_higher",     # total higher-order SOs (SOi, i>=1)
    "n_sos",           # all semi-organizing sets found (every order)
    "n_leaves",        # DFS dead-ends (req≠0, no extension possible)
]

_PERF = [
    "dfs_explored",            # total DFS states processed
    "dfs_branching",           # explored − SSMs − dead-ends
    "dfs_canonical_pruned",    # extensions skipped by canonical ordering
]

_TIMING = [
    "t_load_ms", "t_erc_ms", "t_hier_ms", "t_syn_ms",
    "t_comp_ms", "t_gen_ms", "t_elem_ms", "t_so_ms", "t_total_ms",
]

_END = ["status"]


# ── Dynamic column name helpers ───────────────────────────────────────────────

def _elem_col(n: int)   -> str: return f"elem_len{n}"
def _so_col(k: int)     -> str: return f"so_o{k}"
def _leaves_col(n: int) -> str: return f"leaves_len{n}"


def _scan_maximums(rows: list[dict]) -> tuple[int, int, int]:
    """Return (max_elem_len, max_so_order, max_leaves_len) across all rows."""
    max_elem = max_leaves = 0
    max_so = -1
    for row in rows:
        for col in row:
            if col.startswith("elem_len"):
                try: max_elem = max(max_elem, int(col[8:]))
                except ValueError: pass
            elif col.startswith("so_o"):
                try: max_so = max(max_so, int(col[4:]))
                except ValueError: pass
            elif col.startswith("leaves_len"):
                try: max_leaves = max(max_leaves, int(col[10:]))
                except ValueError: pass
    return max_elem, max_so, max_leaves


def _column_order(max_elem_len: int, max_so_order: int, max_leaves_len: int) -> list[str]:
    cols = list(_FRONT)
    cols += list(_TOTALS)
    cols += [_elem_col(n) for n in range(1, max_elem_len + 1)]
    if max_so_order >= 0:
        cols += [_so_col(k) for k in range(0, max_so_order + 1)]
    cols += list(_PERF)
    cols += [_leaves_col(n) for n in range(1, max_leaves_len + 1)]
    cols += list(_TIMING)
    cols += list(_END)
    return cols


# ── Row builder ───────────────────────────────────────────────────────────────

def make_row(
    network_id: str,
    rn,
    ercs_list,
    syn_fund,
    comp,
    elem,
    so_hier,
    timing_ms: dict,
    status: str = "ok",
) -> dict:
    """
    Build a flat dict row from pipeline outputs.

    Any of rn / ercs_list / syn_fund / comp / elem / so_hier may be None
    (stages that were skipped or failed).

    timing_ms keys: load, erc, hier, syn, comp, gen, elem, so  (values in ms).
    """
    row: dict = {}

    # ── Identity ──────────────────────────────────────────────────────────────
    row["network"]   = network_id
    row["species"]   = rn.n_species   if rn else ""
    row["reactions"] = rn.n_reactions if rn else ""
    if ercs_list is not None:
        row["ercs"]  = len(ercs_list)
        row["p_ercs"] = sum(1 for e in ercs_list if e.is_persistent())
    else:
        row["ercs"] = row["p_ercs"] = ""

    row["fund_syn"]  = len(syn_fund.fundamental) if syn_fund else 0
    row["fund_comp"] = len(comp.fundamental)      if comp     else 0

    # ── Elementary-SO counts ─────────────────────────────────────────────────
    if elem is not None:
        row["n_elem"]   = len(elem.all_elementary_masks)
        row["n_leaves"] = len(elem.leaf_masks)

        # single-ERC elementary SOs always have length 1
        row[_elem_col(1)] = len(elem.single_erc_indices)

        # multi-ERC elementary SOs: from DFS ssm_by_length (keys are erc_set sizes ≥ 2)
        for n, cnt in elem.stats.get("ssm_by_length", {}).items():
            col = _elem_col(int(n))
            row[col] = row.get(col, 0) + cnt

        # DFS performance
        n_expl = elem.stats.get("states_explored", 0)
        n_ssm  = elem.stats.get("ssm_found",       0)
        n_dead = elem.stats.get("leaves_found",     0)
        row["dfs_explored"]         = n_expl
        row["dfs_branching"]        = max(0, n_expl - n_ssm - n_dead)
        row["dfs_canonical_pruned"] = elem.stats.get("canonical_pruned", 0)

        # Dead-ends by length
        for n, cnt in elem.stats.get("dead_by_length", {}).items():
            row[_leaves_col(int(n))] = cnt
    else:
        row["n_elem"] = row["n_leaves"] = 0
        row["dfs_explored"] = row["dfs_branching"] = row["dfs_canonical_pruned"] = 0

    # ── Higher-order SO counts ───────────────────────────────────────────────
    if so_hier is not None:
        row["n_so_higher"] = so_hier.total_so()
        row["n_sos"]        = len(so_hier.all_so_masks)
        for k, masks in so_hier.so_by_order.items():
            row[_so_col(int(k))] = len(masks)
    else:
        # No higher-order search: the only SOs we know are the elementary
        # ones themselves (order 0).
        n_ep = len(elem.all_elementary_masks) if elem else 0
        row["n_so_higher"] = n_ep
        row["n_sos"]         = n_ep
        if n_ep:
            row[_so_col(0)] = n_ep

    # ── Timing ────────────────────────────────────────────────────────────────
    key_map = {
        "load": "t_load_ms", "erc":  "t_erc_ms",  "hier": "t_hier_ms",
        "syn":  "t_syn_ms",  "comp": "t_comp_ms",  "gen":  "t_gen_ms",
        "elem": "t_elem_ms", "so":   "t_so_ms",
    }
    for k, col in key_map.items():
        row[col] = round(timing_ms.get(k, 0.0), 1)
    row["t_total_ms"] = round(sum(timing_ms.values()), 1)

    row["status"] = status
    return row


# ── CSV writer ────────────────────────────────────────────────────────────────

def update_results_csv(csv_path: str, new_row: dict) -> None:
    """
    Read the existing CSV (if any), replace or append new_row (matched by
    'network'), expand dynamic columns if the new row introduces larger
    lengths or orders, and rewrite the file.

    Calling this after every network means partial batch runs are recoverable.
    """
    # ── Read existing rows ────────────────────────────────────────────────────
    existing: list[dict] = []
    if os.path.exists(csv_path) and os.path.getsize(csv_path) > 0:
        with open(csv_path, newline="", encoding="utf-8") as fh:
            existing = list(csv.DictReader(fh))

    # ── Replace matching row or append ────────────────────────────────────────
    replaced = False
    for i, row in enumerate(existing):
        if row.get("network") == new_row.get("network"):
            existing[i] = new_row
            replaced = True
            break
    if not replaced:
        existing.append(new_row)

    # ── Discover max dynamic indices across ALL rows (including the new one) ──
    max_elem, max_so, max_leaves = _scan_maximums(existing)

    # ── Build full ordered column list ───────────────────────────────────────
    cols = _column_order(max_elem, max(max_so, 0), max_leaves)

    # ── Write ─────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        _str_cols = set(_FRONT) | set(_END)
        for row in existing:
            out = {}
            for col in cols:
                if col in row:
                    out[col] = row[col]
                elif col in _str_cols:
                    out[col] = ""
                else:
                    out[col] = 0
            writer.writerow(out)
