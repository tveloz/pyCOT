"""
results_io.py — Persistent CSV results store for the COT pipeline.

One row per network run.  Dynamic column groups grow automatically when a
new run introduces a larger length / order than any previously stored row:

  epms_len{N}   — EPMs whose erc_set has exactly N members
  espm_o{K}     — SOs of order K (order 0 = EPMs, order 1 = first extension layer, …)
  leaves_len{N} — DFS dead-ends (req≠0, no extension) with erc_set size N

All other columns are fixed.

Public API
----------
make_row(network_id, rn, ercs_list, syn_fund, comp, epm, espm,
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
    "n_epms",          # total EPMs
    "n_espms",         # total ESPMs (all orders)
    "n_sos",           # all semi-organizing sets found
    "n_leaves",        # DFS dead-ends (req≠0, no extension possible)
]

_PERF = [
    "dfs_explored",            # total DFS states processed
    "dfs_branching",           # explored − SSMs − dead-ends
    "dfs_canonical_pruned",    # extensions skipped by canonical ordering
]

_TIMING = [
    "t_load_ms", "t_erc_ms", "t_hier_ms", "t_syn_ms",
    "t_comp_ms", "t_gen_ms", "t_epm_ms", "t_espm_ms", "t_total_ms",
]

_END = ["status"]


# ── Dynamic column name helpers ───────────────────────────────────────────────

def _epms_col(n: int)   -> str: return f"epms_len{n}"
def _espm_col(k: int)   -> str: return f"espm_o{k}"
def _leaves_col(n: int) -> str: return f"leaves_len{n}"


def _scan_maximums(rows: list[dict]) -> tuple[int, int, int]:
    """Return (max_epm_len, max_espm_order, max_leaves_len) across all rows."""
    max_epm = max_leaves = 0
    max_espm = -1
    for row in rows:
        for col in row:
            if col.startswith("epms_len"):
                try: max_epm = max(max_epm, int(col[8:]))
                except ValueError: pass
            elif col.startswith("espm_o"):
                try: max_espm = max(max_espm, int(col[6:]))
                except ValueError: pass
            elif col.startswith("leaves_len"):
                try: max_leaves = max(max_leaves, int(col[10:]))
                except ValueError: pass
    return max_epm, max_espm, max_leaves


def _column_order(max_epm_len: int, max_espm_order: int, max_leaves_len: int) -> list[str]:
    cols = list(_FRONT)
    cols += list(_TOTALS)
    cols += [_epms_col(n)  for n in range(1, max_epm_len   + 1)]
    if max_espm_order >= 0:
        cols += [_espm_col(k) for k in range(0, max_espm_order + 1)]
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
    epm,
    espm,
    timing_ms: dict,
    status: str = "ok",
) -> dict:
    """
    Build a flat dict row from pipeline outputs.

    Any of rn / ercs_list / syn_fund / comp / epm / espm may be None
    (stages that were skipped or failed).

    timing_ms keys: load, erc, hier, syn, comp, gen, epm, espm  (values in ms).
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

    # ── EPM counts ────────────────────────────────────────────────────────────
    if epm is not None:
        row["n_epms"]   = len(epm.all_epm_masks)
        row["n_leaves"] = len(epm.leaf_masks)

        # single-ERC EPMs always have length 1
        row[_epms_col(1)] = len(epm.single_epm_indices)

        # multi-ERC EPMs: from DFS ssm_by_length (keys are erc_set sizes ≥ 2)
        for n, cnt in epm.stats.get("ssm_by_length", {}).items():
            col = _epms_col(int(n))
            row[col] = row.get(col, 0) + cnt

        # DFS performance
        n_expl = epm.stats.get("states_explored", 0)
        n_ssm  = epm.stats.get("ssm_found",       0)
        n_dead = epm.stats.get("leaves_found",     0)
        row["dfs_explored"]         = n_expl
        row["dfs_branching"]        = max(0, n_expl - n_ssm - n_dead)
        row["dfs_canonical_pruned"] = epm.stats.get("canonical_pruned", 0)

        # Dead-ends by length
        for n, cnt in epm.stats.get("dead_by_length", {}).items():
            row[_leaves_col(int(n))] = cnt
    else:
        row["n_epms"] = row["n_leaves"] = 0
        row["dfs_explored"] = row["dfs_branching"] = row["dfs_canonical_pruned"] = 0

    # ── ESPM / SO counts ──────────────────────────────────────────────────────
    if espm is not None:
        row["n_espms"] = espm.total_espm()
        row["n_sos"]   = len(espm.all_so_masks)
        for k, masks in espm.espm_by_order.items():
            row[_espm_col(int(k))] = len(masks)
    else:
        # No ESPM computed: the only SOs we know are the EPMs themselves (order 0)
        n_ep = len(epm.all_epm_masks) if epm else 0
        row["n_espms"] = n_ep
        row["n_sos"]   = n_ep
        if n_ep:
            row[_espm_col(0)] = n_ep

    # ── Timing ────────────────────────────────────────────────────────────────
    key_map = {
        "load": "t_load_ms", "erc":  "t_erc_ms",  "hier": "t_hier_ms",
        "syn":  "t_syn_ms",  "comp": "t_comp_ms",  "gen":  "t_gen_ms",
        "epm":  "t_epm_ms",  "espm": "t_espm_ms",
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
    max_epm, max_espm, max_leaves = _scan_maximums(existing)

    # ── Build full ordered column list ───────────────────────────────────────
    cols = _column_order(max_epm, max(max_espm, 0), max_leaves)

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
