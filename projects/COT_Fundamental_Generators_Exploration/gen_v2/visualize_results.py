"""
visualize_results.py -- plots for compare_old_new.py's results
(outputs/gen_v2_comparison/old_vs_new_*.csv).

Reads EVERY old_vs_new_*.csv found under OUT_DIR (not just one), so sweeps
run at different MIN_REACTIONS/MAX_REACTIONS windows (e.g. a future
800-2000 run targeting larger networks) all show up on the same plot. If
a (network, engine) pair appears in more than one file -- a network
re-run at a wider window, or re-run after a checkpoint resume -- the
later file (by filename sort) wins, on the assumption that a re-run
reflects more completed work, not less.

Produces two panels in outputs/gen_v2_comparison/comparison.png:
  1. old vs new total-SO count per network (symlog, since counts span
     0 to several thousand) -- the headline correctness/coverage check.
  2. wall-clock time vs network size per engine, with the per-run time
     budget marked -- shows where the scaling wall actually is.
Also prints a short agreement/regression tally to stdout.

Run: python -m gen_v2.visualize_results
"""
from __future__ import annotations

import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
OUT_DIR = os.path.join(_proj, "outputs", "gen_v2_comparison")
TIME_BUDGET_S = 300  # reference line only -- see compare_old_new.py's TIMEOUT_S

NUMERIC_COLS = ("n_reactions", "n_species", "n_ercs", "wall_time_s",
                "n_so_total", "n_elementary", "max_order_reached", "states_explored")


def load_all(csv_glob: str | None = None) -> pd.DataFrame:
    paths = sorted(glob.glob(csv_glob or os.path.join(OUT_DIR, "old_vs_new_*.csv")))
    if not paths:
        raise SystemExit(f"No old_vs_new_*.csv files found under {OUT_DIR} -- "
                          f"run gen_v2.compare_old_new first.")
    frames = []
    for p in paths:
        df = pd.read_csv(p)
        df["_source"] = os.path.basename(p)
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    df = df.drop_duplicates(subset=["network", "engine"], keep="last").reset_index(drop=True)
    for col in NUMERIC_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    print(f"Loaded {len(df)} rows from {len(paths)} file(s): {[os.path.basename(p) for p in paths]}")
    return df


def plot_so_counts(df: pd.DataFrame, ax) -> None:
    piv = df.pivot_table(index="network", columns="engine", values="n_so_total", aggfunc="first")
    status = df.pivot_table(index="network", columns="engine", values="status", aggfunc="first")
    sizes = df.groupby("network")["n_reactions"].max()
    order = sizes.sort_values(na_position="last").index
    piv, status = piv.reindex(order), status.reindex(order)

    x = np.arange(len(piv))
    width = 0.38
    colors = {"old": "#4C72B0", "new": "#DD8452"}
    for i, engine in enumerate(("old", "new")):
        if engine not in piv.columns:
            continue
        vals = piv[engine].fillna(0).values
        ax.bar(x + (i - 0.5) * width, vals, width, label=engine, color=colors[engine])
        for j, st in enumerate(status.get(engine, pd.Series(index=piv.index)).values):
            if st not in ("ok", None) and not (isinstance(st, float) and np.isnan(st)):
                ax.text(x[j] + (i - 0.5) * width, 1, str(st), rotation=90,
                         ha="center", va="bottom", fontsize=7, color="dimgray")
    ax.set_yscale("symlog")
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("semi-organizations found (symlog)")
    ax.set_title("Old (production) vs new (gen_v2) engine: total SOs per network")
    ax.legend()


def plot_scaling(df: pd.DataFrame, ax) -> None:
    markers = {"old": "o", "new": "s"}
    colors = {"old": "#4C72B0", "new": "#DD8452"}
    for engine in ("old", "new"):
        sub = df[(df["engine"] == engine) & df["n_reactions"].notna() & df["wall_time_s"].notna()]
        ok = sub[sub["status"] == "ok"]
        other = sub[sub["status"] != "ok"]
        ax.scatter(ok["n_reactions"], ok["wall_time_s"], color=colors[engine],
                   marker=markers[engine], label=f"{engine} (complete)", alpha=0.85)
        if len(other):
            ax.scatter(other["n_reactions"], other["wall_time_s"], color=colors[engine],
                       marker="x", s=80, label=f"{engine} (partial/timeout)")
    ax.axhline(TIME_BUDGET_S, color="red", linestyle="--", linewidth=1,
               label=f"{TIME_BUDGET_S}s time budget")
    ax.set_xlabel("reactions")
    ax.set_ylabel("wall time (s, log)")
    ax.set_yscale("log")
    ax.set_title("Wall-clock scaling vs network size")
    ax.legend(fontsize=8, loc="upper left")


def print_summary(df: pd.DataFrame) -> None:
    print("=" * 70)
    nets = sorted(df["network"].unique())
    agree = new_superset = old_superset = incomplete = 0
    regressions = []
    for net in nets:
        rows = {r.engine: r for r in df[df.network == net].itertuples()}
        old, new = rows.get("old"), rows.get("new")
        if old is None or new is None or old.status != "ok" or new.status != "ok":
            incomplete += 1
            continue
        if old.n_so_total == new.n_so_total:
            agree += 1
        elif new.n_so_total > old.n_so_total:
            new_superset += 1
        else:
            old_superset += 1
            regressions.append(net)
    print(f"{len(nets)} networks total  |  "
          f"agree={agree}  new-superset={new_superset}  "
          f"old-superset={old_superset}  incomplete={incomplete}")
    if regressions:
        print(f"*** POSSIBLE REGRESSION -- old found MORE than new on: {regressions} ***")
    print("=" * 70)


def main():
    df = load_all()
    print_summary(df)
    n_net = df["network"].nunique()
    fig, axes = plt.subplots(2, 1, figsize=(max(10, 0.45 * n_net), 11))
    plot_so_counts(df, axes[0])
    plot_scaling(df, axes[1])
    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, "comparison.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
