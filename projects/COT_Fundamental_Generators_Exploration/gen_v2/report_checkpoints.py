"""
report_checkpoints.py -- human-readable summary of every gen_v2
checkpoint under outputs/gen_v2_checkpoints/ (written by run_large.py
and/or compare_old_new.py -- they share the same directory and format).

This is the answer to "I ran run_large.py, where are my results": that
script only prints to the terminal (gone once it scrolls away/closes)
and writes raw .pkl state (not human-readable on its own). This script
reads every .pkl, reports what each one found so far, and saves it
somewhere persistent -- a stdout table, a CSV, and a bar-chart PNG.

Run: python -m gen_v2.report_checkpoints
Writes:
  outputs/gen_v2_checkpoints/summary.csv
  outputs/gen_v2_checkpoints/summary.png
"""
from __future__ import annotations

import csv
import glob
import os
import re
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

CKPT_DIR = os.path.join(_proj, "outputs", "gen_v2_checkpoints")
OUT_CSV = os.path.join(CKPT_DIR, "summary.csv")
OUT_PNG = os.path.join(CKPT_DIR, "summary.png")

_LABEL_RE = re.compile(r"\((\d+)\s*rxn,\s*(\d+)\s*ERCs\)")

FIELDNAMES = ["network", "n_reactions", "n_ercs", "status",
              "n_so_total", "n_elementary", "max_order_reached",
              "states_resolved", "checkpoint_size_mb", "last_modified"]


def main():
    from gen_v2.engine import load_checkpoint, build_result

    paths = sorted(p for p in glob.glob(os.path.join(CKPT_DIR, "*.pkl")) if not p.endswith(".tmp"))
    orphan_tmps = sorted(glob.glob(os.path.join(CKPT_DIR, "*.pkl.tmp")))
    if not paths:
        print(f"No checkpoints found under {CKPT_DIR} -- nothing to report yet.")
        return

    rows = []
    print(f"{'network':28s} {'status':10s} {'SOs':>9s} {'elem':>6s} {'order':>6s} "
          f"{'resolved':>11s} {'size(MB)':>9s}  reactions/ERCs")
    print("-" * 110)
    for p in paths:
        name = os.path.splitext(os.path.basename(p))[0]
        try:
            ck = load_checkpoint(p)
        except Exception as exc:
            print(f"{name:28s}  FAILED TO LOAD ({exc}) -- possibly mid-write when the process stopped")
            continue
        result = build_result(ck.discovered_sos, ck.stats, ck.parent_of, complete=ck.complete)
        label = getattr(ck, "network_label", "") or name
        m = _LABEL_RE.search(label)
        n_reactions, n_ercs = (m.group(1), m.group(2)) if m else ("?", "?")
        status = "complete" if ck.complete else "partial"
        max_order = max(result.so_by_order.keys(), default=0)
        size_mb = os.path.getsize(p) / (1024 * 1024)
        mtime = os.path.getmtime(p)
        print(f"{name:28s} {status:10s} {len(result.all_so_masks):>9d} "
              f"{len(result.elementary_masks):>6d} {max_order:>6d} "
              f"{result.stats.get('states_resolved', 0):>11d} {size_mb:>9.1f}  "
              f"{n_reactions} rxn, {n_ercs} ERCs")
        rows.append({
            "network": name, "n_reactions": n_reactions, "n_ercs": n_ercs, "status": status,
            "n_so_total": len(result.all_so_masks), "n_elementary": len(result.elementary_masks),
            "max_order_reached": max_order, "states_resolved": result.stats.get("states_resolved", 0),
            "checkpoint_size_mb": round(size_mb, 1),
            "last_modified": __import__("datetime").datetime.fromtimestamp(mtime).isoformat(timespec="seconds"),
        })

    n_complete = sum(1 for r in rows if r["status"] == "complete")
    print("-" * 110)
    print(f"{len(rows)} checkpoints -- {n_complete} complete, {len(rows) - n_complete} partial (re-run "
          f"gen_v2.run_large on the same names to continue those)")
    if orphan_tmps:
        print(f"\n{len(orphan_tmps)} orphaned .pkl.tmp file(s) (a save that didn't finish replacing the "
              f"real file -- the .pkl itself is still the last GOOD save, these are safe to delete):")
        for t in orphan_tmps:
            print(f"    {t}")

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {OUT_CSV}")

    try:
        _plot(rows)
        print(f"Wrote {OUT_PNG}")
    except ImportError:
        print("(matplotlib not available -- skipped the plot, CSV/table above still complete)")


def _plot(rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    def sort_key(r):
        try:
            return int(r["n_reactions"])
        except (ValueError, TypeError):
            return 10 ** 9
    rows = sorted(rows, key=sort_key)

    names = [r["network"] for r in rows]
    vals = [r["n_so_total"] for r in rows]
    colors = ["#55A868" if r["status"] == "complete" else "#C44E52" for r in rows]

    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(rows)), 6))
    x = np.arange(len(rows))
    ax.bar(x, vals, color=colors)
    for i, r in enumerate(rows):
        if r["status"] != "complete":
            ax.text(x[i], vals[i], "partial", rotation=90, ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=60, ha="right", fontsize=8)
    ax.set_yscale("symlog")
    ax.set_ylabel("semi-organizations found so far (symlog)")
    ax.set_title("gen_v2 checkpoint status (green=complete, red=partial)")
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150)


if __name__ == "__main__":
    main()
