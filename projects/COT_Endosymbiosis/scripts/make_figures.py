"""
make_figures.py -- summary figures for the endosymbiosis project.

Fig 1: before/after EPM count + max/mean size, grouped by case (toy,
mitochondrial-type, chloroplast-type dark, chloroplast-type light), each
showing host-alone / symbiont-alone / merged as three bars, with hybrid
fraction annotated.

Run AFTER all endosymbiosis_analysis.py runs are complete; numbers are
hardcoded here from the logged results (see outputs/*.log for the raw
source of each number) rather than re-running the pipeline, since this is
purely a plotting step -- keep it that way, don't silently recompute with
different parameters than what's in the logs.
"""
from __future__ import annotations

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_fig_dir = os.path.join(_here, "..", "figures")
os.makedirs(_fig_dir, exist_ok=True)

# name -> {host: (n_epms, mean, min, max), symb: (...), merged: (n_epms, mean, min, max, hybrid_frac)}
# Numbers transcribed from outputs/*.log (raw pipeline output) -- see
# PROGRESS.md for the full provenance of each case.
CASES = {
    "toy model": {
        "host": (1, 8.0, 8, 8),
        "symb": (1, 12.0, 12, 12),
        "merged": (1, 22.0, 22, 22, 1.0),
    },
    "mitochondrial-type\n(E. coli core,\nself-derived)": {
        "host": (15, 15.1, 13, 17),
        "symb": (15, 15.1, 13, 17),
        "merged": (27, 28.2, 26, 30, 1.0),
    },
    "chloroplast-type\n(yeast + Synecho.,\ndark)": {
        "host": (152, 42.9, 40, 50),
        "symb": (63, 24.0, 21, 33),
        "merged": (221, 66.0, 63, 75, 1.0),
    },
    "chloroplast-type\n(yeast + Synecho.,\n+light)": {
        "host": (152, 42.9, 40, 50),
        "symb": (63, 27.0, 24, 36),
        "merged": (221, 69.0, 66, 78, 1.0),
    },
}


def plot_cases(cases: dict, out_name: str = "fig_complexification_summary"):
    names = list(cases.keys())
    x = np.arange(len(names))
    w = 0.25

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))

    for i, key in enumerate(["host", "symb", "merged"]):
        counts = [cases[n][key][0] for n in names]
        offset = (i - 1) * w
        color = {"host": "#3498db", "symb": "#9b59b6", "merged": "#27ae60"}[key]
        label = {"host": "host alone", "symb": "symbiont alone", "merged": "merged"}[key]
        ax1.bar(x + offset, counts, width=w * 0.9, color=color, label=label)

    ax1.set_yscale("log")
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, fontsize=7.5)
    ax1.set_ylabel("number of EPMs (log scale)")
    ax1.set_title("(a) EPM count")
    ax1.legend(fontsize=7, frameon=False)
    ax1.spines[['top', 'right']].set_visible(False)

    for i, key in enumerate(["host", "symb", "merged"]):
        means = [cases[n][key][1] for n in names]
        mins = [cases[n][key][2] for n in names]
        maxs = [cases[n][key][3] for n in names]
        offset = (i - 1) * w
        color = {"host": "#3498db", "symb": "#9b59b6", "merged": "#27ae60"}[key]
        for xi, (m, lo, hi) in zip(x + offset, zip(means, mins, maxs)):
            ax2.plot([xi, xi], [lo, hi], color=color, linewidth=1.6, solid_capstyle='round')
            ax2.scatter([xi], [m], color=color, s=26, zorder=3,
                        edgecolor="#2d3436", linewidth=0.4)

    for i, n in enumerate(names):
        merged = cases[n]["merged"]
        frac = merged[4] if len(merged) > 4 else None
        if frac is not None:
            ax2.text(x[i] + w, merged[3] + 2, f"{frac*100:.0f}% hybrid",
                      fontsize=6.5, ha="center", color="#27ae60")

    ax2.set_xticks(x)
    ax2.set_xticklabels(names, fontsize=7.5)
    ax2.set_ylabel("EPM size (species); dot=mean, line=min-max")
    ax2.set_title("(b) EPM size, and hybrid fraction of merged EPMs")
    ax2.spines[['top', 'right']].set_visible(False)

    fig.tight_layout()
    out_pdf = os.path.join(_fig_dir, out_name + ".pdf")
    out_png = os.path.join(_fig_dir, out_name + ".png")
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    print(f"wrote {out_pdf} and {out_png}")


if __name__ == "__main__":
    plot_cases(CASES)
