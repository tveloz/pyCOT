"""
make_fig_organisms.py -- two-panel figure comparing EPM count and size range
across five organisms (E. coli iML1515 as reference + four new BiGG
reconstructions) under matched core environments plus each organism's own
full native medium.

Reads:
  outputs/ecoli_epm_comparison/results.csv          (for iML1515 rows)
  outputs/cross_organism_epm_comparison/results.csv (for the 4 new organisms)
"""
from __future__ import annotations

import os
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
ECOLI_CSV = os.path.join(_repo_root, 'projects', 'COT_Fundamental_Generators_Exploration',
                          'outputs', 'ecoli_epm_comparison', 'results.csv')
CROSS_CSV = os.path.join(_repo_root, 'projects', 'COT_Fundamental_Generators_Exploration',
                          'outputs', 'cross_organism_epm_comparison', 'results.csv')

ORGANISM_LABEL = {
    'iML1515': 'E. coli\n(iML1515)',
    'iJN678': 'Synechocystis\n(iJN678)',
    'iYO844': 'B. subtilis\n(iYO844)',
    'iND750': 'S. cerevisiae\n(iND750)',
    'iAF987': 'G. metallireducens\n(iAF987)',
}
ORGANISMS = ['iML1515', 'iJN678', 'iYO844', 'iND750', 'iAF987']
SCENARIOS = ['aerobic_core', 'anaerobic_core', 'carbon_starvation_core', 'full_native']
SCEN_LABEL = {'aerobic_core': 'aerobic', 'anaerobic_core': 'anaerobic',
              'carbon_starvation_core': 'C-starved', 'full_native': 'full native'}
SCEN_COLOR = {'aerobic_core': '#3498db', 'anaerobic_core': '#9b59b6',
              'carbon_starvation_core': '#e74c3c', 'full_native': '#27ae60'}

data = {}
for csv_path, models in [(ECOLI_CSV, {'iML1515'}), (CROSS_CSV, {'iJN678', 'iYO844', 'iND750', 'iAF987'})]:
    with open(csv_path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            if row['model'] in models:
                data[(row['model'], row['scenario'])] = row

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.6))

n_org = len(ORGANISMS)
n_scen = len(SCENARIOS)
bar_w = 0.19
x = np.arange(n_org)

for si, scen in enumerate(SCENARIOS):
    counts = []
    for org in ORGANISMS:
        row = data.get((org, scen))
        counts.append(int(row['n_epms']) if row else np.nan)
    offset = (si - (n_scen - 1) / 2.0) * bar_w
    axA.bar(x + offset, counts, width=bar_w * 0.92, color=SCEN_COLOR[scen],
            label=SCEN_LABEL[scen])

axA.set_xticks(x)
axA.set_xticklabels([ORGANISM_LABEL[o] for o in ORGANISMS], fontsize=7.5)
axA.set_ylabel("number of EPMs")
axA.set_title("(a) EPM count", fontsize=10)
axA.set_ylim(0, 260)
axA.legend(fontsize=7, frameon=False, ncol=4, loc="upper center",
           bbox_to_anchor=(0.5, 1.14), columnspacing=1.0, handletextpad=0.4)
axA.spines[['top', 'right']].set_visible(False)

for si, scen in enumerate(SCENARIOS):
    offset = (si - (n_scen - 1) / 2.0) * bar_w
    for oi, org in enumerate(ORGANISMS):
        row = data.get((org, scen))
        if row is None:
            continue
        lo, mean, hi = int(row['epm_size_min']), float(row['epm_size_mean']), int(row['epm_size_max'])
        xpos = oi + offset
        axB.plot([xpos, xpos], [lo, hi], color=SCEN_COLOR[scen], linewidth=1.6, solid_capstyle='round')
        axB.scatter([xpos], [mean], color=SCEN_COLOR[scen], edgecolor="#2d3436",
                    s=16, zorder=3, linewidth=0.5)

axB.set_xticks(x)
axB.set_xticklabels([ORGANISM_LABEL[o] for o in ORGANISMS], fontsize=7.5)
axB.set_ylabel("EPM size (species), min–max, dot = mean")
axB.set_title("(b) EPM size range", fontsize=10)
axB.spines[['top', 'right']].set_visible(False)

fig.suptitle("")
fig.tight_layout()

out_pdf = os.path.join(_here, "fig4_cross_organism_comparison.pdf")
out_png = os.path.join(_here, "fig4_cross_organism_comparison.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"wrote {out_pdf} and {out_png}")

# Print a quick text summary for writing the Results paragraph accurately.
print("\nSummary table:")
for org in ORGANISMS:
    for scen in SCENARIOS:
        row = data.get((org, scen))
        if row:
            print(f"  {org:8s} {scen:24s} n_epms={row['n_epms']:>4s}  "
                  f"size={row['epm_size_min']}-{row['epm_size_max']} (mean {row['epm_size_mean']})  "
                  f"n_ercs={row['n_ercs']}  food_n={len(row['food_species'].split(','))}")
