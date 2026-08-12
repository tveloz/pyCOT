"""
make_fig_illustrative.py -- two-panel static figure explaining the ERC /
fundamental synergy / fundamental complementarity vocabulary used in the
paper, on a small, REAL (not invented) fragment of E. coli central carbon
metabolism: PGI, ADK1 (reverse), PFK, FBA, TPI, GAPD, taken verbatim from
e_coli_core's own BiGG reaction list (illustrative_glycolysis_candidate.txt).

Panel A: the bipartite reaction network (species circles, reaction boxes).
Panel B: the resulting ERC hierarchy -- containment (gray), the one
fundamental synergy (orange, E1+E2 -> E6), and the fundamental
complementarity chain (blue dashed, labeled by the species supplied).

Both panels are computed from the real pyCOT pipeline output, not hand-
drawn -- see the printed diagnostic this was built from for the exact
E0/ERC/relation numbers reproduced here.
"""
from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
import networkx as nx

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations import (
    build_rndata, compute_ercs, build_hierarchy,
    compute_synergies_basis_first, compute_complementarities,
)

NET_PATH = os.path.join(_here, 'illustrative_glycolysis_candidate.txt')

# ---------------------------------------------------------------------------
# Recompute (small, ~instant) -- keeps the figure exactly in sync with the
# network file rather than hardcoding numbers that could drift.
# ---------------------------------------------------------------------------
rn = read_txt(NET_PATH, exact_names=True)
rn_data = build_rndata(rn, network_id='illustrative')
ercs = compute_ercs(rn_data, verify=False)
hier = build_hierarchy(ercs)
syn = compute_synergies_basis_first(ercs, hier)
comp = compute_complementarities(ercs, hier, syn)
names = rn_data.species_names

print(f"[make_fig_illustrative] {len(ercs)} ERCs, {len(syn.fundamental)} fundamental "
      f"synerg{'y' if len(syn.fundamental)==1 else 'ies'}, "
      f"{len(comp.fundamental)} fundamental complementarities")

# ---------------------------------------------------------------------------
# Panel A -- the bipartite reaction network
# ---------------------------------------------------------------------------
REACTIONS = [
    ("PGI",       ["g6p_c"],                 ["f6p_c"]),
    ("ADK1_rev",  ["adp_c"],                 ["amp_c", "atp_c"]),
    ("PFK",       ["atp_c", "f6p_c"],        ["adp_c", "fdp_c", "h_c"]),
    ("FBA",       ["fdp_c"],                 ["dhap_c", "g3p_c"]),
    ("TPI",       ["dhap_c"],                ["g3p_c"]),
    ("GAPD",      ["g3p_c", "nad_c", "pi_c"],["13dpg_c", "h_c", "nadh_c"]),
]
FOOD = {"nad_c", "pi_c"}

G = nx.DiGraph()
for sp in names:
    G.add_node(f"sp:{sp}", kind="species", food=(sp in FOOD))
for rxn, subs, prods in REACTIONS:
    G.add_node(f"rx:{rxn}", kind="reaction")
    for s in subs:
        G.add_edge(f"sp:{s}", f"rx:{rxn}")
    for p in prods:
        G.add_edge(f"rx:{rxn}", f"sp:{p}")

# Manual layered layout (left->right along the real pathway order), with an
# explicit y-slot per node (not alphabetical) so the two upstream branches
# (PGI/glycolysis on top, ADK1/adenylate-kinase on bottom) stay visually
# parallel and edges don't needlessly cross.
pos_A = {
    "sp:g6p_c":    (0.0,  1.6),
    "sp:adp_c":    (0.0, -1.6),
    "sp:nad_c":    (0.0,  0.55),
    "sp:pi_c":     (0.0, -0.55),

    "rx:PGI":      (1.35,  1.6),
    "rx:ADK1_rev": (1.35, -1.6),

    "sp:f6p_c":    (2.7,  1.6),
    "sp:atp_c":    (2.7, -1.1),
    "sp:amp_c":    (2.7, -2.1),

    "rx:PFK":      (4.05, 0.6),

    "sp:fdp_c":    (5.4,  0.9),
    "sp:h_c":      (5.4, -0.35),

    "rx:FBA":      (6.75, 0.9),

    "sp:dhap_c":   (8.1,  1.5),
    "sp:g3p_c":    (8.1,  0.3),

    "rx:TPI":      (9.45, 1.5),

    "rx:GAPD":     (10.8, 0.3),

    "sp:13dpg_c":  (12.15, 0.9),
    "sp:nadh_c":   (12.15, -0.3),
}

fig = plt.figure(figsize=(11.5, 5.4))
axA = fig.add_axes([0.03, 0.08, 0.46, 0.86])
axB = fig.add_axes([0.55, 0.08, 0.43, 0.86])

SPECIES_COLOR = "#dfe6e9"
FOOD_COLOR = "#74b9ff"
REACTION_COLOR = "#ffeaa7"
EDGE_COLOR = "#636e72"

for node, (x, y) in pos_A.items():
    kind = G.nodes[node]["kind"]
    label = node.split(":", 1)[1]
    if kind == "species":
        color = FOOD_COLOR if G.nodes[node].get("food") else SPECIES_COLOR
        axA.add_patch(Circle((x, y), 0.33, facecolor=color, edgecolor="#2d3436",
                              linewidth=1.0, zorder=3))
        axA.text(x, y, label.replace("_c", ""), ha="center", va="center",
                  fontsize=6.6, zorder=4)
    else:
        axA.add_patch(FancyBboxPatch((x - 0.36, y - 0.20), 0.72, 0.40,
                                      boxstyle="round,pad=0.02,rounding_size=0.06",
                                      facecolor=REACTION_COLOR, edgecolor="#2d3436",
                                      linewidth=1.0, zorder=3))
        axA.text(x, y, label, ha="center", va="center", fontsize=6.6,
                  fontweight="bold", zorder=4)

for u, v in G.edges():
    x1, y1 = pos_A[u]
    x2, y2 = pos_A[v]
    arrow = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                             mutation_scale=9, color=EDGE_COLOR, linewidth=0.9,
                             shrinkA=13, shrinkB=13, zorder=2,
                             connectionstyle="arc3,rad=0.08")
    axA.add_patch(arrow)

axA.set_xlim(-0.6, 9.0 * 1.35 + 0.6)
axA.set_ylim(-2.1, 2.1)
axA.axis("off")
axA.set_title("(a) Illustrative reaction network\n(6 real E. coli reactions: PGI, ADK1, PFK, FBA, TPI, GAPD)",
               fontsize=9)
axA.scatter([], [], marker='o', color=FOOD_COLOR, edgecolor="#2d3436", s=90, label="food species")
axA.scatter([], [], marker='o', color=SPECIES_COLOR, edgecolor="#2d3436", s=90, label="species")
axA.scatter([], [], marker='s', color=REACTION_COLOR, edgecolor="#2d3436", s=90, label="reaction")
axA.legend(loc="lower center", bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=7,
           frameon=False, handletextpad=0.3, columnspacing=1.0)

# ---------------------------------------------------------------------------
# Panel B -- the ERC hierarchy: containment + fundamental synergy/complementarity
# ---------------------------------------------------------------------------
n = len(ercs)
levels = [0] * n
children_of = [[] for _ in range(n)]
for i in range(n):
    for p in hier.parents[i]:
        children_of[p].append(i)
order_by_desc = sorted(range(n), key=lambda i: len(hier.descendants[i]))
for i in order_by_desc:
    if children_of[i]:
        levels[i] = 1 + max(levels[c] for c in children_of[i])

by_level: dict[int, list[int]] = {}
for i in range(n):
    by_level.setdefault(levels[i], []).append(i)

pos_B = {}
for lvl, idxs in by_level.items():
    m = len(idxs)
    for k, i in enumerate(sorted(idxs)):
        x = (k - (m - 1) / 2.0) * 2.0
        pos_B[i] = (x, lvl * 1.9)

COLOR_CONTAINMENT = "#7f8c8d"
COLOR_SYNERGY = "#e67e22"
COLOR_COMPLEMENTARITY = "#2980b9"
COLOR_NODE = "#dfe6e9"
COLOR_NODE_E0 = "#74b9ff"

for i in range(n):
    x, y = pos_B[i]
    sp = [names[j] for j in range(rn_data.n_species) if (ercs[i].species_mask >> j) & 1]
    is_e0 = (ercs[i].species_mask | rn_data.E0_mask) == rn_data.E0_mask and ercs[i].is_persistent() and \
            ercs[i].species_mask == (rn_data.E0_mask & ercs[i].species_mask | ercs[i].species_mask)
    color = COLOR_NODE_E0 if ercs[i].species_mask == rn_data.E0_mask else COLOR_NODE
    axB.add_patch(Circle((x, y), 0.34, facecolor=color, edgecolor="#2d3436",
                          linewidth=1.1, zorder=3))
    axB.text(x, y, f"E{i}", ha="center", va="center", fontsize=8, fontweight="bold", zorder=4)
    axB.text(x, y - 0.52, f"({ercs[i].size()} sp.)", ha="center", va="top", fontsize=6, zorder=4)

for i in range(n):
    for p in hier.parents[i]:
        x1, y1 = pos_B[i]
        x2, y2 = pos_B[p]
        arrow = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=10, color=COLOR_CONTAINMENT,
                                 linewidth=1.1, shrinkA=15, shrinkB=15, zorder=1)
        axB.add_patch(arrow)

for st in syn.fundamental:
    xi, yi = pos_B[st.i]
    xj, yj = pos_B[st.j]
    xk, yk = pos_B[st.k]
    mx, my = (xi + xj) / 2.0, (yi + yj) / 2.0 - 0.15
    axB.scatter([mx], [my], marker="D", s=90, color=COLOR_SYNERGY, edgecolor="#2d3436",
                zorder=3)
    axB.text(mx, my, "+", ha="center", va="center", fontsize=8, color="white",
             fontweight="bold", zorder=4)
    for (xa, ya) in [(xi, yi), (xj, yj)]:
        axB.add_patch(FancyArrowPatch((xa, ya), (mx, my), arrowstyle="-", color=COLOR_SYNERGY,
                                       linewidth=1.3, shrinkA=15, shrinkB=6, zorder=2))
    axB.add_patch(FancyArrowPatch((mx, my), (xk, yk), arrowstyle="-|>", color=COLOR_SYNERGY,
                                   linewidth=1.3, mutation_scale=10, shrinkA=6, shrinkB=15, zorder=2))

comp_seen = set()
comp_list = [fc for fc in comp.fundamental if (fc.prod_idx, fc.cons_idx) not in comp_seen
             and not comp_seen.add((fc.prod_idx, fc.cons_idx))]
_rads = [-0.30, 0.30, -0.45, 0.15, -0.15, 0.45]
for idx, fc in enumerate(comp_list):
    x1, y1 = pos_B[fc.prod_idx]
    x2, y2 = pos_B[fc.cons_idx]
    rad = _rads[idx % len(_rads)]
    arrow = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                             mutation_scale=10, color=COLOR_COMPLEMENTARITY, linewidth=1.2,
                             linestyle=(0, (4, 2)), shrinkA=17, shrinkB=17, zorder=1,
                             connectionstyle=f"arc3,rad={rad}")
    axB.add_patch(arrow)
    # Label placed along the curved path (approximate midpoint of the arc),
    # offset perpendicular to the chord so it doesn't sit on the line itself.
    mx, my = (x1 + x2) / 2.0, (y1 + y2) / 2.0
    dx, dy = (x2 - x1), (y2 - y1)
    perp = (-dy, dx)
    plen = max((perp[0] ** 2 + perp[1] ** 2) ** 0.5, 1e-6)
    mx += perp[0] / plen * rad * 1.1
    my += perp[1] / plen * rad * 1.1
    sp_name = names[fc.species].replace("_c", "")
    axB.text(mx, my, sp_name, fontsize=6.2, color=COLOR_COMPLEMENTARITY,
              ha="center", va="center", zorder=5,
              bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.9))

axB.text(pos_B[0][0], pos_B[0][1] - 0.85, "E0 = food closure\n(nad, pi) — a separate\npersistent ERC", fontsize=6, ha="center",
         va="top", color="#2d3436")

xs = [p[0] for p in pos_B.values()]
ys = [p[1] for p in pos_B.values()]
axB.set_xlim(min(xs) - 1.2, max(xs) + 1.2)
axB.set_ylim(min(ys) - 1.7, max(ys) + 0.9)
axB.axis("off")
axB.set_title("(b) Fundamental ERC hierarchy\n(containment, 1 synergy, complementarity chain)", fontsize=9)

axB.add_patch(FancyArrowPatch((0, 0), (0, 0), color=COLOR_CONTAINMENT, linewidth=1.1))
h1 = axB.plot([], [], color=COLOR_CONTAINMENT, linewidth=1.4, label="containment")[0]
h2 = axB.scatter([], [], marker="D", color=COLOR_SYNERGY, edgecolor="#2d3436", s=70,
                  label="fundamental synergy")
h3 = axB.plot([], [], color=COLOR_COMPLEMENTARITY, linewidth=1.4, linestyle=(0, (4, 2)),
              label="fundamental complementarity")[0]
axB.legend(handles=[h1, h2, h3], loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=1,
           fontsize=7, frameon=False, handletextpad=0.5)

out_pdf = os.path.join(_here, "fig0_illustrative_example.pdf")
out_png = os.path.join(_here, "fig0_illustrative_example.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"[make_fig_illustrative] wrote {out_pdf} and {out_png}")
