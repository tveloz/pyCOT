"""
make_fig_toy.py -- static figure for the toy endosymbiosis model: host
network (left), symbiont network (right), interface reactions (crossing
arrows), and the species reachable only once merged highlighted.

Positions are hand-placed (not auto-laid-out) -- the network is small and
fixed, and explicit placement avoids the overlap/off-canvas problems a
generic layout produced on the first attempt.
"""
from __future__ import annotations

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch

_here = os.path.dirname(os.path.abspath(__file__))
_fig_dir = os.path.join(_here, "..", "figures")
os.makedirs(_fig_dir, exist_ok=True)

FOOD_COLOR = "#74b9ff"
SPECIES_COLOR = "#dfe6e9"
NEW_COLOR = "#fdcb6e"
REACTION_COLOR = "#ffeaa7"
EDGE_COLOR = "#636e72"
INTERFACE_COLOR = "#e17055"

# (x, y) hand layout. Host occupies x in [0, 4.4], symbiont x in [6.4, 11.2].
SPECIES_POS = {
    # host food (left column)
    "glc_h": (0.0, 4.6), "adp_h": (0.0, 3.6), "pi_h": (0.0, 2.6), "nad_h": (0.0, 1.6),
    # host internal
    "pyr_h": (2.2, 3.4), "atp_h": (2.2, 2.4), "nadh_h": (2.2, 1.4),
    "lac_h": (4.4, 1.0), "fescluster_h": (4.4, 3.8), "macromol_h": (4.4, 2.6),
    # symbiont food (right column, mirrored)
    "glc_s": (11.2, 5.4), "o2_s": (11.2, 4.4), "adp_s": (11.2, 3.4),
    "pi_s": (11.2, 2.4), "nad_s": (11.2, 1.4), "fe_s": (11.2, 0.4), "s_s": (11.2, -0.6),
    # symbiont internal
    "pyr_s": (9.0, 3.8), "atp_s": (9.0, 2.8), "nadh_s": (9.0, 1.8), "co2_s": (9.0, 0.8),
    "fescluster_s": (6.8, 1.6),
}
REACTIONS = [
    ("GLYC", ["glc_h", "adp_h", "pi_h", "nad_h"], ["pyr_h", "atp_h", "nadh_h"], (1.1, 4.0)),
    ("FERM", ["pyr_h", "nadh_h"], ["lac_h", "nad_h"], (3.3, 1.9)),
    ("MAINT", ["atp_h"], ["adp_h", "pi_h"], (1.1, 3.0)),
    ("BIOSYN", ["atp_h", "fescluster_h", "nad_h"], ["macromol_h"], (3.3, 3.3)),
    ("GLYC_S", ["glc_s", "adp_s", "pi_s", "nad_s"], ["pyr_s", "atp_s", "nadh_s"], (10.1, 4.0)),
    ("OXID", ["pyr_s", "nad_s"], ["co2_s", "nadh_s"], (7.9, 2.4)),
    ("ETC", ["nadh_s", "o2_s", "adp_s", "pi_s"], ["nad_s", "atp_s"], (10.1, 2.4)),
    ("ISC", ["fe_s", "s_s", "atp_s"], ["fescluster_s", "adp_s"], (7.9, 0.3)),
]
INTERFACE = [("pyr_h", "pyr_s"), ("atp_s", "atp_h"), ("fescluster_s", "fescluster_h")]
HOST_FOOD = {"glc_h", "adp_h", "pi_h", "nad_h"}
SYMB_FOOD = {"glc_s", "o2_s", "adp_s", "pi_s", "nad_s", "fe_s", "s_s"}
NEW_ONLY_MERGED = {"fescluster_h", "macromol_h"}

fig, ax = plt.subplots(figsize=(12.5, 6.5))

for name, subs, prods, (rx, ry) in REACTIONS:
    ax.add_patch(FancyBboxPatch((rx - 0.42, ry - 0.17), 0.84, 0.34,
                                 boxstyle="round,pad=0.02,rounding_size=0.05",
                                 facecolor=REACTION_COLOR, edgecolor="#2d3436",
                                 linewidth=0.9, zorder=3))
    ax.text(rx, ry, name, ha="center", va="center", fontsize=6.6, fontweight="bold", zorder=4)
    for s in subs:
        x, y = SPECIES_POS[s]
        ax.add_patch(FancyArrowPatch((x, y), (rx, ry), arrowstyle="-|>", mutation_scale=8,
                                      color=EDGE_COLOR, linewidth=0.75, shrinkA=13, shrinkB=10, zorder=2))
    for s in prods:
        x, y = SPECIES_POS[s]
        ax.add_patch(FancyArrowPatch((rx, ry), (x, y), arrowstyle="-|>", mutation_scale=8,
                                      color=EDGE_COLOR, linewidth=0.75, shrinkA=10, shrinkB=13, zorder=2))

for s, (x, y) in SPECIES_POS.items():
    if s in NEW_ONLY_MERGED:
        color = NEW_COLOR
    elif s in HOST_FOOD or s in SYMB_FOOD:
        color = FOOD_COLOR
    else:
        color = SPECIES_COLOR
    ax.add_patch(Circle((x, y), 0.30, facecolor=color, edgecolor="#2d3436", linewidth=1.1, zorder=5))
    label = s.replace("_h", "").replace("_s", "").replace("cluster", "-clus")
    ax.text(x, y, label, ha="center", va="center", fontsize=6.4, zorder=6)

for sp_a, sp_b in INTERFACE:
    # Draw in the REAL direction (sp_a -> sp_b, source -> destination, as
    # given in INTERFACE) -- an earlier version normalised to always draw
    # host->symbiont regardless of the tuple's actual order, which silently
    # reversed the arrowhead on 2 of the 3 real interface reactions (ATP and
    # the Fe-S cluster both flow symbiont->host, only pyruvate flows
    # host->symbiont). Caught by visually inspecting the rendered figure.
    x1, y1 = SPECIES_POS[sp_a]
    x2, y2 = SPECIES_POS[sp_b]
    arrow = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=12,
                             color=INTERFACE_COLOR, linewidth=1.8, linestyle=(0, (5, 2)),
                             connectionstyle="arc3,rad=0.12", shrinkA=16, shrinkB=16, zorder=7)
    ax.add_patch(arrow)

ax.text(2.2, 5.9, "HOST\n(glycolysis + fermentation)", ha="center", fontsize=11, fontweight="bold")
ax.text(9.0, 6.0, "SYMBIONT\n(oxidative phosphorylation + Fe-S biogenesis)", ha="center", fontsize=11, fontweight="bold")

ax.scatter([], [], marker='o', color=FOOD_COLOR, edgecolor="#2d3436", s=100, label="food species")
ax.scatter([], [], marker='o', color=SPECIES_COLOR, edgecolor="#2d3436", s=100, label="internal species")
ax.scatter([], [], marker='o', color=NEW_COLOR, edgecolor="#2d3436", s=100,
           label="reachable ONLY once merged")
ax.plot([], [], color=INTERFACE_COLOR, linewidth=2.0, linestyle=(0, (5, 2)), label="interface (transport) reaction")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.1), ncol=4, fontsize=8.5, frameon=False)

ax.set_title("Toy endosymbiosis model: host alone = 1 EPM (8 sp.), symbiont alone = 1 EPM (12 sp.),\n"
              "merged = 1 EPM (22 sp., 100% hybrid) -- interface unlocks 2 species neither partner reaches alone",
              fontsize=10)
ax.set_xlim(-0.8, 12.0)
ax.set_ylim(-1.3, 6.6)
ax.axis("off")

out_pdf = os.path.join(_fig_dir, "fig_toy_model.pdf")
out_png = os.path.join(_fig_dir, "fig_toy_model.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"wrote {out_pdf} and {out_png}")
