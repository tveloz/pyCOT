"""
make_fig_fusigenic.py -- bar chart of fusigenic_inflow_search.py results:
(a) toy host exhaustive sweep, (b) e_coli_core exhaustive sweep,
(c) iAF1260 targeted trace-metal sweep -- numbers transcribed from
outputs/*.log.
"""
from __future__ import annotations

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_here = os.path.dirname(os.path.abspath(__file__))
_fig_dir = os.path.join(_here, "..", "figures")
os.makedirs(_fig_dir, exist_ok=True)

TOY_HOST = [("fescluster_h", 2.00), ("macromol_h", 1.00), ("atp_h", 0.0),
            ("lac_h", 0.0), ("nadh_h", 0.0), ("pyr_h", 0.0)]

ECOLI_CORE = [("2pg_c", 41.15), ("3pg_c", 41.15), ("pep_c", 41.15),
              ("acon_C_c", 23.13), ("cit_c", 23.13), ("icit_c", 23.13),
              ("adp_c", 14.80), ("atp_c", 14.80), ("accoa_c", 13.42),
              ("dhap_c", 10.07), ("fdp_c", 10.07), ("g3p_c", 10.07)]

IAF1260 = [("fe2_e", 10.02), ("ni2_e", 6.00), ("cobalt2_e", 6.00),
           ("mn2_e", 6.00), ("mg2_e", 6.00), ("zn2_e", 6.00), ("k_e", 6.00),
           ("na1_e", 4.12), ("mobd_e", 4.01), ("tungs_e", 4.01),
           ("ca2_e", 4.01), ("so4_e", 4.01)]

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3))

panels = [
    (axes[0], TOY_HOST, "(a) toy host, exhaustive\n(6 candidates)", "#3498db"),
    (axes[1], ECOLI_CORE, "(b) e_coli_core, exhaustive\n(top 12 of 65 candidates)", "#27ae60"),
    (axes[2], IAF1260, "(c) iAF1260, targeted trace-metal sweep\n(12 of 15 candidates shown)", "#9b59b6"),
]

for ax, data, title, color in panels:
    names = [d[0] for d in data]
    scores = [d[1] for d in data]
    y = range(len(names))
    ax.barh(list(y), scores, color=color)
    ax.set_yticks(list(y))
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("fusion score", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.spines[['top', 'right']].set_visible(False)

fig.tight_layout()

out_pdf = os.path.join(_fig_dir, "fig_fusigenic_search.pdf")
out_png = os.path.join(_fig_dir, "fig_fusigenic_search.png")
fig.savefig(out_pdf, bbox_inches="tight")
fig.savefig(out_png, dpi=200, bbox_inches="tight")
print(f"wrote {out_pdf} and {out_png}")
