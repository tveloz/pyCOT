"""
make_fig_native_networks.py -- regenerate the toy-model reaction-network
figures using pyCOT's OWN native visualization tool
(src/pyCOT/visualization/rn_visualize.py: create_bipartite_graph_from_rn +
graphviz), replacing the earlier hand-rolled matplotlib figure
(make_fig_toy.py) for the PNAS paper.

Uses create_bipartite_graph_from_rn (native pyCOT graph construction)
directly, then renders with graphviz using rankdir=LR (left-to-right) --
the print-suitable orientation -- rather than calling
rn_visualize_png_in_out verbatim, whose default top-to-bottom layout
produces an unusably tall/narrow image for the ~20-node merged network.
The graph CONSTRUCTION (which nodes/edges exist, from the real
ReactionNetwork object) is 100% native pyCOT; only the rendering
orientation is adjusted for the print page.

Produces one PNG+PDF per network: host alone, symbiont alone, merged.
"""
from __future__ import annotations

import os
import sys

os.environ["PATH"] = os.environ.get("PATH", "") + r";C:\Program Files\Graphviz\bin"

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))

from graphviz import Digraph
from pyCOT.io.functions import read_txt
from pyCOT.visualization.rn_visualize import create_bipartite_graph_from_rn

_toy = os.path.join(_here, '..', 'toy_model')
_fig_dir = os.path.join(_here, '..', 'figures')
os.makedirs(_fig_dir, exist_ok=True)

HOST_FOOD = {"glc_h", "adp_h", "pi_h", "nad_h"}
SYMB_FOOD = {"glc_s", "o2_s", "adp_s", "pi_s", "nad_s", "fe_s", "s_s"}
NEW_ONLY_MERGED = {"fescluster_h", "macromol_h"}

SPECIES_COLOR = "#dfe6e9"
FOOD_COLOR = "#74b9ff"
NEW_COLOR = "#fdcb6e"
REACTION_COLOR = "#ffeaa7"
EDGE_COLOR = "#636e72"


def render(net_path: str, out_name: str, food: set[str], highlight: set[str] | None = None,
           rankdir: str = "LR"):
    highlight = highlight or set()
    rn = read_txt(net_path, exact_names=True)
    # Native pyCOT bipartite-graph construction (species/reaction nodes,
    # reactant/product edges) -- unchanged from rn_visualize.py.
    graph, met_nodes, rxn_nodes = create_bipartite_graph_from_rn(rn)

    dot = Digraph(comment=out_name)
    dot.attr(rankdir=rankdir, nodesep="0.25", ranksep="0.45")

    for idx, (tipo, nombre) in enumerate(graph.nodes()):
        if tipo == 'specie':
            color = NEW_COLOR if nombre in highlight else (FOOD_COLOR if nombre in food else SPECIES_COLOR)
            dot.node(str(idx), nombre, shape="circle", style="filled", fillcolor=color,
                     fontsize="11", width="0.55", fixedsize="false")
        else:
            dot.node(str(idx), nombre, shape="box", style="filled", fillcolor=REACTION_COLOR,
                     fontsize="11")

    for src, dst in graph.edge_list():
        data = graph.get_edge_data(src, dst)
        label = str(data)
        if label == "1":
            dot.edge(str(src), str(dst), color=EDGE_COLOR)
        else:
            dot.edge(str(src), str(dst), label=label, color=EDGE_COLOR, fontsize="9")

    out_base = os.path.join(_fig_dir, out_name)
    dot.render(out_base, format="png", cleanup=True)
    dot.render(out_base, format="pdf", cleanup=True)
    print(f"wrote {out_base}.png / .pdf")


if __name__ == "__main__":
    render(os.path.join(_toy, "host_alone.txt"), "native_host", HOST_FOOD, rankdir="LR")
    render(os.path.join(_toy, "symbiont_alone.txt"), "native_symbiont", SYMB_FOOD, rankdir="LR")
    render(os.path.join(_toy, "merged.txt"), "native_merged", HOST_FOOD | SYMB_FOOD,
           NEW_ONLY_MERGED, rankdir="LR")
