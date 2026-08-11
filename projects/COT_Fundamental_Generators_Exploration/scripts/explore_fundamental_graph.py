"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  explore_fundamental_graph.py — Interactive manual exploration of the       ║
║                                  ERC fundamental graph                       ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT IT DOES
------------
Loads one reaction network, computes its ERC hierarchy, fundamental
synergies and fundamental complementarities, then drops you into an
interactive Python console with a `gx` object (a GraphExplorer,
cot_gen/explorer.py) already loaded and ready to use.

From that console you can, by hand:
  - Pick an ERC and inspect it, or see its local neighborhood (containment /
    synergy / complementarity) up to a chosen radius.
  - Build a generator step by step: add ERCs, see which species/reactions
    it reaches, what it still requires, whether it has become
    semi-self-maintaining (SSM) -- i.e. formed a semi-organization.
  - See exactly which other ERCs are connected to your current generator
    and *why* (producer of an open requirement, consumer of something
    already produced, synergy partner, or hierarchy ancestor).
  - Remove an ERC to roll back, or save/restore named checkpoints to try
    alternative extension paths from the same point.
  - Draw the current generator + its local neighborhood (small, always
    readable), or the whole network with the generator highlighted in
    context (size-capped, since full hierarchies can be large).

WHAT IT DOES NOT DO
-------------------
This is a manual exploration tool, not the validated search. Nothing here
enforces canonical ordering, minimality, or fundamentality -- you can add
ERCs in any order, including ones that would never appear in an
irreducible generator. For the exhaustive, validated computation of all
EPMs/ESPMs, use run_network.py / cot_gen.epm.compute_epms.

HOW TO RUN
----------
  1. Edit the CONFIGURATION section below (pick a network).
  2. Press ▶ (play) in VS Code, or run:
       python projects/COT_Fundamental_Generators_Exploration/scripts/explore_fundamental_graph.py
  3. At the ">>>" prompt, try (for example):
       gx.describe(0)
       gx.neighbors(0, radius=2)
       gx.local_view(0, radius=2)
       gx.add(0); gx.candidates(); gx.summary()
       gx.context_view()
       help_gx()          # reprints this cheat sheet

Or, from a Jupyter notebook / IPython session, skip this script entirely:
    from cot_gen.explorer import GraphExplorer
    gx = GraphExplorer.from_network("e_coli_core")
"""

# ── Path setup (do not edit) ──────────────────────────────────────────────────
from __future__ import annotations
import os, sys, code

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — edit this block, then press Play                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# Short catalogue name (e.g. "e_coli_core", "BIOMD0000000237") or a direct
# path to a .txt reaction-network file.
NETWORK = "BIOMD0000000237"
NETWORK   = "e_coli_core"
NETWORK   = "data\\Examples_tests\\LUCA\\KEGG_data.txt"


# Start the generator with these ERC indices already added (empty = start
# from nothing and pick your first ERC interactively with gx.add(...)).
INITIAL_GENERATOR: list[int] = []

# ╔══════════════════════════════════════════════════════════════════════════════╝

from cot_gen.explorer import GraphExplorer

SEP = "=" * 72
print(SEP)
print(f"Fundamental Graph Explorer  —  {NETWORK}")
print(SEP)

gx = GraphExplorer.from_network(NETWORK)

for _idx in INITIAL_GENERATOR:
    gx.add(_idx)


def help_gx():
    print(__doc__)


print("""
Ready. `gx` is your GraphExplorer. Try:
    gx.describe(0)                 gx.candidates()
    gx.neighbors(0, radius=2)      gx.summary()
    gx.add(0)                      gx.degree_stats()
    gx.remove(0)                   gx.history()
    gx.checkpoint("a")             gx.restore("a")
    gx.synergy_level_distribution()   -- which hierarchy depth has the most synergies
    gx.generator_view()               -- exactly the generator's subgraph (+ candidates)
    gx.local_view(0, radius=2)        -- one ERC's neighborhood
    gx.context_view()                 -- whole-network hierarchy, generator lit up
    gx.context_view(synergy_layout="interleaved")  -- synergies banded by hierarchy depth
    help_gx()   -- reprint the full cheat sheet
    exit()      -- leave the console
""")

if __name__ == "__main__":
    code.interact(banner="", local=dict(globals(), **locals()), exitmsg="")
