"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  maxraf_analysis.py — RAF-style sustainability pruning for the ERC graph    ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT IT COMPUTES
----------------
An ERC is "sustainable" only if every species it requires can eventually be
produced by SOME combination of other sustainable ERCs. This is exactly the
Hordijk-Steel RAF (Reflexively Autocatalytic, Food-generated) condition,
applied to the ERC hierarchy with fundamental complementarity as the
producer/consumer relation -- with one correction the plain RAF definition
doesn't need: adding an ERC here can trigger a CHAIN of fundamental
SYNERGIES (not just complementarity), so "what a candidate pool can produce"
must be computed via the real synergy closure (`erc_syn_close`), not just
the raw union of each member's own product mask.

Algorithm (maxRAF, adapted)
----------------------------
  pool_0 = every ERC
  repeat:
      implied, _ = erc_syn_close(0, pool)   # full synergy closure of pool
      prod = union of prod_mask over `implied`
      pool_next = { i in implied : req_mask[i] subset of prod }
  until pool_next == pool
  return pool  (the maximal sustainable set, "maxRAF")

Termination and correctness
----------------------------
Let Phi(pool) = the filtering step above. Phi is monotone (bigger pool ->
bigger-or-equal `implied` -> bigger-or-equal `prod` -> the req-subset test
can only get EASIER, never harder). Starting from pool_0 = everyone (the
top of the subset lattice), Phi(pool_0) subset pool_0 trivially, and
monotonicity then forces the whole sequence pool_0 >= pool_1 >= pool_2 >= ...
to be non-increasing, so it converges in at most n_ercs iterations to the
GREATEST fixed point of Phi. Because RAF-type fixed points are closed under
union (a structural fact carried over from RAF theory: the union of two
self-sustaining sets is self-sustaining), maxRAF contains every genuine
EPM/ESPM's own ERC-set -- so any ERC excluded from maxRAF can be safely
excluded from the whole EPM/ESPM search: it could never appear in a real
persistent module.

WHAT IT DOES NOT DO
--------------------
Does not modify cot_gen/epm.py. This is a standalone, read-only
characterization tool -- see maxraf_filter.py in this same folder for the
version wired into the actual search (once validated).

HOW TO RUN
----------
  1. Edit NETWORK below.
  2. Run: python projects/COT_fundamental_Generators/scripts/maxraf_analysis.py
"""

from __future__ import annotations
import os, sys

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
for _s in (sys.stdout, sys.stderr):
    if hasattr(_s, "reconfigure"):
        _s.reconfigure(encoding="utf-8")

# ╔══════════════════════════════════════════════════════════════════════════════╗
NETWORK = "BIOMD0000000237"
# ╚══════════════════════════════════════════════════════════════════════════════╝

from pyCOT.io.functions import read_txt
from cot_gen.io_pyCOT import build_rndata
from cot_gen.erc import compute_ercs
from cot_gen.hierarchy import build_hierarchy
from cot_gen.synergy import compute_synergies_basis_first
from cot_gen.complementarity import compute_complementarities
from cot_gen.fundamental_graph import FundamentalGraph


def compute_max_sustainable_pool(g: FundamentalGraph) -> tuple[set[int], set[int]]:
    """Returns (surviving_indices, removed_indices)."""
    n = g.n
    pool = set(range(n))
    while True:
        implied, _ = g.erc_syn_close(0, pool)
        prod = 0
        for i in implied:
            prod |= g.prod_mask[i]
        next_pool = {i for i in implied if not (g.req_mask[i] & ~prod)}
        if next_pool == pool:
            break
        pool = next_pool
    return pool, set(range(n)) - pool


def _load_network(name: str):
    if os.path.isfile(name):
        path = name
    else:
        candidates = []
        for root, _, files in os.walk(os.path.join(_repo, "data", "biomodels")):
            for f in files:
                if f in (f"{name}.txt", f"bigg_{name}.txt"):
                    candidates.append(os.path.join(root, f))
        if not candidates:
            raise FileNotFoundError(name)
        path = candidates[0]
    return build_rndata(read_txt(path), network_id=name), path


def analyze(name: str, path_or_name: str | None = None):
    print(f"\n{'='*72}\n{name}\n{'='*72}")
    rn, path = _load_network(path_or_name or name)
    print(f"loaded: {path}")
    ercs = compute_ercs(rn)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn_result=syn)
    g = FundamentalGraph(ercs, hier, syn, comp)
    print(f"n_ercs={len(ercs)}  n_fundamental_syn={len(syn.fundamental)}  "
          f"n_fundamental_comp={len(comp.fundamental)}")

    survivors, removed = compute_max_sustainable_pool(g)
    print(f"\nmaxRAF survivors: {len(survivors)}/{len(ercs)} "
          f"({100*len(survivors)/max(len(ercs),1):.1f}%)  "
          f"removed: {len(removed)}")

    if not ercs:
        return

    sizes_survivors = sorted(bin(g.species_mask[i]).count('1') for i in survivors)
    sizes_removed = sorted(bin(g.species_mask[i]).count('1') for i in removed)
    if sizes_survivors:
        print(f"  survivor ERC sizes (species count): min={sizes_survivors[0]} "
              f"median={sizes_survivors[len(sizes_survivors)//2]} max={sizes_survivors[-1]}")
    if sizes_removed:
        print(f"  removed  ERC sizes (species count): min={sizes_removed[0]} "
              f"median={sizes_removed[len(sizes_removed)//2]} max={sizes_removed[-1]}")

    n_persistent_survivors = sum(1 for i in survivors if ercs[i].is_persistent())
    n_persistent_removed = sum(1 for i in removed if ercs[i].is_persistent())
    print(f"  persistent ERCs: {n_persistent_survivors} among survivors, "
          f"{n_persistent_removed} among removed "
          f"(sanity: removed should be 0, persistent ERCs are trivially sustainable)")

    # Hierarchy structure among survivors: do they form their own coherent
    # sub-hierarchy, or are they scattered leaves with no internal structure?
    surv_ancestors_within = 0
    surv_descendants_within = 0
    for i in survivors:
        surv_ancestors_within += len(hier.ancestors[i] & survivors)
        surv_descendants_within += len(hier.descendants[i] & survivors)
    avg_anc = surv_ancestors_within / max(len(survivors), 1)
    avg_desc = surv_descendants_within / max(len(survivors), 1)
    print(f"  hierarchy density among survivors: avg #ancestors-within-survivors="
          f"{avg_anc:.2f}  avg #descendants-within-survivors={avg_desc:.2f}")

    # Compare against overall hierarchy density for context
    all_anc = sum(len(hier.ancestors[i]) for i in range(len(ercs))) / max(len(ercs), 1)
    print(f"  (overall hierarchy avg #ancestors per ERC: {all_anc:.2f}, for comparison)")

    # Removed ERCs: how many are "isolated" (no fundamental comp/syn edges at all)
    # vs genuinely disqualified (had edges, but into an unsustainable chain)?
    isolated_removed = 0
    for i in removed:
        has_comp = any(fc.prod_idx == i or fc.cons_idx == i for fc in comp.fundamental)
        has_syn = bool(g.syn_from.get(i))
        if not has_comp and not has_syn:
            isolated_removed += 1
    print(f"  removed ERCs with NO fundamental comp/syn edges at all (fully isolated): "
          f"{isolated_removed}/{len(removed)}")


biomodels_dir = os.path.join(_repo, "data", "biomodels", "biomodels_all_txt")
bigg_dir = os.path.join(_repo, "data", "biomodels", "BiGG")
other_dir = os.path.join(_repo, "data", "biomodels", "BioMD_other")

if __name__ == "__main__":
    analyze("BIOMD0000000091", os.path.join(biomodels_dir, "BIOMD0000000091.txt"))
    analyze("BIOMD0000000237", os.path.join(other_dir, "BIOMD0000000237.txt"))
    analyze("e_coli_core", os.path.join(bigg_dir, "bigg_e_coli_core.txt"))
    analyze("iNJ661", os.path.join(bigg_dir, "bigg_iNJ661.txt"))
