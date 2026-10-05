"""
validate_small.py -- correctness check for gen_v2 against the SAME
oracle-validated small-network corpus so_search.py is checked against
(gold networks, the worked example, BIOMD0000000091, BIOMD0000000999).

Must pass before gen_v2 is trusted on anything larger, where no
brute-force oracle is feasible and old-vs-new (compare_old_new.py) is the
only available cross-check.

Checks, per network:
  1. gen_v2's elementary-SO set matches the oracle exactly.
  2. gen_v2's full SO set matches the oracle exactly (these small networks
     have no pure latent joins, so no reconciliation needed here -- see
     so_search.py's latent_join() docstring for what that means).
  3. gen_v2's result matches the PRODUCTION engine's result exactly (the
     two should agree completely on this corpus, since the only
     behavioral difference -- the closure-complete synergy fix -- was
     already shown to change nothing on networks of this kind).

Run: python -m gen_v2.validate_small   (from the project root, or via
     projects/COT_Fundamental_Generators_Exploration/scripts conventions)
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyCOT.analysis.organizations.cot_types import RNData
from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity import compute_complementarities
from pyCOT.analysis.organizations.so_search import compute_elementary_sos, compute_so_hierarchy
from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph

from gen_v2.engine import explore

from tests.gold_networks import ALL_GOLD
from worked_example import worked_rndata
from oracles.so_oracle import elementary_so_oracle, so_hierarchy_oracle


def _build_rndata_from_gold(net) -> RNData:
    supp_q = tuple(s & ~net.E0_mask for s in net.supp)
    prod_q = tuple(p & ~net.E0_mask for p in net.prod)
    inv = tuple(tuple(r for r, s in enumerate(supp_q) if (s >> i) & 1)
                for i in range(net.n_species))
    return RNData(
        n_species=net.n_species, species_names=tuple(net.species),
        species_index=tuple((n, i) for i, n in enumerate(net.species)),
        n_reactions=len(supp_q), reaction_names=tuple(f"r{i}" for i in range(len(supp_q))),
        supp_raw=tuple(net.supp), prod_raw=tuple(net.prod),
        E0_mask=net.E0_mask, supp_q=supp_q, prod_q=prod_q, species_to_reactions=inv,
    )


def check_network(name: str, rn, verbose=True) -> bool:
    ercs = compute_ercs(rn, verify=False)
    if not ercs:
        if verbose:
            print(f"  {name}: skip (no ERCs)")
        return True
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)

    # production
    elem_old = compute_elementary_sos(rn, ercs, hier, syn, comp, verbose=False)
    hr_old = compute_so_hierarchy(rn, ercs, hier, syn, comp, elem_old,
                                   use_vertical_lift=True, max_order=15, verbose=False)
    old_all = set(hr_old.all_so_masks)
    old_elem = set(elem_old.all_elementary_masks)

    # gen_v2
    g = FundamentalGraph(ercs, hier, syn, comp)
    result_new = explore(g)
    new_all = set(result_new.all_so_masks)
    new_elem = set(result_new.elementary_masks)

    e0 = rn.E0_mask
    old_all_nz = {m for m in old_all if m != e0}
    new_all_nz = {m for m in new_all if m != e0}
    old_elem_nz = {m for m in old_elem if m != e0}
    new_elem_nz = {m for m in new_elem if m != e0}

    ok = True
    if old_all_nz != new_all_nz:
        ok = False
        print(f"  {name}: FAIL old-vs-new mismatch  "
              f"old-only={sorted(hex(m) for m in old_all_nz - new_all_nz)}  "
              f"new-only={sorted(hex(m) for m in new_all_nz - old_all_nz)}")

    # Oracle, where feasible. The threshold is lower than the 24 used
    # elsewhere: BIOMD0000000109 (19 ERCs) was timed directly and its
    # brute-force oracle alone exceeds 90s (unrelated to gen_v2 -- some
    # networks' structure makes the oracle's own per-subset work heavier
    # than others at the same ERC count). The new-vs-old comparison above
    # is still a strong, independent check for networks skipped here.
    if len(ercs) <= 15:
        oracle_elem = set(elementary_so_oracle(rn, ercs))
        oracle_all = set(so_hierarchy_oracle(rn, ercs)["all_sos"])
        oracle_elem_nz = {m for m in oracle_elem if m != e0}
        oracle_all_nz = {m for m in oracle_all if m != e0}
        if new_elem_nz != oracle_elem_nz:
            ok = False
            print(f"  {name}: FAIL gen_v2 elementary vs oracle  "
                  f"ours-only={sorted(hex(m) for m in new_elem_nz - oracle_elem_nz)}  "
                  f"oracle-only={sorted(hex(m) for m in oracle_elem_nz - new_elem_nz)}")
        if new_all_nz != oracle_all_nz:
            ok = False
            print(f"  {name}: FAIL gen_v2 full-SO vs oracle  "
                  f"ours-only={sorted(hex(m) for m in new_all_nz - oracle_all_nz)}  "
                  f"oracle-only={sorted(hex(m) for m in oracle_all_nz - new_all_nz)}")

    if ok and verbose:
        print(f"  {name}: OK  (n_ercs={len(ercs)}  n_so={len(new_all_nz)})")
    return ok


def main():
    print("gen_v2 validation against the oracle-checked small-network corpus")
    print("=" * 70)
    all_ok = True

    for net in ALL_GOLD:
        rn = _build_rndata_from_gold(net)
        all_ok &= check_network(f"gold_{net.name}", rn)

    all_ok &= check_network("worked_example", worked_rndata())

    biomd_cases = ["BIOMD0000000091", "BIOMD0000000999", "BIOMD0000000109"]
    categories = {"BIOMD0000000091": "BioMD_other", "BIOMD0000000999": "BioMD_signaling",
                  "BIOMD0000000109": "BioMD_cell_cycle"}
    for name in biomd_cases:
        path = os.path.join(_repo, "data", "biochemical_databases", categories[name], f"{name}.txt")
        if not os.path.exists(path):
            print(f"  {name}: skip (file not found)")
            continue
        from pyCOT.io.functions import read_txt
        from pyCOT.analysis.organizations.io_pyCOT import build_rndata
        rn = build_rndata(read_txt(path), network_id=name)
        all_ok &= check_network(name, rn)

    print("=" * 70)
    print("ALL PASSED" if all_ok else "SOME FAILED -- do not trust gen_v2 on larger networks yet")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
