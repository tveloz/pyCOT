"""
validate_decomposition.py — Correctness checks for the decomp package.

Run directly: python projects/Decomposition_Theorem/tests/validate_decomposition.py

Three independent checks, on every EPM + a sample of ESPMs of a test network:

  1. Oracle check: decomp's is_organization must agree with a DIRECT,
     decomposition-free self-maintenance LP on the full closed species set
     (Persistent_Modules_Generator.check_self_maintenance on X_full as a
     whole). This is exactly Theorem 2.16's claim, tested empirically.
  2. Consistency check: decompose_hierarchy's incrementally-shortcut result
     must be IDENTICAL (same E, F, circuits, is_organization) to calling
     the standalone decompose() fresh on every node — i.e. the monotonicity
     / circuit-reuse shortcuts in hierarchy.py never change the answer,
     only how much LP work is spent getting there.
  3. Structural sanity: E, F, and every circuit are pairwise disjoint and
     their union is exactly X_full; every circuit is non-empty.
"""
from __future__ import annotations

import os
import sys
import time

_here = os.path.dirname(os.path.abspath(__file__))
_cot_gen_proj = os.path.normpath(os.path.join(_here, "..", "..", "COT_Fundamental_Generators_Exploration"))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_cot_gen_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations import check_self_maintenance

from pyCOT.analysis.organizations import (
    build_rndata,
    compute_ercs,
    build_hierarchy,
    compute_synergies_basis_first,
    compute_complementarities,
    compute_epms, compute_espm,
)
from cot_gen.deep_report import build_so_lattice

from pyCOT.analysis.decomposition import build_full_stoich, decompose, decompose_hierarchy

NETWORK = "e_coli_core"
_DATA_ROOT = os.path.join(_repo, "data", "biochemical_databases", "biomodels_all_txt")


def _find_network(name: str) -> str:
    for root, _dirs, files in os.walk(_DATA_ROOT):
        for fname in files:
            if fname == f"{name}.txt" or fname == f"bigg_{name}.txt":
                return os.path.join(root, fname)
    raise FileNotFoundError(name)


def main():
    net_path = _find_network(NETWORK)
    print(f"Loading {NETWORK} ...")
    rn_pycot = read_txt(net_path)
    rn_data = build_rndata(rn_pycot, network_id=NETWORK)

    ercs = compute_ercs(rn_data, verify=True)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)
    epm_result = compute_epms(rn_data, ercs, hier, syn, comp)
    espm_result = compute_espm(rn_data, ercs, hier, syn, comp, epm_result, max_order=6)

    print(f"ERCs={len(ercs)}  EPMs={len(epm_result.all_epm_masks)}  "
          f"SOs(order<=6)={len(espm_result.all_so_masks)}")

    so_order = {sp: 0 for sp in espm_result.epm_masks}
    for order, masks in espm_result.espm_by_order.items():
        for sp in masks:
            so_order[sp] = order
    so_lattice = build_so_lattice(espm_result.all_so_masks, so_order)

    S_full = build_full_stoich(rn_pycot)

    # ---- Check 2 setup: hierarchy vs standalone -----------------------------
    t0 = time.perf_counter()
    hier_results = decompose_hierarchy(so_lattice, rn_data, S_full, verbose=False)
    t_hier = time.perf_counter() - t0

    t0 = time.perf_counter()
    standalone_results = {sp: decompose(sp, rn_data, S_full) for sp in so_lattice.nodes}
    t_standalone = time.perf_counter() - t0

    mismatches = []
    for sp in so_lattice.nodes:
        a, b = hier_results[sp], standalone_results[sp]
        if (a.E_mask, a.F_mask, a.is_organization) != (b.E_mask, b.F_mask, b.is_organization):
            mismatches.append(sp)
            continue
        a_circ = sorted((c.species_mask, c.reaction_ids, c.is_self_maintaining) for c in a.circuits)
        b_circ = sorted((c.species_mask, c.reaction_ids, c.is_self_maintaining) for c in b.circuits)
        if a_circ != b_circ:
            mismatches.append(sp)

    print(f"\n[Check 2] hierarchy vs standalone: {len(mismatches)} mismatches "
          f"out of {len(so_lattice.nodes)} nodes")
    print(f"          timing: hierarchy={t_hier*1000:.1f}ms  standalone={t_standalone*1000:.1f}ms")
    if mismatches:
        print(f"          MISMATCHED sp_masks (first 5): {mismatches[:5]}")

    # ---- Check 1: oracle agreement ------------------------------------------
    oracle_mismatches = []
    for sp in so_lattice.nodes:
        result = hier_results[sp]
        x_full_mask = sp | rn_data.E0_mask
        names = rn_data.bitset_to_names(x_full_mask)
        species_objs = [rn_pycot.get_species(n) for n in names]
        is_sm, _flux, _prod = check_self_maintenance(species_objs, rn_pycot)
        if bool(is_sm) != bool(result.is_organization):
            oracle_mismatches.append((sp, is_sm, result.is_organization))

    print(f"\n[Check 1] decomp vs direct-LP oracle: {len(oracle_mismatches)} mismatches "
          f"out of {len(so_lattice.nodes)} nodes")
    if oracle_mismatches:
        for sp, oracle, decomp_verdict in oracle_mismatches[:5]:
            print(f"          sp={sp}  oracle={oracle}  decomp={decomp_verdict}")

    # ---- Check 3: structural sanity -----------------------------------------
    sanity_fail = []
    for sp, result in hier_results.items():
        x_full = sp | rn_data.E0_mask
        union = result.E_mask | result.F_mask
        for c in result.circuits:
            if c.species_mask == 0:
                sanity_fail.append((sp, "empty circuit"))
            if union & c.species_mask:
                sanity_fail.append((sp, "circuit overlaps E/F"))
            union |= c.species_mask
        if union != x_full:
            sanity_fail.append((sp, f"union {union:#x} != X_full {x_full:#x}"))

    print(f"\n[Check 3] structural sanity: {len(sanity_fail)} failures "
          f"out of {len(so_lattice.nodes)} nodes")
    if sanity_fail:
        for sp, msg in sanity_fail[:5]:
            print(f"          sp={sp}  {msg}")

    n_org = sum(1 for r in hier_results.values() if r.is_organization)
    print(f"\nSummary: {n_org}/{len(hier_results)} SOs are full organizations.")

    ok = not mismatches and not oracle_mismatches and not sanity_fail
    print("\n" + ("ALL CHECKS PASSED" if ok else "FAILURES DETECTED — see above"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
