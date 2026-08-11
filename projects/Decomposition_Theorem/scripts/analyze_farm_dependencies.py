"""
analyze_farm_dependencies.py — Deep-dive report on every semi-organization
of a network: E/F/circuit decomposition, the food-dependency DAG among
fragile circuits (decomp/dependency.py), a diagnosis of WHY any
non-self-maintaining circuit fails, and -- for every SO that IS a full
organization -- a concrete witness flux vector plus net production per
species (decomp/witness.py).

Built for data/Examples_tests/FarmVariants/Farm_r7_diff_fixed.txt (the r7-diff farm
variant with its R8 species-tokenization bug fixed -- see the sibling
Farm_r7_diff.txt for the original), to check whether the network actually
supports the three intended operating regimes:
  - "chickenless": grain/straw economy running, no chickens/eggs.
  - "grainless"   : chicken/egg economy running, no grain/straw (possibly
                    zero net production on that cycle -- self-maintaining
                    but not overproducing -- which is still a valid
                    organization, just not a growing one).
  - "full"        : both economies running together.

HOW TO RUN
----------
  python projects/Decomposition_Theorem/scripts/analyze_farm_dependencies.py
"""
from __future__ import annotations

import os, sys
_here = os.path.dirname(os.path.abspath(__file__))
_cot_gen_proj = os.path.normpath(os.path.join(_here, "..", "..", "COT_Fundamental_Generators_Exploration"))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_cot_gen_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

from pyCOT.io.functions      import read_txt
from pyCOT.analysis.organizations import (
    build_rndata,
    compute_ercs,
    build_hierarchy,
    compute_synergies_basis_first,
    compute_complementarities,
    compute_epms, compute_espm,
)
from cot_gen.deep_report     import build_so_lattice

from pyCOT.analysis.decomposition import (
    build_full_stoich, so_domain,
    decompose_hierarchy,
    build_dependency_dag, explain_circuit_failure,
    compute_witness,
)

# ╔══════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION                                                       ║
# ╚══════════════════════════════════════════════════════════════════════╝
NETWORK = "data/Examples_tests/FarmVariants/Farm_agro_stages.txt"
ESPM_MAX_ORDER = 100

CHICKEN_SPECIES = {"chickens", "eggs"}
GRAIN_SPECIES = {"grain", "straw"}

# ╔══════════════════════════════════════════════════════════════════════╗
# ║  Pipeline                                                            ║
# ╚══════════════════════════════════════════════════════════════════════╝
net_path = os.path.join(_repo, NETWORK)
NET_ID = os.path.splitext(os.path.basename(NETWORK))[0]

SEP = "=" * 78
print(SEP)
print(f"Semi-org / dependency deep dive — {NET_ID}")
print(SEP)

rn_pycot = read_txt(net_path)
rn_data = build_rndata(rn_pycot, network_id=NET_ID)
ercs = compute_ercs(rn_data, verify=True)
hier = build_hierarchy(ercs)
syn = compute_synergies_basis_first(ercs, hier)
comp = compute_complementarities(ercs, hier, syn)
epm_result = compute_epms(rn_data, ercs, hier, syn, comp, verbose=False)
espm_result = compute_espm(rn_data, ercs, hier, syn, comp, epm_result,
                            max_order=ESPM_MAX_ORDER, verbose=False)

so_order = {sp: 0 for sp in epm_result.all_epm_masks}
for order, masks in espm_result.espm_by_order.items():
    for sp in masks:
        so_order[sp] = order
so_lattice = build_so_lattice(espm_result.all_so_masks, so_order)

S_full = build_full_stoich(rn_pycot)
results = decompose_hierarchy(so_lattice, rn_data, S_full, verbose=False)
print(f"{len(results)} semi-organizations total\n")

nodes_sorted = sorted(results.keys(),
                       key=lambda sp: (so_lattice.order_of.get(sp, 0), bin(sp).count('1')))

regime_hits: dict[str, list[int]] = {"chickenless": [], "grainless": [], "full": [], "neither": []}

for sp in nodes_sorted:
    r = results[sp]
    order = so_lattice.order_of.get(sp, 0)
    names = set(rn_data.bitset_to_names(sp))
    has_chicken = bool(names & CHICKEN_SPECIES)
    has_grain = bool(names & GRAIN_SPECIES)
    if has_chicken and has_grain:
        regime_hits["full"].append(sp)
    elif has_grain:
        regime_hits["chickenless"].append(sp)
    elif has_chicken:
        regime_hits["grainless"].append(sp)
    else:
        regime_hits["neither"].append(sp)

    print(SEP)
    print(f"order {order}  |X|={bin(sp).count('1')}  is_organization={r.is_organization}")
    print(f"  species: {sorted(names)}")
    print(f"  E (catalysts):  {rn_data.bitset_to_names(r.E_mask)}")
    print(f"  F (overprod.):  {rn_data.bitset_to_names(r.F_mask)}")

    for i, c in enumerate(r.circuits):
        print(f"  D{i+1}: {rn_data.bitset_to_names(c.species_mask)}  "
              f"self_maintaining={c.is_self_maintaining}")
        if not c.is_self_maintaining:
            domain = so_domain(sp, rn_data, S_full)
            print(explain_circuit_failure(c, domain, rn_data))

    if r.circuits:
        dag = build_dependency_dag(r, rn_data)
        print("  dependency structure:")
        print(dag.summary(rn_data))

    if r.is_organization:
        domain = so_domain(sp, rn_data, S_full)
        w = compute_witness(domain, r.F_mask)
        if w is None:
            print("  ** expected a witness (is_organization=True) but the global LP "
                  "was infeasible -- investigate, this should not happen. **")
        else:
            print(w.report(rn_data))
    print()

# ── Regime summary ──────────────────────────────────────────────────────
print(SEP)
print("Regime classification (by species membership, not by is_organization)")
print(SEP)
for label, sps in regime_hits.items():
    if not sps:
        continue
    print(f"\n{label}: {len(sps)} SO(s)")
    for sp in sorted(sps, key=lambda s: (so_lattice.order_of.get(s, 0), bin(s).count('1'))):
        r = results[sp]
        tag = "ORGANIZATION" if r.is_organization else "semi-org only"
        print(f"  order {so_lattice.order_of.get(sp, 0):>2}  |X|={bin(sp).count('1'):>2}  "
              f"{tag:14s}  {sorted(rn_data.bitset_to_names(sp))}")
