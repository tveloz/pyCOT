"""
analyze_agro_stages.py — Deep-dive report on data/Examples_tests/FarmVariants/Farm_agro_stages.txt,
the four-stage farm design (pure agro / +cows / +chickens / full) built with
scripts/analyze_farm_dependencies.py's sibling tooling (decomp/dependency.py,
decomp/witness.py).

HOW TO RUN
----------
  python projects/Decomposition_Theorem/scripts/analyze_agro_stages.py
"""
from __future__ import annotations

import os, sys
_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_cot_gen_proj = os.path.normpath(os.path.join(_here, "..", "..", "COT_fundamental_Generators"))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, _cot_gen_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

from pyCOT.io.functions      import read_txt
from cot_gen.io_pyCOT        import build_rndata
from cot_gen.erc             import compute_ercs
from cot_gen.hierarchy       import build_hierarchy
from cot_gen.synergy         import compute_synergies_basis_first
from cot_gen.complementarity import compute_complementarities
from cot_gen.epm             import compute_epms, compute_espm
from cot_gen.deep_report     import build_so_lattice

from decomp.bridge      import build_full_stoich, so_domain
from decomp.hierarchy   import decompose_hierarchy
from decomp.dependency  import build_dependency_dag, explain_circuit_failure
from decomp.witness     import compute_witness

NETWORK = "data/Examples_tests/FarmVariants/Farm_agro_stages.txt"
ESPM_MAX_ORDER = 100

net_path = os.path.join(_repo, NETWORK)
NET_ID = os.path.splitext(os.path.basename(NETWORK))[0]

SEP = "=" * 78
print(SEP)
print(f"Agro-stages deep dive — {NET_ID}")
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

nodes_sorted = sorted(results.keys(), key=lambda sp: bin(sp).count('1'))

for sp in nodes_sorted:
    r = results[sp]
    names = sorted(rn_data.bitset_to_names(sp))
    print(SEP)
    print(f"|X|={bin(sp).count('1')}  is_organization={r.is_organization}")
    print(f"  species: {names}")
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
        if w is not None:
            print(w.report(rn_data))
    print()
