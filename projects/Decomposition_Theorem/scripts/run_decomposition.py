"""
run_decomposition.py — Decompose every semi-organization of a network's
EPM/ESPM lattice (E / F / fragile circuits) and show how the decomposition
evolves order by order, from EPMs up to the maximal semi-organization.

HOW TO RUN
----------
  1. Edit the CONFIGURATION section below.
  2. Press Play in VS Code, or run:
       python projects/Decomposition_Theorem/scripts/run_decomposition.py
"""

# ── Path setup (do not edit) ──────────────────────────────────────────────────
from __future__ import annotations
import os, sys, time
_here = os.path.dirname(os.path.abspath(__file__))
_cot_gen_proj = os.path.normpath(os.path.join(_here, "..", "..", "COT_Fundamental_Generators_Exploration"))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_cot_gen_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")

# ── Imports (do not edit) ─────────────────────────────────────────────────────
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

from pyCOT.analysis.decomposition import build_full_stoich, decompose_hierarchy
from pyCOT.visualization.decomposition_viz import (
    plot_decomposition_evolution,
    plot_organization_hasse, plot_organization_chains,
)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — edit this block, then press Play                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

NETWORK        = "e_coli_core"
#NETWORK  = "BMID000000141754_url"
#NETWORK  = "BIOMD0000000185"
#NETWORK  = "iMM1415" $4000 reactions!
#NETWORK  = "iMM904" #(2072 reaction)
#NETWORK  = "iND750" #(1702 reactions)
#NETWORK  = "iNF517"
#NETWORK  = "iNJ661"     # genome-scale -- ESPM enumeration is not yet fast enough to
#NETWORK  = "iAF692"     # finish at these sizes; the Stage H2 maxSemiOrganization
                          # (a few ms) still works fine even here, but set
                          # COMPUTE_ESPM = False (or a low ESPM_MAX_ORDER) if you pick one.
#NETWORK  = "iSBO_1134"      #(3240 reactions)                      
#NETWORK  = "iSDY_1059"     #(3182 reactions)                          
#NETWORK  = "iSFV_1184"    #(3279 reactions)                           
#NETWORK  = "iSF_1195"    #(3900 reactions)
#NETWORK   = "data\\Examples_tests\\LUCA\\KEGG_data.txt"
#NETWORK = "BIOMD0000000183" #(352 REACTIONS)
NETWORK= "iAB_RBC_283"       #(645 REACTIONS) 
NETWORK = "iIT341" #(737 REACTIONS)
#NETWORK = "data/Examples_tests/FarmVariants/Farm.txt"
#NETWORK = "data/Examples_tests/FarmVariants/Farm_r7_diff.txt"
NETWORK = "data/Examples_tests/FarmVariants/Farm_agro_stages.txt"
ESPM_MAX_ORDER = 100
VERIFY         = True

# ── Visualizations ────────────────────────────────────────────────────────────
MAKE_PLOTS       = True
OUTPUT_DIR       = os.path.join(_here, "..", "outputs")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Script body — no need to edit below this line                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

_DATA_ROOT = os.path.join(_repo, "data", "biochemical_databases", "biomodels_all_txt")


def _discover_networks() -> dict[str, str]:
    seen: dict[str, tuple[str, int]] = {}
    for root, _dirs, files in os.walk(_DATA_ROOT):
        for fname in files:
            if not fname.endswith(".txt"):
                continue
            full = os.path.join(root, fname)
            name = os.path.splitext(fname)[0]
            if name.startswith("bigg_"):
                name = name[5:]
            if name not in seen or len(full) < seen[name][1]:
                seen[name] = (full, len(full))
    return {name: path for name, (path, _) in seen.items()}


_NETWORKS = _discover_networks()
if os.path.isfile(NETWORK):
    net_path, NET_ID = os.path.abspath(NETWORK), os.path.splitext(os.path.basename(NETWORK))[0]
elif NETWORK in _NETWORKS:
    net_path, NET_ID = _NETWORKS[NETWORK], NETWORK
else:
    print(f"Network '{NETWORK}' not found.")
    sys.exit(1)

SEP = "=" * 72
print(SEP)
print(f"Decomposition Theorem — {NET_ID}")
print(SEP)

print("\n[Stage 0] Loading network + running cot_gen pipeline...")
t0 = time.perf_counter()
rn_pycot = read_txt(net_path)
rn_data = build_rndata(rn_pycot, network_id=NET_ID)
ercs = compute_ercs(rn_data, verify=VERIFY)
hier = build_hierarchy(ercs)
syn = compute_synergies_basis_first(ercs, hier)
comp = compute_complementarities(ercs, hier, syn)
epm_result = compute_epms(rn_data, ercs, hier, syn, comp, verbose=False)
espm_result = compute_espm(rn_data, ercs, hier, syn, comp, epm_result,
                            max_order=ESPM_MAX_ORDER, verbose=False)
pipeline_ms = (time.perf_counter() - t0) * 1000
print(f"  ERCs={len(ercs)}  EPMs={len(epm_result.all_epm_masks)}  "
      f"total SOs={len(espm_result.all_so_masks)}  ({pipeline_ms:.0f} ms)")

so_order = {sp: 0 for sp in epm_result.all_epm_masks}
for order, masks in espm_result.espm_by_order.items():
    for sp in masks:
        so_order[sp] = order
so_lattice = build_so_lattice(espm_result.all_so_masks, so_order)

print("\n[Stage D] Decomposing every SO in the lattice (E / F / fragile circuits)...")
t0 = time.perf_counter()
S_full = build_full_stoich(rn_pycot)
results = decompose_hierarchy(so_lattice, rn_data, S_full, verbose=False)
decomp_ms = (time.perf_counter() - t0) * 1000
print(f"  Decomposed {len(results)} SOs in {decomp_ms:.0f} ms")

# ── Per-order evolution table ───────────────────────────────────────────────
print(f"\n{SEP}\nDecomposition structure vs. hierarchy order\n{SEP}")
header = f"{'order':>5} {'#SOs':>6} {'avg|E|':>7} {'avg|F|':>7} {'avg #circ':>10} {'avg circ size':>14} {'#organizations':>15}"
print(header)
print("-" * len(header))

by_order: dict[int, list] = {}
for sp, r in results.items():
    by_order.setdefault(so_lattice.order_of.get(sp, 0), []).append(r)

for order in sorted(by_order):
    rs = by_order[order]
    n = len(rs)
    avg_e = sum(bin(r.E_mask).count('1') for r in rs) / n
    avg_f = sum(bin(r.F_mask).count('1') for r in rs) / n
    avg_nc = sum(len(r.circuits) for r in rs) / n
    all_sizes = [c.size() for r in rs for c in r.circuits]
    avg_cs = sum(all_sizes) / len(all_sizes) if all_sizes else 0.0
    n_org = sum(1 for r in rs if r.is_organization)
    print(f"{order:>5} {n:>6} {avg_e:>7.1f} {avg_f:>7.1f} {avg_nc:>10.2f} "
          f"{avg_cs:>14.1f} {n_org:>8}/{n:<6}")

n_org_total = sum(1 for r in results.values() if r.is_organization)
print(f"\nTotal: {n_org_total}/{len(results)} SOs are full organizations "
      f"(every fragile circuit self-maintains).")

# ── The maximal SO's decomposition, in detail ───────────────────────────────
max_sp = max(results.keys(), key=lambda sp: bin(sp).count('1'))
r_max = results[max_sp]
print(f"\n{SEP}\nLargest SO reached (order {so_lattice.order_of.get(max_sp, 0)}, "
      f"{bin(max_sp).count('1')} species): {r_max.summary()}\n{SEP}")
for i, c in enumerate(sorted(r_max.circuits, key=lambda c: -c.size())):
    names = rn_data.bitset_to_names(c.species_mask)
    tag = "SELF-MAINTAINING" if c.is_self_maintaining else "FAILS self-maintenance"
    preview = ", ".join(names[:6]) + (", ..." if len(names) > 6 else "")
    print(f"  D{i+1}: {c.size()} species, {len(c.reaction_ids)} reactions -> {tag}")
    print(f"       [{preview}]")

# ── Stage V: visualizations ─────────────────────────────────────────────────
if MAKE_PLOTS:
    print(f"\n[Stage V] Rendering visualizations...")
    net_dir = os.path.join(OUTPUT_DIR, NET_ID)
    os.makedirs(net_dir, exist_ok=True)

    evo_path = plot_decomposition_evolution(
        results, so_lattice, os.path.join(net_dir, "decomposition_evolution.png"),
        title=f"{NET_ID} — decomposition structure vs. hierarchy order",
    )
    print(f"  wrote {evo_path}")

    # ── Organization-only views (is_organization == True subset) ────────────
    org_results = {sp: r for sp, r in results.items() if r.is_organization}
    print(f"\n  {len(org_results)}/{len(results)} SOs are full organizations "
          f"-- building organization-specific views...")

    if org_results:
        hasse_path = plot_organization_hasse(
            org_results, so_lattice, os.path.join(net_dir, "organization_hasse.png"),
            title=f"{NET_ID} — Hasse diagram of organizations",
        )
        print(f"  wrote {hasse_path}")

        chains_path = plot_organization_chains(
            org_results, so_lattice, os.path.join(net_dir, "organization_chains.png"),
            n_chains=5, title=f"{NET_ID} — organization chains",
        )
        print(f"  wrote {chains_path}")
    else:
        print("  (no full organizations found -- skipping organization_hasse / organization_chains)")

    print(f"\nAll deep-report files in: {os.path.relpath(net_dir, _repo)}")
