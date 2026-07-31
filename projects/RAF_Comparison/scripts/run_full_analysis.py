"""
run_full_analysis.py -- ONE script, RAF structure + COT decomposition
structure, for one network, end to end, with PNG plots. Edit the
CONFIGURATION block below and press Play; everything else is generic and
should run unmodified on any other network in the corpus.

WHAT IT COMPUTES
----------------
RAF layer:
  - CRS induction (two modes -- pick whichever matches the network):
      'net_zero'       : for kinetic/BIOMD-style networks with SBML-declared
                         modifiers (catalyst = species with net-zero effect
                         within one reaction). See raf/biomodel_crs.py.
      'cofactor_pools' : for BiGG/flux-balance-style networks, where the
                         above heuristic structurally cannot find anything
                         (cofactor pairs like nad_c/nadh_c are never the
                         same species token on both sides of one reaction).
                         See raf/cofactor_pools.py.
  - Global maxRAF, gen(F0,maxRAF).

COT decomposition layer:
  - The FULL semi-organization lattice (cot_gen's own ERC/hierarchy/
    synergy/complementarity/EPM/ESPM pipeline -- no catalysis notion at
    all), decomposed order by order (E/F/fragile circuits/is_organization
    for every semi-organization, via Decomposition_Theorem's engine).

Cross-comparison, PER ORGANIZATION (not just the largest one):
  - For every organization X in the hierarchy, the reaction set actually
    triggered (R_X) is used to build a RESTRICTED CRS (same food, same
    catalyst assignment, reaction universe cut down to R_X), and maxRAF is
    recomputed WITHIN that restriction. This answers, at every level of the
    hierarchy separately: if a RAF-based analysis could only ever see the
    reactions this particular organization actually uses, how much of it
    would it reach, and what fraction of R_X even has an identifiable
    catalyst? Comparing this across all organizations (not just the top
    one) is the point -- the RAF/COT gap need not be constant across
    scales.

PLOTS (PNG, written to OUTPUT_DIR)
----------------------------------
  1. <net>_organization_hasse.png -- the FULL organization Hasse diagram
     (decomp/viz.py's plot_organization_hasse, reused as-is): every
     organization drawn as its own horizontal E/F/circuit stacked bar,
     positioned by hierarchy order, with containment arrows -- shows
     multiple PARALLEL circuits explicitly, never averaged away.
  2. <net>_organization_chains.png -- the same diagram restricted to a
     handful of maximally-divergent root-to-top chains (less cluttered
     reading of the same structure on larger lattices).
  3. <net>_raf_vs_organization_per_org.png -- for EVERY organization: |X|
     vs. the LOCAL gen(F0,maxRAF) reachable using only that organization's
     own R_X, plus the catalyzed/uncatalyzed reaction split, both ordered
     by hierarchy level.
  4. <net>_top_organization_breakdown.png -- E / F / per-circuit species
     counts for the largest organization found, one fixed color per role.

HOW TO RUN
----------
  1. Edit the CONFIGURATION section below.
  2. python projects/RAF_Comparison/scripts/run_full_analysis.py
"""
from __future__ import annotations

import os
import sys
import time
import tempfile

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_decomp_proj = os.path.normpath(os.path.join(_here, "..", "..", "Decomposition_Theorem"))
_cot_gen_proj = os.path.normpath(os.path.join(_here, "..", "..", "COT_fundamental_Generators"))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, _decomp_proj, _cot_gen_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8", line_buffering=True)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION -- edit this block, then press Play                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# NETWORK_CHOICE selects EXACTLY ONE of the two networks reported in
# Sec. 4 of the paper -- this replaces the earlier pattern of stacking
# several assignments where a later one silently overrides an earlier one
# (easy to lose track of which config is actually active; this is how a
# stale FOOD=[] from a previous edit went unnoticed). Add new networks as
# new elif branches, never by re-assigning NETWORK_PATH/CATALYST_MODE/
# FOOD/OUTFLOW a second time.
NETWORK_CHOICE = "heteropolymer"   # "heteropolymer" (Sec. 4.1) | "ecoli_none" (Sec. 4.2)

# FOOD: leave as None (not []) to use a network's OWN native inflow
# reactions unmodified. FOOD and None are NOT interchangeable: setting
# FOOD to ANY list, including an empty one, STRIPS every native inflow
# reaction from the network and replaces it with one inflow per token in
# the list -- FOOD=[] therefore means "replace native food with nothing"
# (a fully food-less, closed network), not "leave native food alone".
# OUTFLOW is different: it only ever ADDS outflow reactions on top of
# whatever's already there (Def:translation's Omega), so OUTFLOW=[] really
# does mean "no outflow added", safely, every time.

if NETWORK_CHOICE == "heteropolymer":
    # Sec. 4.1's network: catalysts are written literally as net-zero
    # species on both sides of each reaction line. Regenerate this file
    # with scripts/generate_heteropolymer_network.py if it's missing.
    NETWORK_PATH = os.path.join(_proj, "data", "heteropolymer", "heteropolymer_L3_uniform_cleave0.5.txt")
    CATALYST_MODE = "net_zero"
    FOOD: list[str] | None = ["0", "1"]
    OUTFLOW: list[str] = []

elif NETWORK_CHOICE == "ecoli_none":
    # Sec. 4.2's network: no catalyst induction at all (C(r)=empty for
    # every r, taking the data at face value -- isolates the COT
    # decomposition layer, which never reads `catalysts`, from any
    # catalyst-induction guess).
    NETWORK_PATH = os.path.join(_proj, "data", "bigg", "bigg_e_coli_core.txt")
    CATALYST_MODE = "none"
    # bigg_e_coli_core.txt's native inflows are 7 species, not 6 -- it
    # includes co2_e (R74_rev) alongside the six aerobic-glucose-medium
    # species. FOOD must be set explicitly to exclude CO2 and match the
    # paper's scenario.
    FOOD = ["glc__D_e", "o2_e", "nh4_e", "pi_e", "h2o_e", "h_e"]
    OUTFLOW = []

else:
    raise ValueError(f"unknown NETWORK_CHOICE: {NETWORK_CHOICE!r}")

ESPM_MAX_ORDER = 60          # cap on how large a semi-organization ESPM search explores
OUTPUT_DIR = os.path.join(_here, "..", "figures")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Script body -- no need to edit below this line                            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pyCOT.io.functions import read_txt  # noqa: E402
from cot_gen.io_pyCOT import build_rndata  # noqa: E402
from cot_gen.erc import compute_ercs  # noqa: E402
from cot_gen.hierarchy import build_hierarchy  # noqa: E402
from cot_gen.synergy import compute_synergies_basis_first  # noqa: E402
from cot_gen.complementarity import compute_complementarities  # noqa: E402
from cot_gen.epm import compute_epms, compute_espm  # noqa: E402
from cot_gen.deep_report import build_so_lattice  # noqa: E402

from decomp.bridge import build_full_stoich  # noqa: E402
from decomp.hierarchy import decompose_hierarchy  # noqa: E402
from decomp.viz import plot_organization_hasse, plot_organization_chains  # noqa: E402

from raf.biomodel_crs import crs_from_biomodel, crs_from_biomodel_cofactor_pools, crs_from_biomodel_no_catalysis  # noqa: E402
from raf.raf_algo import compute_maxRAF  # noqa: E402
from raf.crs import gen, CRS  # noqa: E402
from raf.inflow_scenarios import build_scenario_txt  # noqa: E402

# Fixed categorical colors, one per structural role -- assigned by IDENTITY
# (E is always this color, F is always that color, circuit i is always the
# i-th color in the sequence) and never reassigned/cycled based on which
# roles happen to be present in a given plot.
COLOR_E = "#8a8f98"       # enabler -- neutral gray (rare, background role)
COLOR_F = "#3b82c4"       # overproducible / food-like -- blue
CIRCUIT_COLORS = ["#c4553b", "#c49b3b", "#5ba86b", "#8a5bc4", "#3bc4b8", "#c45ba0"]

# NET_DIR is OUTPUT_DIR itself (flat, no per-network subfolder): the
# figure filenames below already embed NET_ID, and OUTPUT_DIR is the
# figures/ folder referenced directly (e.g. figures/bigg_e_coli_core_
# organization_hasse) by the paper's own \includegraphics paths.
NET_ID = os.path.splitext(os.path.basename(NETWORK_PATH))[0]
NET_DIR = OUTPUT_DIR
os.makedirs(NET_DIR, exist_ok=True)

SEP = "=" * 100
print(f"{SEP}\n{NET_ID}  (catalyst_mode={CATALYST_MODE})\n{SEP}")

base_txt = open(NETWORK_PATH, encoding="utf-8", errors="ignore").read()
if FOOD is not None or OUTFLOW:
    scenario_txt = build_scenario_txt(base_txt, FOOD or [], OUTFLOW)
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as f:
        f.write(scenario_txt)
        _tmp_path = f.name
    rn_pycot = read_txt(_tmp_path, exact_names=True)
    os.unlink(_tmp_path)
else:
    rn_pycot = read_txt(NETWORK_PATH, exact_names=True)

rn_data = build_rndata(rn_pycot, network_id=NET_ID)
print(f"species={rn_data.n_species}  reactions={rn_data.n_reactions}  "
      f"food(E0 closure)={bin(rn_data.E0_mask).count('1')}")
if rn_data.E0_mask == 0:
    print(f"  !! WARNING: food set is EMPTY after applying FOOD={FOOD!r}. If FOOD is an "
          f"empty list (not None), this is almost certainly a mistake -- an empty list "
          f"STRIPS every native inflow reaction and adds none back (see FOOD/OUTFLOW "
          f"comment above), unlike FOOD=None which leaves native inflows untouched. "
          f"Everything downstream (RAF layer AND semi-organization lattice) will be "
          f"degenerate with no food to bootstrap from.")

# ── RAF layer (global) ──────────────────────────────────────────────────────
print(f"\n[RAF layer -- {CATALYST_MODE}]")
if CATALYST_MODE == "cofactor_pools":
    crs, crs_stats = crs_from_biomodel_cofactor_pools(rn_pycot, rn_data)
elif CATALYST_MODE == "none":
    crs, crs_stats = crs_from_biomodel_no_catalysis(rn_pycot, rn_data)
else:
    crs, crs_stats = crs_from_biomodel(rn_pycot, rn_data)
if FOOD is not None:
    crs.food = frozenset(FOOD) & crs.species
print(f"  {crs_stats}")
maxraf = compute_maxRAF(crs)
X_raf = gen(crs.food, [crs.reactions[n] for n in maxraf]) if maxraf else frozenset(crs.food)
print(f"  global maxRAF={len(maxraf)} reactions  gen(F0,maxRAF)={len(X_raf)} species")

# ── COT decomposition layer: full semi-organization hierarchy ──────────────
print(f"\n[COT layer -- full semi-organization hierarchy]")
t0 = time.perf_counter()
ercs = compute_ercs(rn_data, verify=False)
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
n_org = sum(1 for r in results.values() if r.is_organization)
print(f"  ERCs={len(ercs)}  SOs={len(results)}  organizations={n_org}  "
      f"({time.perf_counter()-t0:.1f}s)")

# ── Per-organization RAF comparison (NOT just the largest one) ─────────────
def local_maxraf(r) -> tuple[list[str], frozenset[str]]:
    """maxRAF computed using ONLY the reactions this organization's own R_X
    triggers -- same food, same catalyst assignment, smaller reaction
    universe -- rather than the network-wide maxRAF."""
    r_x_names = {rn_data.reaction_name(idx) for idx in r.R_X}
    restricted = {name: rxn for name, rxn in crs.reactions.items() if name in r_x_names}
    restricted_crs = CRS(species=crs.species, reactions=restricted, food=crs.food)
    local = compute_maxRAF(restricted_crs)
    local_gen = gen(restricted_crs.food, [restricted_crs.reactions[n] for n in local]) if local \
        else frozenset(restricted_crs.food)
    return local, local_gen


orgs = [(sp, r) for sp, r in results.items() if r.is_organization]
orgs.sort(key=lambda kv: bin(kv[0] | kv[1].E_mask | kv[1].F_mask).count('1'))

print(f"\n[Per-organization RAF vs COT comparison]")
header = (f"{'order':>5} {'|X|':>5} {'|R_X|':>6} {'catalyzed':>9} {'uncatalyzed':>11} "
          f"{'local maxRAF rx':>15} {'local gen(F0,maxRAF)':>20}   |E| |F|  circuit sizes")
print(header)
per_org_rows = []
for sp, r in orgs:
    n_species = bin(sp | r.E_mask | r.F_mask).count('1')
    n_cat = n_uncat = 0
    for r_idx in r.R_X:
        rname = rn_data.reaction_name(r_idx)
        if rname not in crs.reactions:
            continue
        if crs.reactions[rname].catalysts:
            n_cat += 1
        else:
            n_uncat += 1
    local, local_gen = local_maxraf(r)
    n_e, n_f = bin(r.E_mask).count("1"), bin(r.F_mask).count("1")
    circuit_sizes = sorted((c.size() for c in r.circuits), reverse=True)
    print(f"{so_lattice.order_of.get(sp,0):>5} {n_species:>5} {len(r.R_X):>6} "
          f"{n_cat:>9} {n_uncat:>11} {len(local):>15} {len(local_gen):>20}   "
          f"{n_e:>3} {n_f:>3}  {circuit_sizes}")
    per_org_rows.append(dict(order=so_lattice.order_of.get(sp, 0), n_species=n_species,
                              n_cat=n_cat, n_uncat=n_uncat,
                              local_maxraf=len(local), local_gen=len(local_gen)))

top_sp, top_r = orgs[-1] if orgs else (None, None)
top_n_species = per_org_rows[-1]["n_species"] if per_org_rows else 0
if top_r is not None:
    X_raf_mask = 0
    sp_index = {rn_data.species_name(i): i for i in range(rn_data.n_species)}
    for name in X_raf:
        if name in sp_index:
            X_raf_mask |= 1 << sp_index[name]
    is_valid_so = X_raf_mask in results
    print(f"\nglobal gen(F0,maxRAF) ({len(X_raf)} species) is itself a real semi-organization: {is_valid_so}")

# ═══════════════════════════════════════════════════════════════════════════
# PLOTS
# ═══════════════════════════════════════════════════════════════════════════

org_results = {sp: r for sp, r in results.items() if r.is_organization}

# Plot 1 & 2: the actual organization Hasse diagram (every organization as
# its own E/F/circuit stacked bar -- parallel circuits are never averaged
# away here) and the divergent-chains reading of the same structure.
p1 = plot_organization_hasse(org_results, so_lattice,
                              os.path.join(NET_DIR, f"{NET_ID}_organization_hasse.png"),
                              title=f"{NET_ID} -- organization Hasse diagram")
if p1:
    print(f"\nwrote {p1}")
p2 = plot_organization_chains(org_results, so_lattice,
                               os.path.join(NET_DIR, f"{NET_ID}_organization_chains.png"),
                               title=f"{NET_ID} -- organization chains")
if p2:
    print(f"wrote {p2}")

# Plot 3: per-organization RAF-vs-COT comparison, ordered by hierarchy level.
if per_org_rows:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    idx = list(range(len(per_org_rows)))
    xticks = [f"ord.{row['order']}\n|X|={row['n_species']}" for row in per_org_rows]

    width = 0.35
    ax1.bar([i - width / 2 for i in idx], [row["n_species"] for row in per_org_rows],
            width=width, label="organization |X|", color=COLOR_F)
    ax1.bar([i + width / 2 for i in idx], [row["local_gen"] for row in per_org_rows],
            width=width, label="local gen(F0,maxRAF)", color=CIRCUIT_COLORS[0])
    ax1.set_xticks(idx)
    ax1.set_xticklabels(xticks, fontsize=7, rotation=0)
    ax1.set_ylabel("species count")
    ax1.set_title("RAF reach vs. true organization,\nper organization")
    ax1.legend(frameon=False, fontsize=8)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax2.bar(idx, [row["n_cat"] for row in per_org_rows], width=0.6,
            label="catalyzed", color=CIRCUIT_COLORS[0])
    ax2.bar(idx, [row["n_uncat"] for row in per_org_rows], width=0.6,
            bottom=[row["n_cat"] for row in per_org_rows], label="uncatalyzed", color=COLOR_E)
    ax2.set_xticks(idx)
    ax2.set_xticklabels(xticks, fontsize=7, rotation=0)
    ax2.set_ylabel("reaction count (R_X)")
    ax2.set_title("RAF-catalyst coverage of R_X,\nper organization")
    ax2.legend(frameon=False, fontsize=8)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.suptitle(f"{NET_ID}", fontsize=11)
    fig.tight_layout()
    p3 = os.path.join(NET_DIR, f"{NET_ID}_raf_vs_organization_per_org.png")
    fig.savefig(p3, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {p3}")

# Plot 4: E/F/circuit breakdown of the largest organization.
if top_r is not None:
    fig, ax = plt.subplots(figsize=(6, 4))
    labels4 = ["E\n(enablers)", "F\n(overproducible)"] + [f"D{i}" for i in range(len(top_r.circuits))]
    vals4 = [bin(top_r.E_mask).count('1'), bin(top_r.F_mask).count('1')] + \
            [c.size() for c in top_r.circuits]
    colors4 = [COLOR_E, COLOR_F] + CIRCUIT_COLORS[:len(top_r.circuits)]
    ax.bar(labels4, vals4, color=colors4)
    ax.set_ylabel("species count")
    ax.set_title(f"{NET_ID} -- largest organization ({top_n_species} species) breakdown")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    p4 = os.path.join(NET_DIR, f"{NET_ID}_top_organization_breakdown.png")
    fig.savefig(p4, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {p4}")

print(f"\nAll outputs in: {os.path.relpath(NET_DIR, _repo)}")
