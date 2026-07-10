"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  run_network.py — Full COT generative-structure pipeline on one network     ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT IT DOES
------------
Runs all stages of the COT analysis on a single reaction network:

  Stage 0  Load the network file and build the internal RNData structure.
  Stage 1  Discover all ERCs: closures, minimal bases, req/prod masks.
  Stage H  Build the ERC containment hierarchy (Hasse diagram).
  Stage S  Compute basic → maximal → fundamental synergies.
  Stage C  Compute complementarities (basic, pure, fundamental per paper Defs 22–26).
  Stage G  Identify primitive ERCs and generative basis (synergy reachability).
  Stage E  Compute EPMs (minimal persistent modules, paper Def 29).

WHAT IT OUTPUTS
---------------
  • Per-stage counts and wall-clock time (printed to console).
  • MetaNetwork summary: ERCs, hierarchy, synergy, complementarity, generators, EPMs.
  • Up to 10 sample fundamental synergies (species names shown).
  • Up to 10 sample fundamental complementarities.

HOW TO RUN
----------
  1. Edit the CONFIGURATION section below.
  2. Press ▶ (play) in VS Code, or run:
       python projects/COT_fundamental_Generators/scripts/run_network.py
"""

# ── Path setup (do not edit) ──────────────────────────────────────────────────
from __future__ import annotations
import os, sys, time
_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# ── Imports (do not edit) ─────────────────────────────────────────────────────
from pyCOT.io.functions      import read_txt
from cot_gen.io_pyCOT        import build_rndata
from cot_gen.erc             import compute_ercs
from cot_gen.hierarchy       import build_hierarchy
from cot_gen.synergy         import compute_synergies, compute_synergies_basis_first
from cot_gen.complementarity import compute_complementarities
from cot_gen.generators      import compute_generators
from cot_gen.epm             import compute_epms, compute_espm
from cot_gen.metanetwork     import build_metanetwork
from cot_gen.metrics         import Counters

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — edit this block, then press Play                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── Network ───────────────────────────────────────────────────────────────────
# Short name from the database, e.g.: "e_coli_core", "iJO1366"
# Or a direct path to a .txt file.
# Set SHOW_LIST = True to print all available network names.
NETWORK   = "e_coli_core"
#NETWORK  = "BMID000000141754_url"
#NETWORK  = "BIOMD0000000185"
#NETWORK  = "iMM1415" $4000 reactions!
#NETWORK  = "iMM904" 
NETWORK  = "iND750" 
NETWORK  = "iNF517" 
NETWORK  = "iNJ661" 
SHOW_LIST = False

# ── Synergy algorithm selection ───────────────────────────────────────────────
# USE_ERC_HIERARCHY   : Fast basis-first method using the ERC containment
#                        hierarchy.  Always produces fundamental synergies.
#                        Recommended for all normal runs.
# USE_ERC_BRUTE_FORCE : Slow pair-enumeration over all ERC combinations.
#                        Respects SYNERGY_LEVEL below.  Use for validation.
# Both True           : runs both algorithms and compares results side-by-side.
USE_ERC_HIERARCHY   = True
USE_ERC_BRUTE_FORCE = False

# ── Synergy level (only applies when USE_ERC_BRUTE_FORCE = True) ──────────────
# "basic" | "maximal" | "fundamental"
SYNERGY_LEVEL = "fundamental"

# ── ERC count safety limit (synergy is O(|E|^3)) ─────────────────────────────
# Set to 0 to always run synergy.
ERC_LIMIT = 2500

# ── Compute EPMs ──────────────────────────────────────────────────────────────
COMPUTE_EPM = True

# ── Compute ESPMs (requires COMPUTE_EPM = True and synergy/comp available) ────
COMPUTE_ESPM = True
ESPM_MAX_ORDER = 0   # stop extending after this order

# ── Verification ─────────────────────────────────────────────────────────────
VERIFY = True

# ── Results CSV (set to None to disable) ─────────────────────────────────────
RESULTS_CSV = os.path.join(_here, "..", "outputs", "cot_results.csv")

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Script body — no need to edit below this line                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

_DATA_ROOT = os.path.join(_repo, "data", "biomodels")

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

if SHOW_LIST:
    items = sorted(_NETWORKS.items(), key=lambda kv: os.path.basename(kv[1]))
    print(f"{'Name':<40}  {'Folder'}")
    print("-" * 70)
    for name, path in items:
        folder = os.path.basename(os.path.dirname(path))
        print(f"  {name:<38}  {folder}")
    print(f"\n{len(items)} networks available.")
    sys.exit(0)

# ── Resolve network path ──────────────────────────────────────────────────────
if os.path.isfile(NETWORK):
    net_path = os.path.abspath(NETWORK)
    NET_ID   = os.path.splitext(os.path.basename(net_path))[0]
elif NETWORK in _NETWORKS:
    net_path = _NETWORKS[NETWORK]
    NET_ID   = NETWORK
else:
    lower   = NETWORK.lower()
    matches = [(k, v) for k, v in _NETWORKS.items() if k.lower().startswith(lower)]
    if len(matches) == 1:
        NET_ID, net_path = matches[0]
        print(f"[info] Matched '{NET_ID}'")
    elif len(matches) > 1:
        print(f"Ambiguous network name '{NETWORK}'. Matches:")
        for k, _ in matches[:10]:
            print(f"  {k}")
        sys.exit(1)
    else:
        print(f"Network '{NETWORK}' not found. Set SHOW_LIST = True to see options.")
        sys.exit(1)

# ── Banner ────────────────────────────────────────────────────────────────────
SEP = "=" * 72
print(SEP)
print(f"COT Generative Analysis  —  {NET_ID}")
print(f"File     : {os.path.relpath(net_path, _repo)}")
_syn_methods = "+".join(m for m, f in [("hierarchy", USE_ERC_HIERARCHY), ("brute-force", USE_ERC_BRUTE_FORCE)] if f) or "none"
print(f"Synergy  : {_syn_methods}   Verify: {VERIFY}")
print(SEP)

# ── Stage 0: Load ─────────────────────────────────────────────────────────────
print("\n[Stage 0] Loading network...")
t0 = time.perf_counter()
try:
    rn_pycot = read_txt(net_path)
except Exception as e:
    print(f"ERROR loading {net_path}: {e}")
    sys.exit(1)
load_ms = (time.perf_counter() - t0) * 1000

t0 = time.perf_counter()
rn = build_rndata(rn_pycot, network_id=NET_ID)
build_ms = (time.perf_counter() - t0) * 1000

n_E0 = bin(rn.E0_mask).count('1')
print(f"  Species       : {rn.n_species}")
print(f"  Reactions     : {rn.n_reactions}")
print(f"  E0 (inflow)   : {n_E0} species")
print(f"  Load time     : {load_ms:.1f} ms  /  Build time: {build_ms:.1f} ms")

# ── Stage 1: ERCs ─────────────────────────────────────────────────────────────
note = "  (assertions OFF)" if not VERIFY else ""
print(f"\n[Stage 1] Computing ERCs...{note}")
t0 = time.perf_counter()
ctr1 = Counters()
ercs = compute_ercs(rn, counters=ctr1, verify=VERIFY)
s1_ms = (time.perf_counter() - t0) * 1000

n_pers = sum(1 for e in ercs if e.is_persistent())
print(f"  ERCs found    : {len(ercs)}")
print(f"  Persistent    : {n_pers}")
print(f"  Min-bases     : {ctr1.get('erc.min_bases_total', 0)}")
print(f"  Wall time     : {s1_ms:.1f} ms")

if not ercs:
    print("\nNo ERCs found.")
    sys.exit(0)

# ── Stage H: Hierarchy ────────────────────────────────────────────────────────
print(f"\n[Stage H] Building ERC containment hierarchy...")
t0 = time.perf_counter()
hier = build_hierarchy(ercs)
sh_ms = (time.perf_counter() - t0) * 1000

n_pairs    = len(ercs) * (len(ercs) - 1) // 2
n_comp_p   = sum(len(hier.ancestors[i]) for i in range(hier.n))
hasse_e    = sum(len(hier.parents[i]) for i in range(hier.n))
print(f"  ERC pairs     : {n_pairs}  (comparable: {n_comp_p}, incomparable: {n_pairs - n_comp_p})")
print(f"  Hasse edges   : {hasse_e}")
print(f"  Wall time     : {sh_ms:.1f} ms")

# ── Stage S: Synergies ────────────────────────────────────────────────────────
_n_ercs       = len(ercs)
_skip_synergy = not USE_ERC_HIERARCHY and not USE_ERC_BRUTE_FORCE
ss_ms         = 0.0
ss_bf_ms      = 0.0
syn = syn_fund = None

if _skip_synergy:
    print(f"\n[Stage S] SKIPPED — neither USE_ERC_HIERARCHY nor USE_ERC_BRUTE_FORCE selected.")
elif ERC_LIMIT > 0 and _n_ercs > ERC_LIMIT:
    est_s = _n_ercs ** 3 / 5e6
    print(f"\n[Stage S] SKIPPED — {_n_ercs} ERCs exceeds ERC_LIMIT={ERC_LIMIT}.")
    print(f"  (~{est_s/60:.0f} min at 5M ops/s).  Set ERC_LIMIT = 0 to force.")
    _skip_synergy = True

if not _skip_synergy:
    # ── Hierarchy (fast) method ───────────────────────────────────────────────
    if USE_ERC_HIERARCHY:
        lbl = "ERC-hierarchy" + ("  [primary]" if USE_ERC_BRUTE_FORCE else "")
        print(f"\n[Stage S] Computing synergies ({lbl})...")
        t0       = time.perf_counter()
        syn_fund = compute_synergies_basis_first(ercs, hier)
        ss_ms    = (time.perf_counter() - t0) * 1000
        syn      = syn_fund
        print(f"  Basic         : {len(syn_fund.basic)}")
        print(f"  Maximal       : {len(syn_fund.maximal)}")
        print(f"  Fundamental   : {len(syn_fund.fundamental)}")
        print(f"  Wall time     : {ss_ms:.1f} ms")

    # ── Brute-force method ────────────────────────────────────────────────────
    if USE_ERC_BRUTE_FORCE:
        if USE_ERC_HIERARCHY:
            # Both selected: run brute-force as comparison against hierarchy
            print(f"\n[Stage S-BF] Brute-force comparison (level={SYNERGY_LEVEL})...")
            t0     = time.perf_counter()
            syn_pf = compute_synergies(ercs, hier, level="fundamental", direct=True)
            ss_bf_ms = (time.perf_counter() - t0) * 1000

            set_hf = {(s.i, s.j, s.k) for s in syn_fund.fundamental}
            set_pf = {(s.i, s.j, s.k) for s in syn_pf.fundamental}
            match  = set_hf == set_pf
            print(f"  Basic         : {len(syn_pf.basic)}"
                  f"  (hierarchy: {len(syn_fund.basic)})")
            print(f"  Maximal       : {len(syn_pf.maximal)}"
                  f"  (hierarchy: {len(syn_fund.maximal)})")
            print(f"  Fundamental   : {len(syn_pf.fundamental)}"
                  f"  (hierarchy: {len(syn_fund.fundamental)})")
            print(f"  Wall time     : {ss_bf_ms:.1f} ms  (hierarchy: {ss_ms:.1f} ms)")
            print(f"  Sets match    : {'YES' if match else 'NO — MISMATCH!'}")
            if not match:
                only_hf = set_hf - set_pf
                only_pf = set_pf - set_hf
                if only_hf:
                    print(f"  Only in hierarchy   ({len(only_hf)}): {list(only_hf)[:5]}")
                if only_pf:
                    print(f"  Only in brute-force ({len(only_pf)}): {list(only_pf)[:5]}")
        else:
            # Brute-force only
            print(f"\n[Stage S] Computing {SYNERGY_LEVEL} synergies (brute-force)...")
            t0   = time.perf_counter()
            ctrS = Counters()
            syn  = compute_synergies(ercs, hier, level=SYNERGY_LEVEL, direct=True, counters=ctrS)
            ss_ms = (time.perf_counter() - t0) * 1000
            print(f"  Basic         : {len(syn.basic)}")
            if SYNERGY_LEVEL in ("maximal", "fundamental"):
                print(f"  Maximal       : {len(syn.maximal)}")
            if SYNERGY_LEVEL == "fundamental":
                print(f"  Fundamental   : {len(syn.fundamental)}")
            print(f"  Wall time     : {ss_ms:.1f} ms")

            if SYNERGY_LEVEL != "fundamental":
                print(f"\n[Stage S-F] Recomputing at fundamental level (needed downstream)...")
                t0 = time.perf_counter()
                syn_fund = compute_synergies(ercs, hier, level="fundamental", direct=True)
                print(f"  Fundamental   : {len(syn_fund.fundamental)}  ({(time.perf_counter()-t0)*1000:.1f} ms)")
            else:
                syn_fund = syn

# ── Stage C: Complementarity ──────────────────────────────────────────────────
print(f"\n[Stage C] Computing complementarities (paper Defs 22–26)...")
t0   = time.perf_counter()
comp = compute_complementarities(ercs, hier, syn_result=syn_fund)
sc_ms = (time.perf_counter() - t0) * 1000

print(f"  Basic         : {len(comp.basic)}")
print(f"  Pure          : {len(comp.pure)}  (basic and not synergetic)")
print(f"  Fundamental   : {len(comp.fundamental)}")
print(f"  Wall time     : {sc_ms:.1f} ms")

# ── Stage G: Generators ───────────────────────────────────────────────────────
sg_ms = 0.0
gen   = None
if syn_fund is not None:
    print(f"\n[Stage G] Computing primitive ERCs and generative basis...")
    t0    = time.perf_counter()
    gen   = compute_generators(ercs, hier, syn_fund)
    sg_ms = (time.perf_counter() - t0) * 1000
    print(f"  Primitives    : {len(gen.primitive_indices)} of {len(ercs)} ERCs")
    print(f"  Basis reach   : {len(gen.basis_reach)} ERCs")
    print(f"  Coverage      : {gen.coverage:.1%}  (complete={gen.is_complete()})")
    print(f"  Wall time     : {sg_ms:.1f} ms")
else:
    print(f"\n[Stage G] Skipped (synergy not computed).")

# ── Stage E: EPMs ─────────────────────────────────────────────────────────────
se_ms = 0.0
epm   = None
if COMPUTE_EPM:
    print(f"\n[Stage E] Computing EPMs (paper Def 29, adjacency traversal)...")
    t0  = time.perf_counter()
    epm = compute_epms(rn, ercs, hier, syn_result=syn_fund, comp_result=comp, verbose=True)
    se_ms = (time.perf_counter() - t0) * 1000

    st1 = epm.stats
    print(f"  Single-ERC EPMs  : {len(epm.single_epm_indices)}  (minimal persistent ERCs)")
    print(f"  Multi-ERC EPMs   : {len(epm.multi_epm_masks)}  (from adjacency traversal)")
    print(f"  Total EPMs       : {len(epm.all_epm_masks)}")
    print(f"  Leaves (non-SSM) : {len(epm.leaf_masks)}")
    print(f"  States explored  : {st1.get('states_explored', 0)}")
    print(f"  Comp extensions  : {st1.get('comp_extensions', 0)}")
    print(f"  Syn  extensions  : {st1.get('syn_extensions', 0)}")
    print(f"  Wall time        : {se_ms:.1f} ms")
else:
    print(f"\n[Stage E] Skipped (COMPUTE_EPM = False).")

# ── Stage ES: ESPMs ───────────────────────────────────────────────────────────
ses_ms = 0.0
espm   = None
if COMPUTE_ESPM and epm is not None and syn_fund is not None:
    print(f"\n[Stage ES] Computing ESPMs (paper Def 31, Mode-2 extension)...")
    print(f"  EPMs (order 0)   : {len(epm.all_epm_masks)}")
    t0   = time.perf_counter()
    espm = compute_espm(rn, ercs, hier, syn_fund, comp, epm,
                        max_order=ESPM_MAX_ORDER, verbose=True)
    ses_ms = (time.perf_counter() - t0) * 1000

    total_espm = espm.total_espm()
    max_ord    = espm.max_order()
    print(f"  --- order totals ---")
    for ord_k, masks in sorted(espm.espm_by_order.items()):
        lk = espm.leaf_masks_by_order.get(ord_k, [])
        print(f"  ESPMs order {ord_k:<3}  : {len(masks):5d}  (leaves: {len(lk)})")
    print(f"  Total ESPMs      : {total_espm}  (max order: {max_ord})")
    print(f"  All SOs          : {len(espm.all_so_masks)}")
    print(f"  Wall time        : {ses_ms:.1f} ms")
elif COMPUTE_ESPM:
    print(f"\n[Stage ES] Skipped (EPM or synergy not available).")

# ── Sample fundamental synergies ──────────────────────────────────────────────
if syn_fund is not None and syn_fund.fundamental:
    sample = syn_fund.fundamental[:10]
    print(f"\n{SEP}")
    print(f"SAMPLE FUNDAMENTAL SYNERGIES (first {len(sample)} of {len(syn_fund.fundamental)})")
    print(SEP)
    for s in sample:
        ei = ercs[s.i]; ej = ercs[s.j]; ek = ercs[s.k]
        print(f"  E{s.i}(sz={ei.size()}, req={bin(ei.req_mask).count('1')}) + "
              f"E{s.j}(sz={ej.size()}, req={bin(ej.req_mask).count('1')})"
              f"  ->  E{s.k}(sz={ek.size()}, req={bin(ek.req_mask).count('1')})")

# ── Sample fundamental complementarities ─────────────────────────────────────
if comp.fundamental:
    sample = comp.fundamental[:10]
    print(f"\n{SEP}")
    print(f"SAMPLE FUNDAMENTAL COMPLEMENTARITIES (first {len(sample)} of {len(comp.fundamental)})")
    print(SEP)
    for fc in sample:
        ep = ercs[fc.prod_idx]; ec = ercs[fc.cons_idx]
        sp_name = rn.species_name(fc.species)
        print(f"  E{fc.prod_idx}(sz={ep.size()}) --[{sp_name}]--> E{fc.cons_idx}(sz={ec.size()})")

# ── MetaNetwork summary (ESPM line appended) ──────────────────────────────────
mn = build_metanetwork(ercs, hier, syn_fund, comp, gen, epm)
print(f"\n{SEP}")
print("METANETWORK SUMMARY")
print(SEP)
mn.print_summary(prefix="  ")
if espm is not None:
    total_espm = espm.total_espm()
    order_str  = "  ".join(f"o{k}:{len(v)}" for k, v in sorted(espm.espm_by_order.items()))
    print(f"  ESPMs:             {total_espm} total  ({order_str or 'none'})")
    print(f"  All SOs:           {len(espm.all_so_masks)}")

# ── Timing summary ────────────────────────────────────────────────────────────
total_ms = s1_ms + sh_ms + ss_ms + ss_bf_ms + sc_ms + sg_ms + se_ms + ses_ms
_s_timing = (f"S-hier {ss_ms:.0f}ms | S-brute {ss_bf_ms:.0f}ms"
             if USE_ERC_HIERARCHY and USE_ERC_BRUTE_FORCE
             else f"S {ss_ms:.0f}ms")
print(f"\n{SEP}")
print(f"Total pipeline: {total_ms/1000:.2f}s  "
      f"(ERC {s1_ms:.0f}ms | H {sh_ms:.0f}ms | "
      f"{_s_timing} | C {sc_ms:.0f}ms | G {sg_ms:.0f}ms | "
      f"E {se_ms:.0f}ms | ES {ses_ms:.0f}ms)")

# ── Save to CSV ───────────────────────────────────────────────────────────────
if RESULTS_CSV:
    from cot_gen.results_io import make_row, update_results_csv
    _timing_ms = {
        "load": load_ms + build_ms,
        "erc":  s1_ms,
        "hier": sh_ms,
        "syn":  ss_ms,
        "comp": sc_ms,
        "gen":  sg_ms,
        "epm":  se_ms,
        "espm": ses_ms,
    }
    _row = make_row(NET_ID, rn, ercs, syn_fund, comp, epm, espm, _timing_ms)
    update_results_csv(RESULTS_CSV, _row)
    print(f"\nResults saved → {os.path.relpath(RESULTS_CSV, _repo)}")
