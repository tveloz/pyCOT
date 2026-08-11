"""
bigg_synergy_explore.py — Test ERC synergy/complementarity on large BiGG networks.

Usage
-----
Edit SELECTED (~line 45) to choose a network, then run:

    python benchmarks/bigg_synergy_explore.py

Or pass arguments on the command line:

    python benchmarks/bigg_synergy_explore.py e_coli_core
    python benchmarks/bigg_synergy_explore.py iAF1260 --synergy fundamental
    python benchmarks/bigg_synergy_explore.py iAF1260 --no-verify
    python benchmarks/bigg_synergy_explore.py Recon3D --no-verify --synergy basic

Stages
------
Stage 0 : Load network file → RNData (cot_gen bridge)
Stage 1 : Compute ERCs (cot_gen bitset, Horn propagation)
            --no-verify skips per-ERC closure assertions (~50% faster on large nets)
Stage H : Build ERC containment Hasse diagram (cot_gen bitset, O(|E|^2))
Stage S : Compute synergies (cot_gen bitset: basic / maximal / fundamental)
Stage C : Compute complementarities (cot_gen bitset: Types 1-3)
Stage G : Compute generators (primitive ERCs + reachability)

Network catalogue
-----------------
e_coli_core   95 r,   72 s   — quick smoke test
iAF692       693 r,  405 s   — small metabolic
iAF987      1646 r, 1117 s   — medium
iAM_Pb448   1554 r,  907 s   — medium
iAF1260     2957 r, 1669 s   — large (good benchmark)
iAF1260b    2966 r, 1669 s
iJN1463     3717 r, 2150 s
iYS1720     4006 r, 2445 s
RECON1      5290 r, 2766 s
iCHOv1      9289 r, 4456 s
Recon3D    15834 r, 5835 s   — very large
"""
from __future__ import annotations

import os, sys, time

# ── path setup ────────────────────────────────────────────────────────────────
_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# ══════════════════════════════════════════════════════════════════════════════
#  SELECTION — edit these two lines, or pass CLI args
# ══════════════════════════════════════════════════════════════════════════════
SELECTED     = "iAF1260"     # key from catalogue above
SYNERGY_TYPE = "maximal"     # "basic" | "maximal" | "fundamental"
# ══════════════════════════════════════════════════════════════════════════════

CANDIDATES = {
    "e_coli_core" : ("bigg_e_coli_core.txt",    95,    72),
    "iAF692"      : ("bigg_iAF692.txt",         693,   405),
    "iAF987"      : ("bigg_iAF987.txt",        1646,  1117),
    "iAM_Pb448"   : ("bigg_iAM_Pb448.txt",     1554,   907),
    "iAF1260"     : ("bigg_iAF1260.txt",       2957,  1669),
    "iAF1260b"    : ("bigg_iAF1260b.txt",      2966,  1669),
    "iJN1463"     : ("bigg_iJN1463.txt",       3717,  2150),
    "iYS1720"     : ("bigg_iYS1720.txt",       4006,  2445),
    "RECON1"      : ("bigg_RECON1.txt",        5290,  2766),
    "iCHOv1"      : ("bigg_iCHOv1.txt",       9289,  4456),
    "Recon3D"     : ("bigg_Recon3D.txt",      15834,  5835),
}

import argparse
_p = argparse.ArgumentParser(description="BiGG ERC synergy explorer")
_p.add_argument("network", nargs="?", default=None,
                help="Network key (e.g. iAF1260)")
_p.add_argument("--synergy", choices=["basic", "maximal", "fundamental"],
                default=SYNERGY_TYPE)
_p.add_argument("--no-verify", action="store_true",
                help="Skip ERC structural assertions (faster on large nets)")
_p.add_argument("--list", action="store_true",
                help="List all available networks and exit")
args = _p.parse_args()

if args.list:
    print(f"{'key':<16}  {'rxns':>6}  {'species':>7}  file")
    for k, (f, r, s) in sorted(CANDIDATES.items(), key=lambda kv: kv[1][1]):
        print(f"  {k:<16}  {r:>6}  {s:>7}  {f}")
    sys.exit(0)

if args.network:
    SELECTED = args.network
SYNERGY_TYPE = args.synergy
VERIFY_ERCS  = not args.no_verify

if SELECTED not in CANDIDATES:
    print(f"Unknown network '{SELECTED}'. Run with --list to see options.")
    sys.exit(1)

fname = CANDIDATES[SELECTED][0]
NET_PATH = os.path.join(_repo, "data", "biomodels", "BiGG", fname)

# ── imports ───────────────────────────────────────────────────────────────────
from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations.io_pyCOT import build_rndata
from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies
from pyCOT.analysis.organizations.complementarity import compute_complementarities
from pyCOT.analysis.organizations.generators import compute_generators
from cot_gen.metanetwork import build_metanetwork
from pyCOT.analysis.organizations.metrics import StageContext, Counters, print_summary

# ── banner ────────────────────────────────────────────────────────────────────
SEP = "=" * 70
print(SEP)
print(f"BiGG ERC synergy  —  {SELECTED}  ({fname})")
print(f"synergy level: {SYNERGY_TYPE}   verify: {VERIFY_ERCS}")
print(SEP)

# ── Stage 0: load ─────────────────────────────────────────────────────────────
print("\n[Stage 0] Loading network...")
t0 = time.perf_counter()
rn_pycot = read_txt(NET_PATH)
load_ms = (time.perf_counter() - t0) * 1000

with StageContext("stage0", SELECTED, append_log=True) as ctx0:
    rn = build_rndata(rn_pycot, network_id=SELECTED)
ctx0.row.n_species   = rn.n_species
ctx0.row.n_reactions = rn.n_reactions

n_E0         = bin(rn.E0_mask).count('1')
n_nontrivial = sum(1 for s in rn.supp_q if s != 0)
print(f"  species:          {rn.n_species}")
print(f"  reactions:        {rn.n_reactions}")
print(f"  E0 species:       {n_E0}")
print(f"  non-trivial rxns: {n_nontrivial}")
print(f"  file load:        {load_ms:.1f} ms")
print(f"  Stage-0 build:    {ctx0.row.wall_s*1000:.1f} ms")

# ── Stage 1: ERCs ─────────────────────────────────────────────────────────────
verify_note = "" if VERIFY_ERCS else "  (assertions OFF)"
print(f"\n[Stage 1] Computing ERCs...{verify_note}")
ctr1 = Counters()
with StageContext("stage1", SELECTED,
                  n_species=rn.n_species, n_reactions=rn.n_reactions,
                  counters=ctr1, append_log=True) as ctx1:
    ercs = compute_ercs(rn, counters=ctr1, verify=VERIFY_ERCS)
ctx1.row.n_ercs = len(ercs)

n_persistent = sum(1 for e in ercs if e.is_persistent())
print(f"  ERCs found:       {len(ercs)}")
print(f"  persistent ERCs:  {n_persistent}")
print(f"  min-bases total:  {ctr1.get('erc.min_bases_total', 0)}")
print(f"  wall:             {ctx1.row.wall_s*1000:.1f} ms")
print(f"  peak mem:         {ctx1.row.peak_mb:.2f} MB")

if not ercs:
    print("\nNo ERCs — nothing to synergize.")
    sys.exit(0)

# ── Stage H: hierarchy ────────────────────────────────────────────────────────
print(f"\n[Stage H] Building ERC containment hierarchy...")
ctrH = Counters()
with StageContext("stage_hier", SELECTED,
                  n_species=rn.n_species, n_reactions=rn.n_reactions,
                  n_ercs=len(ercs), counters=ctrH, append_log=True) as ctxH:
    hier = build_hierarchy(ercs)

n_pairs_total = len(ercs) * (len(ercs) - 1) // 2
n_comparable  = sum(len(hier.ancestors[i]) for i in range(hier.n))
print(f"  ERC pairs total:  {n_pairs_total}")
print(f"  comparable pairs: {n_comparable}")
print(f"  Hasse edges:      {sum(len(hier.parents[i]) for i in range(hier.n))}")
print(f"  wall:             {ctxH.row.wall_s*1000:.1f} ms")

# ── Stage S: synergies ────────────────────────────────────────────────────────
print(f"\n[Stage S] Computing {SYNERGY_TYPE} synergies...")
ctrS = Counters()
with StageContext("stage_syn", SELECTED,
                  n_species=rn.n_species, n_reactions=rn.n_reactions,
                  n_ercs=len(ercs), counters=ctrS, append_log=True) as ctxS:
    result = compute_synergies(ercs, hier, level=SYNERGY_TYPE, counters=ctrS)

n_basic  = len(result.basic)
n_max    = len(result.maximal)
n_fund   = len(result.fundamental)
syn_s    = ctxS.row.wall_s
n_incomp = n_pairs_total - n_comparable   # pairs that can interact

print(f"  incomparable pairs:   {n_incomp}")
print(f"  basic synergies:      {n_basic}")
if SYNERGY_TYPE in ("maximal", "fundamental"):
    print(f"  maximal synergies:    {n_max}")
if SYNERGY_TYPE == "fundamental":
    print(f"  fundamental syn.:     {n_fund}")
print(f"  wall:                 {syn_s*1000:.1f} ms")
if syn_s > 0:
    print(f"  pairs/s:              {n_pairs_total/syn_s:,.0f}")

# ── Stage C: complementarities ────────────────────────────────────────────────
print(f"\n[Stage C] Computing complementarities...")
ctrC = Counters()
with StageContext("stage_comp", SELECTED,
                  n_species=rn.n_species, n_reactions=rn.n_reactions,
                  n_ercs=len(ercs), counters=ctrC, append_log=True) as ctxC:
    comp_result = compute_complementarities(ercs, hier, counters=ctrC)

print(f"  Type 1 (req reduction):  {len(comp_result.type1)}")
print(f"  Type 2 (req overlap):    {len(comp_result.type2)}")
print(f"  Type 3 (prod extension): {len(comp_result.type3)}")
print(f"  total:                   {len(comp_result)}")
print(f"  wall:                    {ctxC.row.wall_s*1000:.1f} ms")

# ── Stage G: generators ───────────────────────────────────────────────────────
print(f"\n[Stage G] Computing generators (primitive ERCs)...")
# Generators use fundamental synergies; re-compute if not already available
if SYNERGY_TYPE != "fundamental":
    print("  (re-computing fundamental synergies for generator analysis)")
    fund_result = compute_synergies(ercs, hier, level="fundamental")
else:
    fund_result = result

ctrG = Counters()
with StageContext("stage_gen", SELECTED,
                  n_species=rn.n_species, n_reactions=rn.n_reactions,
                  n_ercs=len(ercs), counters=ctrG, append_log=True) as ctxG:
    gen_result = compute_generators(ercs, hier, fund_result, counters=ctrG)

print(f"  primitive ERCs:  {len(gen_result.primitive_indices)}"
      f"  of {len(ercs)} total")
print(f"  basis reach:     {len(gen_result.basis_reach)}")
print(f"  coverage:        {gen_result.coverage:.1%}")
print(f"  complete basis:  {gen_result.is_complete()}")
print(f"  wall:            {ctxG.row.wall_s*1000:.1f} ms")

# ── MetaNetwork summary ───────────────────────────────────────────────────────
mn = build_metanetwork(ercs, hier, fund_result, comp_result, gen_result)
print(f"\n{SEP}")
print("METANETWORK SUMMARY")
print(SEP)
mn.print_summary(prefix="  ")

# ── Sample output ─────────────────────────────────────────────────────────────
active = result.fundamental or result.maximal or result.basic
if active:
    print(f"\n{SEP}")
    level_label = ("fundamental" if result.fundamental else
                   "maximal"     if result.maximal     else "basic")
    print(f"SAMPLE {level_label.upper()} SYNERGIES (first 20)")
    print(SEP)
    for s in active[:20]:
        ei = ercs[s.i]; ej = ercs[s.j]; ek = ercs[s.k]
        print(f"  E{s.i}(|{ei.size()}|) + E{s.j}(|{ej.size()}|)"
              f"  ->  E{s.k}(|{ek.size()}|)"
              f"  [req={bin(ek.req_mask).count('1')} sp]")

# ── Metrics summary ───────────────────────────────────────────────────────────
print(f"\n{SEP}")
print("STAGE METRICS")
print(SEP)
print_summary()
