"""
╔════════════════════════════════════════════════════════════════════════════════╗
║  compare_oracles.py — Brute-force oracle vs efficient algorithm comparison   ║
╚════════════════════════════════════════════════════════════════════════════════╝

WHAT IT DOES
------------
For each network in the selected size range, runs BOTH the brute-force oracle
AND the efficient algorithm for four stages and checks that they agree:

  Stage 1 — ERCs             : oracle enumerates 2^|S| subsets; efficient uses Horn closure.
  Stage S — Fundamental syn. : oracle is a corrected 3-step (basic→maximal→fundamental);
                               efficient uses the bitset hierarchy filter.
  Stage C — Fundamental comp.: oracle and efficient both apply Def 26 (minprod/mincons);
                               comparison verifies (prod_idx, cons_idx, species_bit) sets.
  Stage E — EPMs             : oracle enumerates ERC-subsets; efficient uses minimal P-ERCs
                               + pairwise complementary pairs.  (Only run when |ERCs| ≤ EPM_ERC_LIMIT.)

Skipped stages are marked with '--' in the table.

HOW TO RUN
----------
  1. Edit the CONFIGURATION section below.
  2. Press ▶ (play) in VS Code, or run:
       python projects/COT_fundamental_Generators/scripts/compare_oracles.py
"""

# ── Path setup ────────────────────────────────────────────────────────────────
from __future__ import annotations
import csv, os, sys, time
_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyCOT.io.functions              import read_txt
from pyCOT.analysis.organizations.io_pyCOT        import build_rndata
from pyCOT.analysis.organizations.erc              import compute_ercs
from pyCOT.analysis.organizations.hierarchy        import build_hierarchy
from pyCOT.analysis.organizations.synergy          import compute_synergies
from pyCOT.analysis.organizations.complementarity  import compute_complementarities
from pyCOT.analysis.organizations.epm              import compute_epms
from oracles.erc_oracle              import compute_ercs_oracle
from oracles.synergy_oracle          import fundamental_synergy_set as oracle_fund_syn
from oracles.complementarity_oracle  import comp_fund_set           as oracle_fund_comp
from oracles.epm_oracle              import epm_oracle_set

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — edit this block, then press Play                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ── Network filter ────────────────────────────────────────────────────────────
MIN_REACTIONS = 100
MAX_REACTIONS = 200

# Which collection: "all" | "bigg" | "biomd"
NETWORK_FILTER = "all"

# ── EPM oracle size limit ─────────────────────────────────────────────────────
# EPM oracle enumerates 2^|ERCs| subsets — only feasible for small networks.
# Set to 0 to skip EPM comparison entirely.
EPM_ERC_LIMIT = 20

# ── Catalogue ─────────────────────────────────────────────────────────────────
CATALOGUE_CSV = "projects/COT_fundamental_Generators/network_catalogue.csv"

# ── Output ────────────────────────────────────────────────────────────────────
OUTPUT_CSV = "compare_oracles_results.csv"

# ── Verification ─────────────────────────────────────────────────────────────
VERIFY = True

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Script body                                                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

_DATA_ROOT = os.path.join(_repo, "data", "biomodels")
_OUT_PATH  = os.path.join(_here, OUTPUT_CSV)


def _is_bigg(path: str) -> bool:
    p = path.replace("\\", "/")
    return "BiGG" in p or os.path.basename(p).startswith("bigg_")


def _is_biomd(path: str) -> bool:
    p = path.replace("\\", "/")
    fname = os.path.basename(p)
    return ("BioMD" in p or "BIOMD" in p) and not fname.startswith("bigg_")


def _discover_all() -> list[tuple[str, str]]:
    seen: dict[str, tuple[str, int]] = {}
    for root, _, files in os.walk(_DATA_ROOT):
        for fname in sorted(files):
            if not fname.endswith(".txt"):
                continue
            full = os.path.join(root, fname)
            name = os.path.splitext(fname)[0]
            if name.startswith("bigg_"):
                name = name[5:]
            if name not in seen or len(full) < seen[name][1]:
                seen[name] = (full, len(full))
    return [(name, path) for name, (path, _) in sorted(seen.items())]


def _load_catalogue(path: str) -> list[tuple[str, str, int]]:
    seen_names: set[str] = set()
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            name = row["name"]
            if name in seen_names:
                continue
            seen_names.add(name)
            try:
                n = int(row["reactions"])
            except (ValueError, KeyError):
                n = -1
            rows.append((name, row["rel_path"], n))
    return rows


# ── Build candidate list ──────────────────────────────────────────────────────
if CATALOGUE_CSV and os.path.isfile(CATALOGUE_CSV):
    cat = _load_catalogue(CATALOGUE_CSV)
    candidates = [
        (name, os.path.join(_repo, rel_path), n)
        for name, rel_path, n in cat
        if MIN_REACTIONS <= n <= MAX_REACTIONS
    ]
else:
    if CATALOGUE_CSV:
        print(f"[warn] Catalogue not found at '{CATALOGUE_CSV}', scanning all files...")
    all_nets = _discover_all()
    candidates = []
    for name, path in all_nets:
        try:
            rn_tmp = read_txt(path)
            n      = len(rn_tmp.reactions())
        except Exception:
            continue
        if MIN_REACTIONS <= n <= MAX_REACTIONS:
            candidates.append((name, path, n))

if NETWORK_FILTER == "bigg":
    candidates = [(n, p, r) for n, p, r in candidates if _is_bigg(p)]
elif NETWORK_FILTER == "biomd":
    candidates = [(n, p, r) for n, p, r in candidates if _is_biomd(p)]

candidates.sort(key=lambda x: x[2])

print(f"Found {len(candidates)} unique networks with {MIN_REACTIONS}–{MAX_REACTIONS} reactions "
      f"(filter: {NETWORK_FILTER}).")
if not candidates:
    print("Nothing to compare. Increase MAX_REACTIONS or check CATALOGUE_CSV.")
    sys.exit(0)

# ── Column layout ─────────────────────────────────────────────────────────────
COLS = [
    "network", "reactions", "ercs",
    "erc_match",  "t_erc_eff_ms",  "t_erc_orc_ms",
    "fund_syn",   "syn_match",     "t_syn_eff_ms",  "t_syn_orc_ms",
    "fund_comp",  "comp_match",    "t_comp_eff_ms", "t_comp_orc_ms",
    "epms",       "epm_match",     "t_epm_eff_ms",  "t_epm_orc_ms",
    "all_ok", "error",
]


def _ms(v) -> str:
    return f"{v:7.1f}" if v is not None else "      -"


def _ok(v) -> str:
    if v is None:  return "  - "
    return " OK " if v else "FAIL"


# ── Table header ──────────────────────────────────────────────────────────────
SEP = "=" * 130
print(SEP)
print(f"  {'Network':<32} {'Rxn':>4}  "
      f"{'ERCs':>5} {'eff ms':>7} {'orc ms':>7} {'ERC?':>5}  "
      f"{'fSYN':>5} {'eff ms':>7} {'orc ms':>7} {'SYN?':>5}  "
      f"{'fCMP':>5} {'eff ms':>7} {'orc ms':>7} {'CMP?':>5}  "
      f"{'EPMs':>5} {'eff ms':>7} {'orc ms':>7} {'EPM?':>5}  ALL")
print(f"  {'(eff=efficient  orc=brute-force oracle  --=skipped)'}")
print(SEP)

rows_out: list[dict] = []

for net_name, net_path, n_rxn in candidates:
    row: dict = {col: None for col in COLS}
    row["network"]   = net_name
    row["reactions"] = n_rxn
    row["error"]     = ""

    try:
        rn_pycot = read_txt(net_path)
        rn       = build_rndata(rn_pycot, network_id=net_name)

        # ── ERCs ──────────────────────────────────────────────────────────────
        t0 = time.perf_counter()
        ercs = compute_ercs(rn, verify=VERIFY)
        row["t_erc_eff_ms"] = (time.perf_counter() - t0) * 1000

        t0 = time.perf_counter()
        ercs_orc = compute_ercs_oracle(list(rn.supp_q), list(rn.prod_q))
        row["t_erc_orc_ms"] = (time.perf_counter() - t0) * 1000

        masks_eff = sorted(e.species_mask for e in ercs)
        masks_orc = sorted(e["species_mask"] for e in ercs_orc)
        row["ercs"]      = len(ercs)
        row["erc_match"] = (masks_eff == masks_orc)

        hier = build_hierarchy(ercs)

        # ── Fundamental synergies ─────────────────────────────────────────────
        t0 = time.perf_counter()
        syn = compute_synergies(ercs, hier, level="fundamental")
        row["t_syn_eff_ms"] = (time.perf_counter() - t0) * 1000

        t0 = time.perf_counter()
        orc_fsyn = oracle_fund_syn(ercs)
        row["t_syn_orc_ms"] = (time.perf_counter() - t0) * 1000

        eff_fsyn = {(s.i, s.j, s.k) for s in syn.fundamental}
        row["fund_syn"]  = len(syn.fundamental)
        row["syn_match"] = (eff_fsyn == orc_fsyn)

        # ── Fundamental complementarity ───────────────────────────────────────
        t0 = time.perf_counter()
        comp = compute_complementarities(ercs, hier, syn_result=syn)
        row["t_comp_eff_ms"] = (time.perf_counter() - t0) * 1000

        t0 = time.perf_counter()
        orc_fcomp = oracle_fund_comp(ercs)
        row["t_comp_orc_ms"] = (time.perf_counter() - t0) * 1000

        eff_fcomp = {(fc.prod_idx, fc.cons_idx, fc.species) for fc in comp.fundamental}
        row["fund_comp"]  = len(comp.fundamental)
        row["comp_match"] = (eff_fcomp == orc_fcomp)

        # ── EPMs ──────────────────────────────────────────────────────────────
        n_ercs = len(ercs)
        if EPM_ERC_LIMIT > 0 and n_ercs <= EPM_ERC_LIMIT:
            t0 = time.perf_counter()
            epm = compute_epms(rn, ercs, hier, comp_result=comp)
            row["t_epm_eff_ms"] = (time.perf_counter() - t0) * 1000

            t0 = time.perf_counter()
            orc_epms = epm_oracle_set(rn, ercs)
            row["t_epm_orc_ms"] = (time.perf_counter() - t0) * 1000

            eff_epms = set(epm.all_epm_masks)
            row["epms"]      = len(epm.all_epm_masks)
            row["epm_match"] = (eff_epms == orc_epms)
        # else: leave EPM columns as None (shown as --)

        row["all_ok"] = bool(
            row["erc_match"] and row["syn_match"] and row["comp_match"]
            and (row["epm_match"] is None or row["epm_match"])
        )

    except Exception as exc:
        row["error"]  = str(exc)[:120]
        row["all_ok"] = False

    rows_out.append(row)

    e_n = f"{row['ercs']}"      if row["ercs"]       is not None else "?"
    s_n = f"{row['fund_syn']}"  if row["fund_syn"]   is not None else "?"
    c_n = f"{row['fund_comp']}" if row["fund_comp"]  is not None else "?"
    p_n = f"{row['epms']}"      if row["epms"]       is not None else "-"
    all_str = " OK " if row["all_ok"] else "FAIL"

    line = (
        f"  {net_name:<32} {n_rxn:>4}  "
        f"{e_n:>5} {_ms(row['t_erc_eff_ms'])} {_ms(row['t_erc_orc_ms'])} {_ok(row['erc_match'])}  "
        f"{s_n:>5} {_ms(row['t_syn_eff_ms'])} {_ms(row['t_syn_orc_ms'])} {_ok(row['syn_match'])}  "
        f"{c_n:>5} {_ms(row['t_comp_eff_ms'])} {_ms(row['t_comp_orc_ms'])} {_ok(row['comp_match'])}  "
        f"{p_n:>5} {_ms(row['t_epm_eff_ms'])} {_ms(row['t_epm_orc_ms'])} {_ok(row['epm_match'])}  "
        f"{all_str}"
    )
    if row["error"]:
        line += f"  ERROR: {row['error'][:50]}"
    print(line)

# ── Summary ───────────────────────────────────────────────────────────────────
n_ok   = sum(1 for r in rows_out if r["all_ok"])
n_fail = len(rows_out) - n_ok
print(SEP)
print(f"PASSED: {n_ok}/{len(rows_out)}    FAILED: {n_fail}")
print(f"  eff ms = efficient algorithm wall time")
print(f"  orc ms = brute-force oracle wall time (should match but be slower)")
print(f"  --     = stage skipped (EPM oracle only runs when |ERCs| ≤ {EPM_ERC_LIMIT})")
if n_fail:
    print("\nFailed networks:")
    for r in rows_out:
        if not r["all_ok"]:
            print(f"  {r['network']:<40}  {r['error'][:80]}")

# ── Write CSV ─────────────────────────────────────────────────────────────────
with open(_OUT_PATH, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=COLS)
    writer.writeheader()
    writer.writerows(rows_out)
print(f"\nResults saved to: {_OUT_PATH}")
