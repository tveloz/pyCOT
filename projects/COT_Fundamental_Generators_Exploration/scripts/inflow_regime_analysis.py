"""
inflow_regime_analysis.py — Centler-style inflow-regime comparison.

Reconstructs, for any network, the methodology of Centler, Speroni di
Fenizio, Matsumaru & Dittrich (2006): inflow means environment, and
different environments (inflow regimes) reveal different organizational
hierarchies from the SAME underlying reaction network.

Three phases, run and reported together:

  Phase A — ISOLATED structure (no inflow at all)
    Strip every native inflow reaction ("=> X") from the network and
    compute ERCs + the fundamental hierarchy (synergy/complementarity) on
    what remains: the network "closed off from the universe". This is
    the cheap, purely structural req/prod statistics the network's own
    reaction wiring implies, independent of any environmental assumption.

  Phase B — REGIME structure (one or more chosen food sets)
    For each named scenario (a list of species tokens to feed in via a
    fresh "=> X" reaction, replacing whatever native inflows existed),
    recompute the same ERC/hierarchy statistics and compare against
    Phase A: which species/ERCs only become reachable once this regime's
    food is supplied, how req/prod distributions shift, how hierarchy
    density (fundamental synergies/complementarities) changes.

  Phase C — ORGANIZATIONS per regime (optional, costlier)
    Run pyCOT.analysis.organizations.compute_organizations() under each
    regime and report the verified-organization count/sizes -- directly
    the Centler-paper-style comparison ("this environment yields N
    organizations of these sizes"). Capped by max_espm_order for
    genome-scale safety (see this session's iAF692 timing characterization:
    ERC/hierarchy/EPM are fast at any scale, ESPM order>=2 is not yet
    profiled/safe at genome scale -- default here is conservative).

Reuses the same text-level scenario-building primitive already validated
in projects/RAF_Comparison/raf/inflow_scenarios.py (strip native inflows,
append a fresh food-species reaction per scenario) -- ported here rather
than cross-imported, since RAF_Comparison's copy is framed around
CRS/RAF-specific concerns (Def:translation's Omega) that don't apply to
this module's plain COT-hierarchy framing, and the function is small and
stable enough that a second copy is not a maintenance risk.

Usage (from repo root)
-----------------------
    python projects/COT_Fundamental_Generators_Exploration/scripts/inflow_regime_analysis.py \\
        <path_to_network.txt> <network_name> [--organizations] [--max-espm-order N]

    Or import and call analyze_inflow_regimes(path, name, scenarios) from
    another script, where scenarios is a list of (label, food_tokens) pairs.
    An implicit ("isolated", []) scenario is always run first as the
    baseline regardless of what's passed.
"""
from __future__ import annotations

import sys
import os
import argparse

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations import (
    build_rndata, compute_ercs, build_hierarchy,
    compute_synergies_basis_first, compute_complementarities,
    compute_organizations,
)

OUT_ROOT = os.path.join(_repo_root, 'projects', 'COT_Fundamental_Generators_Exploration',
                         'outputs', 'inflow_regimes')


# ═══════════════════════════════════════════════════════════════════════
# Scenario-network construction (text-level, see module docstring)
# ═══════════════════════════════════════════════════════════════════════

def strip_native_inflows(txt: str) -> list[str]:
    """Every reaction line of `txt` EXCEPT its native inflow reactions
    (empty left-hand side). Outflow reactions are left untouched."""
    kept = []
    for line in txt.splitlines():
        line = line.rstrip()
        if not line.strip():
            continue
        body = line.split(";", 1)[0]
        if ":" not in body or "=>" not in body:
            continue
        lhs = body.split(":", 1)[1].split("=>", 1)[0].strip()
        if lhs == "":
            continue  # drop native inflow
        kept.append(line)
    return kept


def build_scenario_txt(base_txt: str, food_tokens: list[str]) -> str:
    """
    The SAME real (non-inflow) reactions, with native inflows replaced by
    one fresh "=> token" reaction per entry in food_tokens. food_tokens
    must be exact species tokens as they already appear in base_txt (e.g.
    "glc__D_e" for bare BiGG convention), so they resolve to the same
    species pyCOT already knows about rather than minting a new one.
    """
    lines = strip_native_inflows(base_txt)
    for i, tok in enumerate(food_tokens):
        lines.append(f"INFLOW_{i}:  => 1 {tok};")
    return "\n".join(lines) + "\n"


# ═══════════════════════════════════════════════════════════════════════
# Per-scenario structural profile (Phases A/B)
# ═══════════════════════════════════════════════════════════════════════

def _profile_scenario(txt: str, label: str) -> dict:
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False,
                                      encoding='utf-8') as f:
        f.write(txt)
        tmp_path = f.name
    try:
        rn = read_txt(tmp_path, exact_names=True)
        rn_data = build_rndata(rn, network_id=label)
        ercs = compute_ercs(rn_data, verify=False)
        hier = build_hierarchy(ercs)
        syn = compute_synergies_basis_first(ercs, hier)
        comp = compute_complementarities(ercs, hier, syn)

        req_sizes = [bin(e.req_mask).count('1') for e in ercs]
        prod_sizes = [bin(e.prod_mask).count('1') for e in ercs]
        n_persistent = sum(1 for e in ercs if e.is_persistent())

        return {
            'label': label,
            'n_species_total': rn_data.n_species,
            'n_reactions': rn_data.n_reactions,
            'e0_size': bin(rn_data.E0_mask).count('1'),
            'n_ercs': len(ercs),
            'n_persistent_ercs': n_persistent,
            'req_mean': sum(req_sizes) / len(req_sizes) if req_sizes else 0.0,
            'req_max': max(req_sizes, default=0),
            'prod_mean': sum(prod_sizes) / len(prod_sizes) if prod_sizes else 0.0,
            'prod_max': max(prod_sizes, default=0),
            'n_synergies': len(syn.fundamental),
            'n_complementarities': len(comp.fundamental),
            '_rn': rn, '_rn_data': rn_data, '_txt_path': tmp_path,
        }
    finally:
        pass  # tmp_path cleaned up by caller after optional Phase C reuse


# ═══════════════════════════════════════════════════════════════════════
# Orchestration
# ═══════════════════════════════════════════════════════════════════════

def analyze_inflow_regimes(
    path: str,
    name: str,
    scenarios: list[tuple[str, list[str]]],
    *,
    run_organizations: bool = False,
    max_espm_order: int = 2,
    verbose: bool = True,
) -> str:
    """
    Run Phase A (isolated) + Phase B (each scenario) + optional Phase C
    (organizations), write a comparative report.md + summary.csv to
    outputs/inflow_regimes/<name>/, return the output directory.
    """
    with open(path, 'r', encoding='utf-8') as f:
        base_txt = f.read()

    out_dir = os.path.join(OUT_ROOT, name)
    os.makedirs(out_dir, exist_ok=True)

    all_scenarios = [('isolated (no inflow)', [])] + list(scenarios)
    results = []
    tmp_paths = []

    for label, food_tokens in all_scenarios:
        if verbose:
            print(f"[inflow_regime_analysis] {label}: "
                  f"{len(food_tokens)} food species...")
        txt = build_scenario_txt(base_txt, food_tokens)
        prof = _profile_scenario(txt, label)
        prof['food_tokens'] = food_tokens
        tmp_paths.append(prof['_txt_path'])

        if run_organizations:
            rn = prof.pop('_rn')
            org_result = compute_organizations(
                rn, network_id=f"{name}__{label}",
                max_espm_order=max_espm_order, verbose=False,
            )
            sizes = sorted(len(o.species_names) for o in org_result.organizations)
            prof['n_organizations'] = len(org_result.organizations)
            prof['organization_sizes'] = sizes
        else:
            prof.pop('_rn', None)

        prof.pop('_rn_data', None)
        results.append(prof)

    # ── Report ──────────────────────────────────────────────────────────
    lines = [f"# Inflow-regime analysis: {name}\n",
             f"Source: `{path}`\n",
             f"Scenarios: {len(all_scenarios)} "
             f"({len(scenarios)} named regime(s) + the isolated baseline)\n"]

    header = ("| scenario | food species | ERCs | P-ERCs | req mean/max | "
              "prod mean/max | synergies | complementarities |")
    if run_organizations:
        header += " organizations (sizes) |"
    lines.append(header)
    sep = "|---" * (9 if run_organizations else 8) + "|"
    lines.append(sep)
    for r in results:
        row = (f"| {r['label']} | {', '.join(r['food_tokens']) or '(none)'} "
               f"| {r['n_ercs']} | {r['n_persistent_ercs']} "
               f"| {r['req_mean']:.2f} / {r['req_max']} "
               f"| {r['prod_mean']:.2f} / {r['prod_max']} "
               f"| {r['n_synergies']} | {r['n_complementarities']} |")
        if run_organizations:
            row += f" {r['n_organizations']} {r['organization_sizes']} |"
        lines.append(row)

    lines.append("\n(No EPM/ESPM/organization computation performed unless "
                  "--organizations was passed; Phases A/B are req/prod "
                  "structural statistics only, matching the network's own "
                  "wiring under each regime.)\n")

    report_path = os.path.join(out_dir, 'report.md')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(lines) + "\n")

    import csv
    csv_path = os.path.join(out_dir, 'summary.csv')
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        fixed_cols = ['label', 'food_tokens', 'n_species_total', 'n_reactions',
                      'e0_size', 'n_ercs', 'n_persistent_ercs', 'req_mean',
                      'req_max', 'prod_mean', 'prod_max', 'n_synergies',
                      'n_complementarities']
        if run_organizations:
            fixed_cols += ['n_organizations', 'organization_sizes']
        w.writerow(fixed_cols)
        for r in results:
            row = [r.get(c) for c in fixed_cols]
            w.writerow(row)

    for p in tmp_paths:
        try:
            os.remove(p)
        except OSError:
            pass

    if verbose:
        print(f"[inflow_regime_analysis] -> {out_dir}")
    return out_dir


# ═══════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('path')
    ap.add_argument('name')
    ap.add_argument('--organizations', action='store_true',
                     help='Also run Phase C (verified organizations per scenario)')
    ap.add_argument('--max-espm-order', type=int, default=2)
    args = ap.parse_args()

    print("No scenarios given on the CLI beyond the isolated baseline -- "
          "use analyze_inflow_regimes(path, name, scenarios) from a script "
          "to pass named food-set regimes. Running isolated-only.")
    analyze_inflow_regimes(
        args.path, args.name, [],
        run_organizations=args.organizations,
        max_espm_order=args.max_espm_order,
    )
