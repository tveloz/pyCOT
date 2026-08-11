"""
network_structure_profile.py

Cheap, staged structural characterization of a reaction network -- BEFORE
any organization / EPM / ESPM computation. Two stages, run and reported
separately, exactly as requested:

  STEP 1 (network_structure_profile.py's "erc" stage):
    - Compute ERCs (via cot_gen, the same validated engine used throughout
      this project -- scales to genome-size networks, unlike src/pyCOT's
      classical-organization machinery used for reproduce_centler2006.py).
    - Report req/prod size distributions per ERC (cheap, informative:
      how "needy" vs. "productive" the network's building blocks are).
    - Identify P-ERCs (persistent, req=0 -- already self-sufficient).
    - Identify INFLOW species (appear as the product of a `-> X` reaction
      with an EMPTY left side in the raw .txt file) and OUTFLOW species
      (appear as the reactant of an `X ->` reaction with an EMPTY right
      side). This is the same signal that exposed the iAF692 conversion
      bug (zero true inflow reactions) -- this script flags that class of
      problem explicitly, on any network, as a first-class check.
    - Also reports "structural" source/sink species (produced by SOME
      reaction but never consumed by any, or vice versa) as a secondary,
      weaker signal -- distinct from the EXPLICIT inflow/outflow reactions.

  STEP 2 ("hierarchy" stage):
    - Build the ERC containment hierarchy.
    - Compute fundamental synergies and complementarities (the generative
      structure: which ERCs can jointly unlock or feed each other).
    - Report degree distributions and the density/entanglement metrics
      already validated in epm_reach_propagate.py's compute_structural_metrics.

Deliberately stops here -- no EPM/ESPM/organization computation. The output
of this script is meant to inform WHICH inflow/outflow scenarios are worth
testing next (that's a separate, later step).

Outputs, per network, in projects/COT_Fundamental_Generators_Complex/outputs/network_profile/<name>/:
  report.md            -- human-readable summary of both steps
  erc_table.csv         -- one row per ERC: req_size, prod_size, erc_size, is_persistent
  req_prod_histogram.png
  synergy_complementarity_degree_histogram.png

Usage (from repo root):
    python projects/COT_Fundamental_Generators_Complex/scripts/network_structure_profile.py <path_to_network.txt> [network_name]

Or import and call profile_network(path, name) from another script.
"""
from __future__ import annotations

import sys, os, re, csv

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations.io_pyCOT import build_rndata
from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity import compute_complementarities

OUT_ROOT = os.path.join(_repo_root, 'projects', 'COT_Fundamental_Generators_Complex', 'outputs', 'network_profile')


# ═══════════════════════════════════════════════════════════════════════
# Raw-file inflow/outflow detection (independent of any internal engine
# representation -- works directly off the reaction text, so it's a
# transparent, engine-agnostic cross-check)
# ═══════════════════════════════════════════════════════════════════════
def _parse_raw_inflow_outflow(path: str) -> tuple[set[str], set[str]]:
    """
    Returns (inflow_species, outflow_species) by scanning the raw .txt file
    for reactions with an EMPTY left side (`-> X` = inflow) or an EMPTY
    right side (`X ->` = outflow). Handles both '=>' and '->' arrows and
    ignores comments after ';'.
    """
    inflow, outflow = set(), set()
    # Coefficient prefix is optional and its separating space is optional
    # too, e.g. "2 UFPT" (space-separated) and "2UFPT" (compact notation,
    # common in SBML-derived BioModels conversions) must both strip to
    # species "UFPT". A \s+ (space required) version was missing the
    # compact case entirely, silently dropping the reactant term and thus
    # misclassifying reactions like "2UFPT (UFPT) => UFPT (UFPT)" (a decay
    # reaction) as if they had an empty left side (a true inflow) --
    # confirmed on BIOMD0000000446, which reported 27 "explicit inflow
    # species" via this scan even though only 13 reactions have a genuinely
    # empty left side; the real computation (build_rndata/compute_ercs,
    # which uses pyCOT's actual parser, not this lightweight text scan) was
    # never affected by this -- only this diagnostic's own species_in().
    term_re = re.compile(r'(?:\d+(?:\.\d+)?\s*)?([A-Za-z_][A-Za-z0-9_\'\[\]]*)')

    def species_in(side: str) -> list[str]:
        side = side.strip()
        if not side:
            return []
        return [m.group(1) for term in side.split('+') if (m := term_re.match(term.strip()))]

    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.split(';', 1)[0].strip()
            if not line:
                continue
            if ':' in line:
                line = line.split(':', 1)[1]
            arrow = '=>' if '=>' in line else ('->' if '->' in line else None)
            if arrow is None:
                continue
            lhs, rhs = line.split(arrow, 1)
            lhs_species = species_in(lhs)
            rhs_species = species_in(rhs)
            if not lhs_species and rhs_species:
                inflow.update(rhs_species)
            if lhs_species and not rhs_species:
                outflow.update(lhs_species)
    return inflow, outflow


# ═══════════════════════════════════════════════════════════════════════
# STEP 1: ERC + req/prod structural profile
# ═══════════════════════════════════════════════════════════════════════
def step1_erc_profile(path: str, name: str, out_dir: str) -> dict:
    rn = read_txt(path, exact_names=True)
    rn_data = build_rndata(rn, network_id=name)
    ercs = compute_ercs(rn_data, verify=False)

    inflow_sp, outflow_sp = _parse_raw_inflow_outflow(path)

    req_sizes = [bin(e.req_mask).count('1') for e in ercs]
    prod_sizes = [bin(e.prod_mask).count('1') for e in ercs]
    erc_sizes = [bin(e.species_mask).count('1') for e in ercs]
    n_perc = sum(1 for e in ercs if e.is_persistent())

    # structural (reaction-graph) source/sink species -- a WEAKER, secondary
    # signal than the explicit inflow/outflow reactions above: species that
    # happen to never be consumed (or never produced) by ANY reaction,
    # regardless of whether an explicit boundary reaction marks them so
    produced_ever, consumed_ever = set(), set()
    for s, p in zip(rn_data.supp_raw, rn_data.prod_raw):
        for bit in _bits(s):
            consumed_ever.add(bit)
        for bit in _bits(p):
            produced_ever.add(bit)
    species_names = rn_data.species_names
    structural_sources = {species_names[b] for b in produced_ever - consumed_ever}
    structural_sinks = {species_names[b] for b in consumed_ever - produced_ever}

    os.makedirs(out_dir, exist_ok=True)

    # --- CSV: one row per ERC ---
    csv_path = os.path.join(out_dir, 'erc_table.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['erc_id', 'req_size', 'prod_size', 'erc_size', 'is_persistent'])
        for i, e in enumerate(ercs):
            w.writerow([i, bin(e.req_mask).count('1'), bin(e.prod_mask).count('1'),
                        bin(e.species_mask).count('1'), e.is_persistent()])

    # --- plot: req/prod/erc size histograms ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, data, title in zip(axes, [req_sizes, prod_sizes, erc_sizes],
                                ['req size per ERC\n(how many species this module still needs)',
                                 'prod size per ERC\n(how many species this module can make)',
                                 'ERC size per ERC\n(total species touched)']):
        ax.hist(data, bins=_int_bins(data))
        ax.set_title(title, fontsize=9)
        ax.set_xlabel('size')
        ax.set_ylabel('number of ERCs')
    fig.suptitle(f'{name}: req/prod/ERC size distributions ({len(ercs)} ERCs)')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'req_prod_histogram.png'), dpi=150)
    plt.close(fig)

    result = {
        'n_species': rn_data.n_species,
        'n_reactions': rn_data.n_reactions,
        'n_ercs': len(ercs),
        'n_perc': n_perc,
        'req_size_mean': sum(req_sizes) / max(len(req_sizes), 1),
        'req_size_max': max(req_sizes, default=0),
        'prod_size_mean': sum(prod_sizes) / max(len(prod_sizes), 1),
        'prod_size_max': max(prod_sizes, default=0),
        'inflow_species': sorted(inflow_sp),
        'outflow_species': sorted(outflow_sp),
        'n_inflow_reactions_species': len(inflow_sp),
        'n_outflow_reactions_species': len(outflow_sp),
        'structural_sources_only': sorted(structural_sources - inflow_sp),
        'structural_sinks_only': sorted(structural_sinks - outflow_sp),
    }
    return result, rn_data, ercs


def _int_bins(data, max_bins=30):
    """Integer-aligned histogram bin edges -- avoids matplotlib's default
    fractional-width bins looking wrong for discrete count data (e.g. a
    req-size distribution over {0,1} showing a spurious split at 0.5)."""
    lo, hi = min(data, default=0), max(data, default=0)
    n = hi - lo + 2
    if n <= max_bins:
        return [x - 0.5 for x in range(lo, hi + 2)]
    return max_bins


def _bits(mask: int):
    m = mask
    while m:
        lsb = m & (-m)
        yield lsb.bit_length() - 1
        m &= m - 1


# ═══════════════════════════════════════════════════════════════════════
# STEP 2: hierarchy + fundamental synergy/complementarity
# ═══════════════════════════════════════════════════════════════════════
def step2_hierarchy_profile(rn_data, ercs, name: str, out_dir: str) -> dict:
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)

    n = len(ercs)
    n_pairs = n * (n - 1) // 2
    n_comparable = sum(len(hier.ancestors[i]) for i in range(hier.n))
    hasse_edges = sum(len(hier.parents[i]) for i in range(hier.n))

    syn_degree = [0] * n
    for st in syn.fundamental:
        syn_degree[st.i] += 1
        syn_degree[st.j] += 1
    comp_out_degree = [0] * n
    comp_in_degree = [0] * n
    for fc in comp.fundamental:
        comp_out_degree[fc.prod_idx] += 1
        comp_in_degree[fc.cons_idx] += 1

    os.makedirs(out_dir, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, data, title in zip(
        axes, [syn_degree, comp_out_degree, comp_in_degree],
        ['synergy out-degree per ERC\n(how many partners unlock something with it)',
         'complementarity out-degree per ERC\n(how many ERCs it can feed)',
         'complementarity in-degree per ERC\n(how many ERCs can feed it)']):
        ax.hist(data, bins=_int_bins(data))
        ax.set_title(title, fontsize=9)
        ax.set_xlabel('degree')
        ax.set_ylabel('number of ERCs')
    fig.suptitle(f'{name}: fundamental-relation degree distributions')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'synergy_complementarity_degree_histogram.png'), dpi=150)
    plt.close(fig)

    return {
        'n_hasse_edges': hasse_edges,
        'n_comparable_pairs': n_comparable,
        'n_incomparable_pairs': n_pairs - n_comparable,
        'n_fundamental_synergies': len(syn.fundamental),
        'n_fundamental_complementarities': len(comp.fundamental),
        'syn_density': len(syn.fundamental) / max(n, 1),
        'comp_density': len(comp.fundamental) / max(n, 1),
        'max_syn_outdegree': max(syn_degree, default=0),
    }


# ═══════════════════════════════════════════════════════════════════════
# Orchestration + report
# ═══════════════════════════════════════════════════════════════════════
def profile_network(path: str, name: str | None = None) -> str:
    if name is None:
        name = os.path.splitext(os.path.basename(path))[0]
    out_dir = os.path.join(OUT_ROOT, name)
    os.makedirs(out_dir, exist_ok=True)

    print(f'=== {name} ===')
    print('[Step 1] Computing ERCs and req/prod structure...')
    s1, rn_data, ercs = step1_erc_profile(path, name, out_dir)
    print(f'  {s1["n_species"]} species, {s1["n_reactions"]} reactions, {s1["n_ercs"]} ERCs '
          f'({s1["n_perc"]} persistent / P-ERCs)')
    print(f'  req size: mean={s1["req_size_mean"]:.2f} max={s1["req_size_max"]}   '
          f'prod size: mean={s1["prod_size_mean"]:.2f} max={s1["prod_size_max"]}')
    print(f'  explicit inflow species: {s1["n_inflow_reactions_species"]}   '
          f'explicit outflow species: {s1["n_outflow_reactions_species"]}')
    if s1['n_inflow_reactions_species'] == 0:
        print('  !!! WARNING: ZERO explicit inflow reactions found -- this network cannot be '
              'self-maintaining in ANY sense as given (see iAF692 case study). Check the source file.')

    print('[Step 2] Computing ERC hierarchy and fundamental relations...')
    s2 = step2_hierarchy_profile(rn_data, ercs, name, out_dir)
    print(f'  Hasse edges: {s2["n_hasse_edges"]}   fundamental synergies: {s2["n_fundamental_synergies"]}   '
          f'fundamental complementarities: {s2["n_fundamental_complementarities"]}')
    print(f'  syn_density={s2["syn_density"]:.2f}  comp_density={s2["comp_density"]:.2f}  '
          f'max_syn_outdegree={s2["max_syn_outdegree"]}')

    # --- report.md ---
    lines = [
        f'# Network structure profile: {name}\n',
        f'Source file: `{path}`\n',
        '## Step 1: ERCs and req/prod structure\n',
        f'- Species: {s1["n_species"]}   Reactions: {s1["n_reactions"]}   ERCs: {s1["n_ercs"]}',
        f'- P-ERCs (persistent, req=0): {s1["n_perc"]}',
        f'- req size per ERC: mean={s1["req_size_mean"]:.2f}, max={s1["req_size_max"]}',
        f'- prod size per ERC: mean={s1["prod_size_mean"]:.2f}, max={s1["prod_size_max"]}',
        f'- Explicit inflow species ({s1["n_inflow_reactions_species"]}): {", ".join(s1["inflow_species"]) or "(none)"}',
        f'- Explicit outflow species ({s1["n_outflow_reactions_species"]}): {", ".join(s1["outflow_species"][:40])}'
        + (' ... (truncated)' if len(s1["outflow_species"]) > 40 else ''),
        f'- Structural-only sources (produced, never consumed, no explicit inflow reaction): '
        f'{", ".join(s1["structural_sources_only"][:20]) or "(none)"}'
        + (' ... (truncated)' if len(s1["structural_sources_only"]) > 20 else ''),
        f'- Structural-only sinks (consumed, never produced, no explicit outflow reaction): '
        f'{", ".join(s1["structural_sinks_only"][:20]) or "(none)"}'
        + (' ... (truncated)' if len(s1["structural_sinks_only"]) > 20 else ''),
        '',
        '## Step 2: ERC hierarchy and fundamental relations\n',
        f'- Hasse (containment) edges: {s2["n_hasse_edges"]}',
        f'- Comparable / incomparable ERC pairs: {s2["n_comparable_pairs"]} / {s2["n_incomparable_pairs"]}',
        f'- Fundamental synergies: {s2["n_fundamental_synergies"]}  (density {s2["syn_density"]:.2f} per ERC)',
        f'- Fundamental complementarities: {s2["n_fundamental_complementarities"]}  (density {s2["comp_density"]:.2f} per ERC)',
        f'- Max synergy out-degree: {s2["max_syn_outdegree"]}',
        '',
        'Files in this directory: `erc_table.csv`, `req_prod_histogram.png`, '
        '`synergy_complementarity_degree_histogram.png`.',
        '',
        '(No organization / EPM / ESPM computation performed -- structural profile only.)',
    ]
    if s1['n_inflow_reactions_species'] == 0:
        lines.insert(1, '**WARNING: zero explicit inflow reactions found in this file.**\n')

    report_path = os.path.join(out_dir, 'report.md')
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print(f'  -> {out_dir}\n')
    return out_dir


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    path = sys.argv[1]
    name = sys.argv[2] if len(sys.argv) > 2 else None
    profile_network(path, name)
