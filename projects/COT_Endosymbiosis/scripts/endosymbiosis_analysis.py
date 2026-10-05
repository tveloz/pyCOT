"""
endosymbiosis_analysis.py -- generic before/after Elementary
Semi-Organization (ESO) / Elementary Organization (EO) comparison for an
endosymbiotic merger (THEORY.md Section 2.3).

Terminology note (nomenclature change from the original version of this
project): what this module previously called "EPMs" are Elementary
Semi-Organizations (ESOs) -- closed sets built directly from ERCs via
fundamental synergy/complementarity, but not yet LP-verified for
self-maintenance. An Elementary Organization (EO) is an ESO that IS
LP-verified self-maintaining (self_maintenance.check_self_maintenance).
See projects/COT_Fundamental_Generators_Exploration/scripts/
eso_eo_analysis.py's module docstring for the full definition and for the
code-level confirmation that the ESO search itself is ERC-based (combining
precomputed ERC species-masks), not a species-level closure recomputation.

Computes ESOs for the host alone, the symbiont alone, and the merged
network, LP-verifies each for self-maintenance (EO), then reports the
complexification metrics: ESO/EO count and size before/after, hybrid-ESO
count and sizes (ESOs of the merged network whose species set intersects
both the host's and the symbiont's own namespaces), how many hybrid ESOs
are ALSO EOs, and the count/mean/max shift relative to the naive union of
the two separate ESO sets.

Usage
-----
    from endosymbiosis_analysis import compare_endosymbiosis

    compare_endosymbiosis(
        host_path="host_alone.txt", host_tag_suffix="_h",
        symbiont_path="symbiont_alone.txt", symbiont_tag_suffix="_s__endo",
        merged_path="merged.txt",
    )

`host_tag_suffix` / `symbiont_tag_suffix` are used only to classify each
species in the merged network as host-side or symbiont-side (for
hybrid-ESO detection) -- they should match how species were actually
named/relabelled.
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))
sys.path.insert(0, os.path.join(_repo_root, 'projects',
                                 'COT_Fundamental_Generators_Exploration', 'scripts'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations import (
    build_rndata, compute_ercs, build_hierarchy,
    compute_synergies_basis_first, compute_complementarities,
)
from pyCOT.analysis.organizations.so_search import compute_elementary_sos
from pyCOT.analysis.organizations.self_maintenance import check_self_maintenance


def compute_eso_summary(path: str, network_id: str, verbose: bool = False):
    """Run the full ERC -> hierarchy -> synergy/complementarity -> ESO
    pipeline on one network file, then LP-verify every ESO for
    self-maintenance (EO). Returns a dict with everything downstream
    analysis needs (rn, rn_data, ercs, per-ESO species sets + is_eo flags,
    species-name lookup)."""
    rn = read_txt(path, exact_names=True)
    rn_data = build_rndata(rn, network_id=network_id)
    ercs = compute_ercs(rn_data, verify=False)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)
    eso_res = compute_elementary_sos(rn_data, ercs, hier, syn, comp, verbose=verbose)

    names = rn_data.species_names
    species_objs = {s.name: s for s in rn.species()}
    eso_species_sets = []
    eso_is_eo = []
    for mask in eso_res.all_elementary_masks:
        full_mask = mask | rn_data.E0_mask
        sp = {names[j] for j in range(rn_data.n_species) if (full_mask >> j) & 1}
        eso_species_sets.append(sp)
        sp_list = [species_objs[n] for n in sp]
        is_eo, _flux, _prod = check_self_maintenance(sp_list, rn)
        eso_is_eo.append(bool(is_eo))

    e0_species = {names[j] for j in range(rn_data.n_species) if (rn_data.E0_mask >> j) & 1}

    return {
        'rn': rn,
        'rn_data': rn_data,
        'ercs': ercs,
        'syn': syn,
        'comp': comp,
        'eso_res': eso_res,
        'eso_species_sets': eso_species_sets,
        'eso_is_eo': eso_is_eo,
        'eo_species_sets': [s for s, is_eo in zip(eso_species_sets, eso_is_eo) if is_eo],
        'e0_species': e0_species,
        'n_species': rn_data.n_species,
        'n_reactions': rn_data.n_reactions,
        'n_ercs': len(ercs),
    }


def _fmt_sizes(sets):
    sizes = sorted(len(s) for s in sets)
    if not sizes:
        return "(none)"
    return f"n={len(sizes)}, sizes={sizes}, max={max(sizes)}, mean={sum(sizes)/len(sizes):.1f}"


def _fmt_eso_eo(summary):
    n_eso = len(summary['eso_species_sets'])
    n_eo = len(summary['eo_species_sets'])
    frac = (n_eo / n_eso * 100) if n_eso else 0.0
    return (f"{_fmt_sizes(summary['eso_species_sets'])}  ->  "
            f"EO: n={n_eo} ({frac:.0f}%), {_fmt_sizes(summary['eo_species_sets'])}")


def compare_endosymbiosis(
    host_path: str,
    symbiont_path: str,
    merged_path: str,
    *,
    symbiont_suffix: str,
    host_suffix: str | None = None,
    verbose: bool = False,
):
    print(f"=== Host alone: {host_path} ===", flush=True)
    host = compute_eso_summary(host_path, "host_alone", verbose=verbose)
    print(f"  species={host['n_species']} reactions={host['n_reactions']} ERCs={host['n_ercs']}")
    print(f"  E0 (food closure) = {sorted(host['e0_species'])}")
    print(f"  ESO: {_fmt_eso_eo(host)}")

    print(f"\n=== Symbiont alone: {symbiont_path} ===", flush=True)
    symb = compute_eso_summary(symbiont_path, "symbiont_alone", verbose=verbose)
    print(f"  species={symb['n_species']} reactions={symb['n_reactions']} ERCs={symb['n_ercs']}")
    print(f"  E0 (food closure) = {sorted(symb['e0_species'])}")
    print(f"  ESO: {_fmt_eso_eo(symb)}")

    print(f"\n=== Merged: {merged_path} ===", flush=True)
    merged = compute_eso_summary(merged_path, "merged", verbose=verbose)
    print(f"  species={merged['n_species']} reactions={merged['n_reactions']} ERCs={merged['n_ercs']}")
    print(f"  E0 (food closure, size {len(merged['e0_species'])}) = {sorted(merged['e0_species'])}")
    print(f"  ESO: {_fmt_eso_eo(merged)}")

    def _side(sp: str) -> str:
        is_symb = sp.endswith(symbiont_suffix)
        if is_symb:
            return 'S'
        # Host species don't necessarily share one distinctive suffix (real
        # BiGG species just keep their original "_c"/"_e" tags) -- anything
        # not carrying the symbiont tag is host-side by construction, since
        # merge_networks.py relabels EVERY symbiont species with that tag
        # and leaves the host untouched. host_suffix, if given, is used only
        # as an extra assertion/sanity check, not as the primary criterion.
        if host_suffix is not None and not sp.endswith(host_suffix):
            return '?'
        return 'H'

    hybrid, host_only, symb_only, other = [], [], [], []
    hybrid_is_eo, host_only_is_eo, symb_only_is_eo = [], [], []
    for s, is_eo in zip(merged['eso_species_sets'], merged['eso_is_eo']):
        sides = {_side(sp) for sp in s}
        if sides == {'H'}:
            host_only.append(s); host_only_is_eo.append(is_eo)
        elif sides == {'S'}:
            symb_only.append(s); symb_only_is_eo.append(is_eo)
        elif 'H' in sides and 'S' in sides:
            hybrid.append(s); hybrid_is_eo.append(is_eo)
        else:
            other.append(s)

    n_hybrid_eo = sum(hybrid_is_eo)
    print(f"\n  Classified merged ESOs:")
    print(f"    host-only:     {_fmt_sizes(host_only)}  (EO: {sum(host_only_is_eo)}/{len(host_only)})")
    print(f"    symbiont-only: {_fmt_sizes(symb_only)}  (EO: {sum(symb_only_is_eo)}/{len(symb_only)})")
    print(f"    HYBRID (both): {_fmt_sizes(hybrid)}  (EO: {n_hybrid_eo}/{len(hybrid)})")
    if other:
        print(f"    unclassified (check suffixes!): {_fmt_sizes(other)}")
        for s in other:
            print(f"      - {sorted(s)}")
    for s, is_eo in zip(hybrid, hybrid_is_eo):
        print(f"    hybrid ESO species (is_EO={is_eo}): {sorted(s)}")

    baseline_sizes = sorted(len(s) for s in host['eso_species_sets']) + \
                      sorted(len(s) for s in symb['eso_species_sets'])
    merged_sizes = sorted(len(s) for s in merged['eso_species_sets'])
    print(f"\n=== Complexification summary ===")
    print(f"  baseline (host ESO ⊔ symbiont ESO, no interface): "
          f"n={len(baseline_sizes)}, max={max(baseline_sizes) if baseline_sizes else 0}, "
          f"mean={sum(baseline_sizes)/len(baseline_sizes):.1f}" if baseline_sizes else "  baseline: (none)")
    print(f"  merged: n={len(merged_sizes)}, max={max(merged_sizes) if merged_sizes else 0}, "
          f"mean={sum(merged_sizes)/len(merged_sizes):.1f}" if merged_sizes else "  merged: (none)")
    if baseline_sizes and merged_sizes:
        print(f"  Delta max size: {max(merged_sizes) - max(baseline_sizes):+d}")
        print(f"  Delta E0 size (host+symbiont E0 vs merged E0, no double count): "
              f"{len(merged['e0_species'])} vs {len(host['e0_species']) + len(symb['e0_species'])} "
              f"({len(merged['e0_species']) - (len(host['e0_species']) + len(symb['e0_species'])):+d})")
    print(f"  Hybrid ESOs found: {len(hybrid)}, of which EO (self-maintaining): {n_hybrid_eo} "
          f"({(n_hybrid_eo/len(hybrid)*100) if hybrid else 0:.0f}%)")

    return {'host': host, 'symb': symb, 'merged': merged, 'hybrid': hybrid,
            'hybrid_is_eo': hybrid_is_eo, 'host_only': host_only, 'symb_only': symb_only}


if __name__ == "__main__":
    _toy = os.path.join(_here, "..", "toy_model")
    compare_endosymbiosis(
        host_path=os.path.join(_toy, "host_alone.txt"),
        symbiont_path=os.path.join(_toy, "symbiont_alone.txt"),
        merged_path=os.path.join(_toy, "merged.txt"),
        host_suffix="_h",
        symbiont_suffix="__endo",
        verbose=False,
    )
