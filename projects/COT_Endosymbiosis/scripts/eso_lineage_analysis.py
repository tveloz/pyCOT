"""
eso_lineage_analysis.py -- traces each HYBRID ESO of a merged host/symbiont
network back to its standalone origins, asking: is this hybrid ESO built
from pieces that were *already* elementary on their own (host part exactly
matches a host-alone ESO, symbiont part exactly matches a symbiont-alone
ESO -- the interface just glued two pre-existing modules together
unchanged), or is it a genuinely NEW joint structure that does not
decompose into any pre-existing standalone piece?

This is the natural next question after endosymbiosis_analysis.py's
host/symbiont/hybrid classification: hybridization (crossing the namespace
boundary) is necessary but not sufficient evidence of structural novelty --
gluing two untouched, already-elementary modules together via a
transporter is a much weaker claim than a joint reorganization that
produces a module neither parent could produce even after closure.

Classification (per hybrid ESO E of the merged network, split into its
host part H = E minus symbiont-tagged species, and symbiont part S = E's
symbiont-tagged species with the tag stripped):

  - "conserved"       : H equals some host-alone ESO AND S equals some
                         symbiont-alone ESO (exact set equality). The
                         interface additively joined two pre-existing,
                         independently-elementary modules.
  - "host-emergent"    : S matches a symbiont-alone ESO exactly, but H does
                         not match any host-alone ESO -- the host side was
                         restructured by the merger.
  - "symbiont-emergent": H matches a host-alone ESO exactly, but S does not
                         match any symbiont-alone ESO.
  - "fully-emergent"   : neither H nor S matches any standalone ESO --
                         genuinely new joint structure, not decomposable
                         into unchanged pre-existing pieces.

For each category we also report the EO fraction (of the hybrid ESOs in
that category, how many are self-maintaining), since the biologically
interesting question is not just "is this a new structure" but "did that
new structure land in a class that could actually persist on its own
merits, or does it require the joint network's own extra machinery/flux
that unaided reachability doesn't capture" -- EO status still requires
literal self-maintenance of the merged system, so this is a genuine,
independently-computed second axis, not a restatement of "emergent".

Usage: import run_lineage_analysis and call with the same host/symbiont/
merged/suffix arguments used by endosymbiosis_analysis.compare_endosymbiosis.
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
from endosymbiosis_analysis import compute_eso_summary


def _closest_jaccard(target: set, candidates: list[set]) -> tuple[float, set | None]:
    best_j, best_c = 0.0, None
    for c in candidates:
        u = len(target | c)
        j = len(target & c) / u if u else 1.0
        if j > best_j:
            best_j, best_c = j, c
    return best_j, best_c


def run_lineage_analysis(
    host_path: str,
    symbiont_path: str,
    merged_path: str,
    *,
    symbiont_suffix: str,
    label: str = "case",
    verbose: bool = False,
):
    print(f"=== ESO lineage analysis: {label} ===", flush=True)
    host = compute_eso_summary(host_path, f"{label}_host", verbose=verbose)
    symb = compute_eso_summary(symbiont_path, f"{label}_symb", verbose=verbose)
    merged = compute_eso_summary(merged_path, f"{label}_merged", verbose=verbose)

    def is_symb(sp: str) -> bool:
        return sp.endswith(symbiont_suffix)

    n = len(symbiont_suffix)
    counts = {'conserved': 0, 'host-emergent': 0, 'symbiont-emergent': 0, 'fully-emergent': 0}
    eo_counts = {'conserved': 0, 'host-emergent': 0, 'symbiont-emergent': 0, 'fully-emergent': 0}
    records = []

    for sp_set, eso_is_eo in zip(merged['eso_species_sets'], merged['eso_is_eo']):
        sides = {is_symb(sp) for sp in sp_set}
        if sides != {True, False}:
            continue  # not hybrid; lineage question doesn't apply
        H = {sp for sp in sp_set if not is_symb(sp)}
        S_raw = {sp for sp in sp_set if is_symb(sp)}
        S = {sp[:-n] if n else sp for sp in S_raw}

        h_match = H in host['eso_species_sets']
        s_match = S in symb['eso_species_sets']

        if h_match and s_match:
            cat = 'conserved'
        elif h_match and not s_match:
            cat = 'symbiont-emergent'
        elif s_match and not h_match:
            cat = 'host-emergent'
        else:
            cat = 'fully-emergent'

        counts[cat] += 1
        if eso_is_eo:
            eo_counts[cat] += 1

        h_j, _ = (1.0, H) if h_match else _closest_jaccard(H, host['eso_species_sets'])
        s_j, _ = (1.0, S) if s_match else _closest_jaccard(S, symb['eso_species_sets'])
        records.append({
            'category': cat, 'is_eo': eso_is_eo, 'size': len(sp_set),
            'host_part_size': len(H), 'symb_part_size': len(S),
            'host_closest_jaccard': round(h_j, 2), 'symb_closest_jaccard': round(s_j, 2),
        })

    total = sum(counts.values())
    print(f"  host ESOs: {len(host['eso_species_sets'])}, symbiont ESOs: {len(symb['eso_species_sets'])}, "
          f"hybrid ESOs analysed: {total}", flush=True)
    print(f"  {'category':<20} {'count':>6} {'%':>6} {'EO':>6} {'EO%':>6}")
    for cat in ['conserved', 'host-emergent', 'symbiont-emergent', 'fully-emergent']:
        c, e = counts[cat], eo_counts[cat]
        pct = (c / total * 100) if total else 0.0
        epct = (e / c * 100) if c else 0.0
        print(f"  {cat:<20} {c:>6} {pct:>5.0f}% {e:>6} {epct:>5.0f}%", flush=True)

    if verbose:
        near_misses = [r for r in records if r['category'] != 'conserved'
                        and max(r['host_closest_jaccard'], r['symb_closest_jaccard']) > 0]
        print(f"\n  Non-conserved records with partial similarity to a standalone ESO "
              f"(near-miss diagnostic, not a match):")
        for r in sorted(near_misses, key=lambda r: -max(r['host_closest_jaccard'], r['symb_closest_jaccard']))[:10]:
            print(f"    {r}")

    return {'host': host, 'symb': symb, 'merged': merged,
            'counts': counts, 'eo_counts': eo_counts, 'records': records}


if __name__ == "__main__":
    _root = os.path.join(_here, '..')
    _rd = os.path.join(_root, 'real_data')
    _toy = os.path.join(_root, 'toy_model')

    run_lineage_analysis(
        host_path=os.path.join(_toy, 'host_alone.txt'),
        symbiont_path=os.path.join(_toy, 'symbiont_alone.txt'),
        merged_path=os.path.join(_toy, 'merged.txt'),
        symbiont_suffix="__endo", label="toy", verbose=True,
    )
    run_lineage_analysis(
        host_path=os.path.join(_rd, 'mito_host_fermentative.txt'),
        symbiont_path=os.path.join(_rd, 'mito_symbiont_aerobic.txt'),
        merged_path=os.path.join(_rd, 'mito_merged.txt'),
        symbiont_suffix="__endo", label="mito", verbose=True,
    )
    run_lineage_analysis(
        host_path=os.path.join(_rd, 'chloro_host_yeast.txt'),
        symbiont_path=os.path.join(_rd, 'chloro_symbiont_synecho_dark.txt'),
        merged_path=os.path.join(_rd, 'chloro_merged_dark.txt'),
        symbiont_suffix="__endo", label="chloro_dark", verbose=True,
    )
    run_lineage_analysis(
        host_path=os.path.join(_rd, 'chloro_host_yeast.txt'),
        symbiont_path=os.path.join(_rd, 'chloro_symbiont_synecho_light.txt'),
        merged_path=os.path.join(_rd, 'chloro_merged_light.txt'),
        symbiont_suffix="__endo", label="chloro_light", verbose=True,
    )
