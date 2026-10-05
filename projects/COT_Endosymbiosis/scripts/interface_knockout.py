"""
interface_knockout.py -- interface-fragility analysis (THEORY.md Section 2.3,
REPORT.md Section 6b).

For a host/symbiont pair with a curated interface reaction list, builds the
FULL merged network plus one variant per interface reaction with THAT
reaction removed ("knocked out"), computes Elementary Semi-Organizations
(ESOs) for each, and reports:

  - how many of the FULL merged network's hybrid ESOs survive (species set
    still separately reachable / still classified hybrid) vs. disappear
    entirely once a given interface reaction is removed
  - the resulting mean/max hybrid ESO size with that reaction gone

This is the mechanical test of "which specific cross-feeding link is
load-bearing": a reaction whose removal collapses many/large hybrid ESOs
back to non-hybrid (or removes them from the ESO set entirely) is, in this
framework's terms, the load-bearing one -- the computable analogue of "the
gene/transporter most under selective pressure to be retained."

Usage: import run_knockout_analysis and call with the same host/
symbiont/interface arguments used to build the full merger.
"""
from __future__ import annotations

import os
import sys
import time

_here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _here)
from merge_networks import merge_networks
from endosymbiosis_analysis import compute_eso_summary


def _classify(eso_species_sets, host_marker, symbiont_suffix):
    hybrid, host_only, symb_only = [], [], []
    for s in eso_species_sets:
        has_h = any((sp.endswith(host_marker) if host_marker else not sp.endswith(symbiont_suffix))
                    for sp in s)
        has_s = any(sp.endswith(symbiont_suffix) for sp in s)
        if has_h and has_s:
            hybrid.append(s)
        elif has_h:
            host_only.append(s)
        elif has_s:
            symb_only.append(s)
    return hybrid, host_only, symb_only


def run_knockout_analysis(
    host_path: str,
    symbiont_path: str,
    interface_reactions: list[str],
    *,
    symbiont_suffix: str,
    host_marker: str | None = None,
    tmp_dir: str,
    label: str = "case",
):
    print(f"=== Interface knockout analysis: {label} ===", flush=True)
    os.makedirs(tmp_dir, exist_ok=True)

    full_path = merge_networks(
        host_path=host_path, symbiont_path=symbiont_path,
        interface_reactions=interface_reactions,
        out_path=os.path.join(tmp_dir, f"{label}_full.txt"),
    )
    t0 = time.perf_counter()
    full = compute_eso_summary(full_path, f"{label}_full", verbose=False)
    full_hybrid, full_host, full_symb = _classify(full['eso_species_sets'], host_marker, symbiont_suffix)
    full_hybrid_sizes = sorted(len(s) for s in full_hybrid)
    print(f"[full interface] {len(interface_reactions)} reactions, "
          f"{len(full['eso_species_sets'])} ESOs total, {len(full_hybrid)} hybrid "
          f"(sizes {full_hybrid_sizes[:3]}...{full_hybrid_sizes[-3:] if len(full_hybrid_sizes)>3 else ''}, "
          f"mean={sum(full_hybrid_sizes)/max(1,len(full_hybrid_sizes)):.1f}) "
          f"[{time.perf_counter()-t0:.1f}s]", flush=True)

    results = []
    for i in range(len(interface_reactions)):
        remaining = interface_reactions[:i] + interface_reactions[i+1:]
        removed = interface_reactions[i]
        ko_path = merge_networks(
            host_path=host_path, symbiont_path=symbiont_path,
            interface_reactions=remaining,
            out_path=os.path.join(tmp_dir, f"{label}_ko{i}.txt"),
        )
        t0 = time.perf_counter()
        ko = compute_eso_summary(ko_path, f"{label}_ko{i}", verbose=False)
        ko_hybrid, ko_host, ko_symb = _classify(ko['eso_species_sets'], host_marker, symbiont_suffix)
        ko_hybrid_sizes = sorted(len(s) for s in ko_hybrid)
        elapsed = time.perf_counter() - t0
        mean_hyb = sum(ko_hybrid_sizes) / max(1, len(ko_hybrid_sizes))
        print(f"[knock out: {removed.split(':')[0]}] remaining={len(remaining)} reactions, "
              f"{len(ko['eso_species_sets'])} ESOs total, {len(ko_hybrid)} hybrid "
              f"(mean={mean_hyb:.1f}, max={max(ko_hybrid_sizes) if ko_hybrid_sizes else 0}) "
              f"-- vs full: {len(full_hybrid)} hybrid (mean={sum(full_hybrid_sizes)/max(1,len(full_hybrid_sizes)):.1f}) "
              f"[{elapsed:.1f}s]", flush=True)
        results.append({
            'removed': removed, 'n_eso': len(ko['eso_species_sets']),
            'n_hybrid': len(ko_hybrid), 'hybrid_mean': mean_hyb,
            'hybrid_max': max(ko_hybrid_sizes) if ko_hybrid_sizes else 0,
        })

    print(f"\n=== Summary: {label} ===")
    print(f"{'removed reaction':<40} {'hybrid ESOs':>12} {'hybrid mean sz':>15} {'hybrid max sz':>14}")
    print(f"{'(none -- full interface)':<40} {len(full_hybrid):>12} "
          f"{sum(full_hybrid_sizes)/max(1,len(full_hybrid_sizes)):>15.1f} "
          f"{max(full_hybrid_sizes) if full_hybrid_sizes else 0:>14}")
    for r in results:
        name = r['removed'].split(':')[0]
        print(f"{name:<40} {r['n_hybrid']:>12} {r['hybrid_mean']:>15.1f} {r['hybrid_max']:>14}")

    return full, results
