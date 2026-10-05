"""
eso_eo_analysis.py -- generic Elementary Semi-Organization (ESO) to
Elementary Organization (EO) verification, shared by both the E. coli
comparative paper and the COT_Endosymbiosis project.

Terminology (per Dittrich & Speroni di Fenizio 2007, and Centler et al.
2008's own usage): a *semi-organization* is a closed set; an
*organization* is a closed AND self-maintaining set (LP-verified: there
exists a flux vector, strictly positive on every triggered reaction, that
keeps stoichiometric production non-negative for every species). An
*elementary* semi-organization/organization is the order-0 case -- built
directly from ERCs via fundamental synergy/complementarity, without
needing the free-species/latent-join extension that produces higher-order
or spurious structures (see the companion paper's Section 3.1 and
pyCOT.analysis.organizations' module docstring).

What this module was previously called "EPM" (elementary persistent
module) throughout this codebase IS the ESO -- and IS what the reworked
companion paper now calls SO0 (elementary semi-organization): the order-0
output of compute_elementary_sos() is closed by construction (a fixed
point under the network's reactions) but was never itself LP-verified for
self-maintenance -- computing that verification, and reporting the ESO/EO
split, is the entire point of this module. The rename is cosmetic
(aligning with standard COT vocabulary); the underlying computation
(pyCOT.analysis.organizations.so_search.compute_elementary_sos) is
unchanged and, per a direct code audit, already builds every closure by combining PRECOMPUTED
ERC species-masks via the fundamental synergy/complementarity graph
(fundamental_graph.py's FundamentalGraph.extend_state /
erc_syn_close -- pure bitmask Horn-propagation over already-closed ERCs),
never by re-running the species-level reaction-firing closure algorithm
per candidate. That species-level closure algorithm (erc.py's
closure_opt) is used exactly once, to build each ERC's own species_mask
in the first place -- the correct and only place it is needed.

Public API
----------
compute_eso_eo_summary(rn, rn_data, elem_res, *, epsilon=1e-6) -> dict
    Runs the LP self-maintenance check (self_maintenance.check_self_maintenance)
    on every ESO found by compute_elementary_sos, and returns per-ESO records plus
    summary statistics (ESO count, EO count, EO fraction, size stats for
    each).
"""
from __future__ import annotations

import time

from pyCOT.analysis.organizations.self_maintenance import check_self_maintenance


def compute_eso_eo_summary(rn, rn_data, elem_res, *, epsilon: float = 1e-6, verbose: bool = False) -> dict:
    """
    Parameters
    ----------
    rn       : ReactionNetwork  -- the original (non-bitset) network object,
               as returned by pyCOT.io.functions.read_txt. Needed because
               check_self_maintenance takes real Species objects and calls
               rn.sub_reaction_network() internally.
    rn_data  : RNData           -- compiled bitset network (for species_names, E0_mask)
    elem_res : ElementarySOResult -- output of
               pyCOT.analysis.organizations.so_search.compute_elementary_sos

    Returns
    -------
    dict with keys:
      'n_eso'         : int  -- total elementary semi-organizations found
      'n_eo'          : int  -- of those, how many are LP-verified organizations
      'eo_fraction'   : float
      'eso_sizes'     : list[int]  -- species count per ESO (all of them)
      'eo_sizes'      : list[int]  -- species count, EOs only
      'non_eo_sizes'  : list[int]  -- species count, ESOs that are NOT EOs
      'records'       : list[dict] -- per-ESO detail: {mask, size, is_eo, species}
      'lp_time_s'     : float -- total wall time spent in the LP checks
    """
    names = rn_data.species_names
    species_objs = {s.name: s for s in rn.species()}

    records = []
    t0 = time.perf_counter()
    for mask in elem_res.all_elementary_masks:
        full_mask = mask | rn_data.E0_mask
        sp_names = [names[j] for j in range(rn_data.n_species) if (full_mask >> j) & 1]
        sp_list = [species_objs[n] for n in sp_names]
        is_eo, flux, prod = check_self_maintenance(sp_list, rn, epsilon=epsilon)
        records.append({
            'mask': full_mask, 'size': len(sp_names), 'is_eo': bool(is_eo),
            'species': sp_names,
        })
        if verbose:
            print(f"    ESO size={len(sp_names):4d}  is_EO={is_eo}", flush=True)
    lp_time = time.perf_counter() - t0

    eso_sizes = [r['size'] for r in records]
    eo_sizes = [r['size'] for r in records if r['is_eo']]
    non_eo_sizes = [r['size'] for r in records if not r['is_eo']]
    n_eso = len(records)
    n_eo = len(eo_sizes)

    return {
        'n_eso': n_eso,
        'n_eo': n_eo,
        'eo_fraction': (n_eo / n_eso) if n_eso else 0.0,
        'eso_sizes': eso_sizes,
        'eo_sizes': eo_sizes,
        'non_eo_sizes': non_eo_sizes,
        'records': records,
        'lp_time_s': lp_time,
    }


def summary_line(label: str, summary: dict) -> str:
    eo_mean = (sum(summary['eo_sizes']) / len(summary['eo_sizes'])) if summary['eo_sizes'] else 0.0
    eso_mean = (sum(summary['eso_sizes']) / len(summary['eso_sizes'])) if summary['eso_sizes'] else 0.0
    return (f"{label}: {summary['n_eso']} ESO ({eso_mean:.1f} mean sz) -> "
            f"{summary['n_eo']} EO ({summary['eo_fraction']*100:.0f}%, {eo_mean:.1f} mean sz)  "
            f"[LP: {summary['lp_time_s']:.2f}s]")
