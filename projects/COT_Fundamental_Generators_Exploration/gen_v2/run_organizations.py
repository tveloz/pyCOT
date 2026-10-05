"""
run_organizations.py -- compute verified ORGANIZATIONS (not just
semi-organizations) for a single, user-supplied reaction network, using
gen_v2's corrected SO-search engine.

Why this script exists
-----------------------
gen_v2.engine.explore() computes the full semi-organization (SO) lattice --
elementary SOs (order 0) and every higher-order SO built from them -- more
correctly than the production engine (see gen_v2/__init__.py: a real bug in
the production generator-vs-closure logic was traced and fixed here; on at
least one real network the production engine silently missed hundreds of
genuine SOs). But explore() does NOT do the LP self-maintenance check that
turns a semi-organization (req==0 bookkeeping only) into a true organization
(Def. 2.4: an actual non-negative flux vector exists, strictly positive on
every triggered reaction). That check already exists --
self_maintenance.check_self_maintenance, the same LP kernel the production
organizations.py uses -- it's just not wired up to gen_v2 yet. This script
adds exactly that: run gen_v2's SO search, then LP-verify every SO found.

Status caveat (read before trusting this on something important)
------------------------------------------------------------------
gen_v2 is explicitly experimental. It has been validated two ways so far:
  1. validate_small.py -- exact match against the brute-force oracle AND
     the production engine, on the gold networks + a handful of small real
     BioModels networks. This passes.
  2. compare_old_new.py / run_large.py -- old-vs-new agreement on larger
     (80-1500 reaction) real networks, where no brute-force oracle is
     feasible. This sweep is the thing currently running on your other
     machine and has NOT finished (as of the last check, every network in
     it was still "partial", none "complete").
So: for a small-to-moderate network (roughly under 100-200 reactions, no
unusual hierarchy shape) gen_v2 and the production engine have so far always
agreed, and gen_v2 is the more defensible of the two (it's the one with the
known bug fixed). For anything larger, or if you want maximum confidence,
also run the SAME network through the production pipeline
(pyCOT.analysis.organizations.compute_organizations) and diff the results --
see cross_check_with_production() below, on by default.

HOW TO RUN
----------
  1. Edit NETWORK below to point at your network's .txt file (pyCOT reaction
     format: "name: reactants => products;" one per line, "=>" with nothing
     on the left for inflow).
  2. python projects/COT_Fundamental_Generators_Exploration/gen_v2/run_organizations.py
"""
from __future__ import annotations

import os
import sys
import time

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations.io_pyCOT import build_rndata
from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity import compute_complementarities
from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph
from pyCOT.analysis.organizations.self_maintenance import check_self_maintenance

from gen_v2.engine import explore

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — edit this block                                        ║
# ╚══════════════════════════════════════════════════════════════════════════╝

#NETWORK = "C:/Users/tvelo/Dropbox/Public/AcademicWork/Europe/CLEA\Postdocs/TempletonPostdoc/sftw/pyCOT/data/Ecological_models/MLM_Afidos/RN_MLM26_V02.txt"   # absolute path, or relative to repo root
NETWORK= "data/biochemical_databases/Other/central_ecoli.txt"
#NETWORK= "data/Examples_tests/autopoietic.txt"
NETWORK= "data\\biochemical_databases\\BiGG\\bigg_e_coli_core.txt"
NETWORK= "data\\biochemical_databases\\BioMD_apoptosis\\BIOMD0000000407.txt"
NETWORK= "data\\biochemical_databases\\BioMD_cell_cycle\\BIOMD0000000056.txt"
NETWORK= "data\\biochemical_databases\\BiGG\\bigg_iAB_RBC_283.txt"
NETWORK_ID = "Test"                         # label for output; blank = derive from filename

VERIFY_ORGANIZATIONS = True             # False = report SOs only, skip the LP stage
CROSS_CHECK_WITH_PRODUCTION = True      # also run the old engine and diff -- see caveat above
TIME_BUDGET_S = None                    # seconds, or None = no limit (small/moderate networks)

PRINT_FUNDAMENTAL_HIERARCHY = True      # print ERCs + fundamental synergies/complementarities
PRINT_MAX_ROWS = 200                    # per section, so genome-scale networks don't flood the console
PLOT_FUNDAMENTAL_HIERARCHY = True       # ERC hierarchy + fundamental relations + level-pair stats panel

PLOT_HIERARCHY = True                   # render the Hasse diagram (semi-orgs vs. organizations)
HIGHLIGHT_INDEX = None                  # int, or None -- pick the [idx] shown in the printed
                                         # listing below to highlight ONE node in the diagram
                                         # (run once with None to see the indices, then rerun)

RESULTS_CSV = os.path.join(_here, "..", "outputs", "gen_v2_organizations", "results.csv")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  Script body                                                            ║
# ╚══════════════════════════════════════════════════════════════════════════╝


def _species_for_mask(rn, rn_data, quotiented_mask: int):
    """Map a quotiented (E0-stripped) species bitmask back to real pyCOT
    Species objects for the FULL set (E0 is always implicitly present)."""
    full_mask = quotiented_mask | rn_data.E0_mask
    names = set(rn_data.bitset_to_names(full_mask))
    return [sp for sp in rn.species() if sp.name in names]


def compute_organizations_gen_v2(
    rn, rn_data, ercs, hier, syn, comp,
    *, verify: bool = True, time_budget_s: float | None = None, verbose: bool = True,
):
    """
    Run gen_v2's SO search, then LP-verify every SO found.

    Returns
    -------
    dict with:
      'so_result'      : gen_v2.engine.ExploreResult
      'semiorganizations' : list[dict]  -- one per SO: {mask, names, order, is_organization, flux}
      'organizations'      : list[dict]  -- the LP-verified subset
    """
    g = FundamentalGraph(ercs, hier, syn, comp)
    if verbose:
        print(f"  [gen_v2] exploring the SO lattice ({len(ercs)} ERCs)...")
    t0 = time.perf_counter()
    so_result = explore(g, verbose=verbose, time_budget_s=time_budget_s)
    t_explore = time.perf_counter() - t0
    if verbose:
        print(f"  [gen_v2] done in {t_explore:.2f}s -- "
              f"{len(so_result.all_so_masks)} semi-organizations "
              f"({len(so_result.elementary_masks)} elementary), "
              f"complete={so_result.complete}")

    semiorgs: list[dict] = []
    orgs: list[dict] = []
    e0 = rn_data.E0_mask
    so_order = {sp: k for k, masks in so_result.so_by_order.items() for sp in masks}

    if verbose and verify:
        print(f"  [gen_v2] LP-verifying {len(so_result.all_so_masks)} semi-organizations...")
    t0 = time.perf_counter()
    for mask in so_result.all_so_masks:
        full_names = tuple(sorted(rn_data.bitset_to_names(mask | e0)))
        order = so_order.get(mask, 0)
        is_org = False
        flux = None
        if verify:
            sp_list = _species_for_mask(rn, rn_data, mask)
            ok, flux_vec, _prod = check_self_maintenance(sp_list, rn)
            if ok:
                is_org = True
                sub_rn = rn.sub_reaction_network(sp_list)
                rxn_names = [r.node.name for r in
                             sorted(sub_rn.reactions(), key=lambda r: r.node.index)]
                flux = {name: float(v) for name, v in zip(rxn_names, flux_vec)}
        record = {'mask': mask, 'names': full_names, 'order': order,
                  'is_organization': is_org, 'flux': flux}
        semiorgs.append(record)
        if is_org:
            orgs.append(record)
    t_lp = time.perf_counter() - t0
    if verbose and verify:
        print(f"  [gen_v2] LP stage done in {t_lp:.2f}s -- "
              f"{len(orgs)}/{len(semiorgs)} are verified organizations")

    return {'so_result': so_result, 'semiorganizations': semiorgs, 'organizations': orgs}


def cross_check_with_production(rn_data, ercs, hier, syn, comp, gen_v2_so_masks: set[int], *, verbose=True):
    """
    Run the SAME network through the production engine
    (pyCOT.analysis.organizations.so_search) and diff its full SO set
    against gen_v2's. Any mismatch here means this network's hierarchy
    shape is one where the two engines currently disagree -- worth
    reporting back rather than trusting either silently.
    """
    from pyCOT.analysis.organizations.so_search import compute_elementary_sos, compute_so_hierarchy

    elem = compute_elementary_sos(rn_data, ercs, hier, syn, comp, verbose=False)
    hr = compute_so_hierarchy(rn_data, ercs, hier, syn, comp, elem,
                               use_vertical_lift=True, max_order=15, verbose=False)
    prod_masks = set(hr.all_so_masks)

    e0 = rn_data.E0_mask
    gen_v2_nz = {m for m in gen_v2_so_masks if m != e0}
    prod_nz = {m for m in prod_masks if m != e0}

    if gen_v2_nz == prod_nz:
        if verbose:
            print(f"  [cross-check] gen_v2 and the production engine AGREE "
                  f"({len(gen_v2_nz)} semi-organizations).")
        return True

    only_gen_v2 = gen_v2_nz - prod_nz
    only_prod = prod_nz - gen_v2_nz
    print(f"  [cross-check] *** MISMATCH *** -- gen_v2={len(gen_v2_nz)}  production={len(prod_nz)}")
    print(f"    only in gen_v2 ({len(only_gen_v2)}): {[hex(m) for m in sorted(only_gen_v2)][:10]}")
    print(f"    only in production ({len(only_prod)}): {[hex(m) for m in sorted(only_prod)][:10]}")
    print(f"    (gen_v2 is the one with the known production-engine bug fixed -- "
          f"'only in production' entries are the ones to be suspicious of, not "
          f"the other way around; but please report this network, it's new data "
          f"for the validation sweep.)")
    return False


def find_e0_index(ercs, rn_data):
    """Index of the inflow ERC E0 (species_mask == E0_mask), or None if no inflow."""
    if not rn_data.E0_mask:
        return None
    return next((i for i, e in enumerate(ercs) if e.species_mask == rn_data.E0_mask), None)


def print_fundamental_hierarchy(ercs, hier, syn, comp, rn_data, stats, *, max_rows: int = 200):
    """
    Print the fundamental hierarchy: every ERC (level, size, full species set
    incl. E0, direct parents), every fundamental synergy and complementarity
    (with the levels of the ERCs involved), and the counts per level triad
    (synergies, target included) and per level pair (complementarities).

    Species sets are shown with the inflow closure E0 added back
    (ERCData.species_mask is E0-stripped), matching the SO listing.
    """
    from pyCOT.visualization.erc_hierarchy_viz import relation_level_counts, syn_triad_label

    e0 = rn_data.E0_mask
    lv = stats.levels

    def _cap(items, render):
        for row in items[:max_rows]:
            print(render(row))
        if len(items) > max_rows:
            print(f"    ... {len(items) - max_rows} more (raise PRINT_MAX_ROWS to see all)")

    print("\n" + "=" * 70)
    print(f"FUNDAMENTAL HIERARCHY  ({stats.n_ercs} ERCs, {stats.n_levels} levels, "
          f"{stats.n_fundamental_syn} fundamental synergies, "
          f"{stats.n_fundamental_comp} fundamental complementarities)")
    print("=" * 70)

    print("\nERCs (level 0 = bottom; E0 = inflow ERC, contained in every ERC)")
    order = sorted(range(len(ercs)), key=lambda i: (lv[i], i))
    def _erc_row(i):
        names = rn_data.bitset_to_names(ercs[i].species_mask | e0)
        parents = ", ".join(f"E{p}" for p in hier.parents[i]) or "-"
        tag = " (inflow)" if ercs[i].species_mask == e0 and e0 else ""
        persist = "persistent" if ercs[i].is_persistent() else "          "
        return (f"  E{i:<4}{tag} L{lv[i]:<3} size={len(names):<4} {persist}  "
                f"parents=[{parents}]  {tuple(names)}")
    _cap(order, _erc_row)

    print("\nFundamental synergies   E_i [L] + E_j [L]  ->  E_k [L]")
    syn_rows = sorted(syn.fundamental, key=lambda t: (min(lv[t.i], lv[t.j]), max(lv[t.i], lv[t.j]), t.i, t.j, t.k))
    _cap(syn_rows, lambda t: f"  E{t.i} [L{lv[t.i]}] + E{t.j} [L{lv[t.j]}]  ->  E{t.k} [L{lv[t.k]}]")

    print("\nFundamental complementarities   producer [L] --species--> consumer [L]"
          "   (chain = hierarchy-comparable pair)")
    comp_rows = sorted(comp.fundamental, key=lambda f: (min(lv[f.prod_idx], lv[f.cons_idx]),
                                                         max(lv[f.prod_idx], lv[f.cons_idx]),
                                                         f.prod_idx, f.cons_idx))
    _cap(comp_rows, lambda f: (f"  E{f.prod_idx} [L{lv[f.prod_idx]}] --{rn_data.species_names[f.species]}--> "
                               f"E{f.cons_idx} [L{lv[f.cons_idx]}]" + ("   (chain)" if f.chain else "")))

    counts = relation_level_counts(stats, syn, comp)
    print("\nFundamental synergies per level triad   (contributor + contributor -> target)")
    for key, v in counts["syn"].items():
        print(f"  {syn_triad_label(key):<12}{v:>6}")
    print("\nFundamental complementarities per level pair (unordered)")
    for (a, b), v in counts["comp"].items():
        print(f"  {f'({a},{b})':<12}{v:>6}")


def plot_fundamental_hierarchy(ercs, hier, syn, comp, rn_data, stats, net_id: str, e0_index):
    """
    Render the ERC hierarchy with fundamental synergies (orange diamonds) and
    complementarities (blue dashed) overlaid, plus a stats panel (summary
    counts, histogram + table of relations per level pair), via pyCOT's
    erc_hierarchy_viz. Saved under outputs/gen_v2_organizations/ and opened
    in the default browser.
    """
    import webbrowser
    from pyCOT.visualization.erc_hierarchy_viz import plot_erc_hierarchy

    out_path = os.path.join(os.path.dirname(RESULTS_CSV), f"{net_id}_fundamental_hierarchy.html")
    plot_erc_hierarchy(
        ercs, hier, syn, comp, out_path,
        stats=stats, e0_index=e0_index, rn_data=rn_data,
        stats_panel=True, title=f"{net_id} -- fundamental hierarchy",
    )
    print(f"  [plot] fundamental hierarchy -> {os.path.normpath(out_path)}")
    webbrowser.open("file://" + os.path.abspath(out_path))


def plot_organization_hierarchy(
    sorted_semiorgs: list[dict], rn_data, net_id: str,
    *, highlight_index: int | None = None, filename: str | None = None,
):
    """
    Render the full semi-organization lattice as an interactive Hasse
    diagram, reusing pyCOT's existing containment-hierarchy visualizer
    (pyCOT.visualization.rn_visualize.hierarchy_visualize_html -- the same
    hierarchical pyvis Hasse-diagram renderer used elsewhere in pyCOT, not
    reimplemented here). Nodes are the discovered semi-organizations;
    cover edges are direct containment (no intermediate semi-organization
    in between), computed from the species sets themselves.

    Color coding
    ------------
      cyan        : semi-organization only (closed + SSM, not LP-verified
                    self-maintaining)
      limegreen   : verified organization (LP self-maintenance holds)
      gold        : the one node at `highlight_index`, if given -- overrides
                    whichever of the above it would otherwise be

    `sorted_semiorgs` must be in the same order used to print the indexed
    listing, so `highlight_index` picks the node the user actually saw.

    Opens automatically in the default browser; saved under
    visualizations/hierarchy_visualize_html/<filename> (relative to CWD --
    pyCOT's existing convention, unchanged here).
    """
    from pyCOT.visualization.rn_visualize import hierarchy_visualize_html

    # Use the FULL species sets (rec['names'] = mask | E0), not the
    # E0-stripped rec['mask']: the inflow closure E0 belongs to every
    # semi-organization, and the E0-only bottom node must be a subset of
    # every other node for the containment edges to be drawn.
    node_sets = [frozenset(rec['names']) for rec in sorted_semiorgs]

    org_sets = [s for s, rec in zip(node_sets, sorted_semiorgs) if rec['is_organization']]
    semiorg_only_sets = [s for s, rec in zip(node_sets, sorted_semiorgs) if not rec['is_organization']]
    lst_color_subsets = []
    if semiorg_only_sets:
        lst_color_subsets.append(("cyan", semiorg_only_sets))
    if org_sets:
        lst_color_subsets.append(("limegreen", org_sets))

    if highlight_index is not None:
        if not (0 <= highlight_index < len(node_sets)):
            print(f"  [plot] HIGHLIGHT_INDEX={highlight_index} out of range "
                  f"(0..{len(node_sets)-1}) -- ignoring.")
        else:
            chosen = node_sets[highlight_index]
            lst_color_subsets = [
                (color, [s for s in subsets if s != chosen])
                for color, subsets in lst_color_subsets
            ]
            lst_color_subsets.append(("gold", [chosen]))

    fname = filename or f"{net_id}_organization_hierarchy.html"
    print(f"  [plot] legend: cyan=semi-organization only, limegreen=verified organization"
          + (", gold=chosen (HIGHLIGHT_INDEX)" if highlight_index is not None else ""))
    hierarchy_visualize_html(
        node_sets,
        lst_color_subsets=lst_color_subsets,
        node_color="lightgray",   # should never actually be used -- every node is colored above
        filename=fname,
    )


def main():
    if not os.path.isabs(NETWORK):
        net_path = os.path.join(_repo, NETWORK)
    else:
        net_path = NETWORK
    if not os.path.exists(net_path):
        raise SystemExit(f"Network file not found: {net_path}\n"
                          f"Edit NETWORK at the top of this script.")

    net_id = NETWORK_ID or os.path.splitext(os.path.basename(net_path))[0]

    print("=" * 70)
    print(f"gen_v2 organization computation -- {net_id}")
    print("=" * 70)

    rn = read_txt(net_path)
    rn_data = build_rndata(rn, network_id=net_id)
    print(f"  species={rn_data.n_species}  reactions={rn_data.n_reactions}  "
          f"E0(inflow)={bin(rn_data.E0_mask).count('1')}")

    ercs = compute_ercs(rn_data, verify=False)
    if not ercs:
        print("  No ERCs found -- nothing to search.")
        return
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)
    print(f"  ERCs={len(ercs)}  fundamental synergies={len(syn.fundamental)}  "
          f"fundamental complementarities={len(comp.fundamental)}")

    if PRINT_FUNDAMENTAL_HIERARCHY or PLOT_FUNDAMENTAL_HIERARCHY:
        from pyCOT.visualization.erc_hierarchy_viz import compute_hierarchy_stats
        e0_index = find_e0_index(ercs, rn_data)
        hstats = compute_hierarchy_stats(ercs, hier, syn, comp, e0_index=e0_index)
        if PRINT_FUNDAMENTAL_HIERARCHY:
            print_fundamental_hierarchy(ercs, hier, syn, comp, rn_data, hstats, max_rows=PRINT_MAX_ROWS)
        if PLOT_FUNDAMENTAL_HIERARCHY:
            plot_fundamental_hierarchy(ercs, hier, syn, comp, rn_data, hstats, net_id, e0_index)

    result = compute_organizations_gen_v2(
        rn, rn_data, ercs, hier, syn, comp,
        verify=VERIFY_ORGANIZATIONS, time_budget_s=TIME_BUDGET_S,
    )

    if CROSS_CHECK_WITH_PRODUCTION:
        cross_check_with_production(
            rn_data, ercs, hier, syn, comp,
            set(result['so_result'].all_so_masks),
        )

    # Single stable index space (0..N-1, sorted by order then size) used
    # both for this listing and for HIGHLIGHT_INDEX below -- so an index
    # the user sees here always means the same node in the diagram.
    all_sorted = sorted(result['semiorganizations'], key=lambda r: (r['order'], len(r['names'])))

    print("\n" + "=" * 70)
    print(f"ALL SEMI-ORGANIZATIONS  (ORG = LP-verified organization, sso = semi-org only)"
          if VERIFY_ORGANIZATIONS else "ALL SEMI-ORGANIZATIONS (not LP-verified -- VERIFY_ORGANIZATIONS=False)")
    print("=" * 70)
    for idx, rec in enumerate(all_sorted):
        tag = "ORG" if rec['is_organization'] else ("sso" if VERIFY_ORGANIZATIONS else "?  ")
        print(f"  [{idx:>3}] {tag}  order={rec['order']:<3}  size={len(rec['names']):<4}  {rec['names']}")
    n_org = sum(1 for r in all_sorted if r['is_organization'])
    print(f"\nTotal: {len(all_sorted)} semi-organizations"
          + (f", {n_org} verified as true organizations" if VERIFY_ORGANIZATIONS else ""))

    if PLOT_HIERARCHY:
        plot_organization_hierarchy(all_sorted, rn_data, net_id, highlight_index=HIGHLIGHT_INDEX)

    # ── Save CSV ──────────────────────────────────────────────────────────
    import csv
    os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)
    with open(RESULTS_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["network", "order", "n_species", "is_organization", "species"])
        for rec in result['semiorganizations']:
            writer.writerow([net_id, rec['order'], len(rec['names']),
                              rec['is_organization'], ";".join(rec['names'])])
    print(f"\nResults saved -> {RESULTS_CSV}")


if __name__ == "__main__":
    main()
