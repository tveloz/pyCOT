"""
organizations.py — Verified Organizations from the EPM/ESPM lattice.

This is THE way this library computes organizations. It replaces the old
brute_force_organizations / compute_all_organizations enumeration in
Persistent_Modules_Generator.py, which is exact but exponential in the
number of ERCs (observed: ~52 minutes on Centler et al. 2006's 92-species,
23-ERC starvation network) and was never intended to reach genome-scale
networks (hundreds of ERCs, thousands of reactions).

Pipeline
--------
  1. ERC computation                    erc.compute_ercs
  2. Fundamental relations              synergy.compute_synergies_basis_first,
                                         complementarity.compute_complementarities
  3. Fundamental hierarchy              hierarchy.build_hierarchy
  4. EPM/ESPM exploration               epm.compute_epms, epm.compute_espm
  5. LP self-maintenance verification   self_maintenance.check_self_maintenance

Why step 5 is not optional
---------------------------
Steps 1-4 (the cot_gen engine) find every SEMI-organization: a closed
species set X with req(X) = 0, meaning every species some active reaction
in X needs is produced by *some* active reaction in X. This is fast
(genome-scale in practice) but it is a strictly weaker condition than the
classical COT definition of self-maintenance (Def 2.4): SSM says nothing
about whether a non-negative flux vector actually exists that keeps every
triggered reaction genuinely active while balancing net production against
consumption for every species. A species set can be SSM and still fail to
self-maintain (e.g. two reactions that only balance if run at exactly
opposite, impossible-to-realize rates).

Step 5 promotes each SSM found by the fast engine into a verified
Organization by running the same LP check
(self_maintenance.check_self_maintenance) that
projects/Decomposition_Theorem/decomp/circuits.py already relies on when
decomposing organizations into E/F/fragile-circuit structure. This module
is the single, documented, reusable place that combination lives for
whole-network organization enumeration (decomp applies the same LP kernel
per-fragile-circuit instead, for a different purpose).

Free-species extension (why it's needed even after EPM/ESPM)
---------------------------------------------------------------
EPM/ESPM, like the old brute_force_organizations, builds every semi-
organization out of ERC closures. A species that never appears as the
*sole* trigger of any reaction (no "s -> ..." decay/consumption reaction
exists for it in isolation) can fail to appear in ANY ERC's closure on its
own, even though it could be freely added to an already-valid organization
without disturbing closure or self-maintenance -- there is simply no
reaction to fire, so nothing changes. Confirmed on Centler et al. 2006's
E. coli sugar model: Glcex/Lacex/Glyex have this property (the paper's own
text: "the remaining species that do not decay are... Glcex, Lacex, and
Glyex"), and without this extension step the pipeline finds only 1 of the
4 published organizations per scenario (verified against the paper -- see
projects/COT_Fundamental_Generators_Exploration/scripts/reproduce_centler2006.py for
the original, network-specific diagnosis this generalizes).

This module detects such species automatically (no reaction has them as
its lone reactant) and, for every semi-organization found, tests adding
every subset of them via a cheap structural closure check (does it trigger
anything new? if not, it's a free addition candidate) before handing the
extended candidate to the same LP verification as everything else.

Public API
----------
compute_organizations(rn, *, network_id="", max_espm_order=10,
                       verify_organizations=True, verbose=False,
                       counters=None) -> OrganizationsResult
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

from .io_pyCOT import build_rndata
from .erc import compute_ercs
from .hierarchy import build_hierarchy
from .synergy import compute_synergies_basis_first
from .complementarity import compute_complementarities
from .epm import compute_epms, compute_espm
from .closure import closure_opt
from .self_maintenance import check_self_maintenance


def _free_candidate_species(rn_data) -> list[int]:
    """
    Species (bit indices, quotiented space) that never appear as the SOLE
    reactant of any reaction -- i.e. no unconditional decay/consumption
    reaction exists for them in isolation. These are the only species that
    can ever be "freely" addable to a semi-organization without a per-SO
    closure check ruling them out first (see module docstring).

    E0 species are excluded (they are already implicit in every set).
    """
    has_solo_consumer = set()
    for s in rn_data.supp_q:
        if bin(s).count('1') == 1:
            has_solo_consumer.add(s.bit_length() - 1)

    candidates = []
    for i in range(rn_data.n_species):
        if (rn_data.E0_mask >> i) & 1:
            continue
        if i in has_solo_consumer:
            continue
        candidates.append(i)
    return candidates


def _extend_with_free_species(so_masks: set[int], rn_data, free_candidates: list[int],
                               *, exhaustive_cap: int = 10) -> set[int]:
    """
    For every semi-organization mask in so_masks, try adding subsets of
    free_candidates not already present; keep those that don't trigger any
    new reaction (structural closure check only -- LP verification happens
    later, in compute_organizations' Stage 6, for every mask this returns).

    Two regimes, chosen per base mask by how many free candidates are
    actually missing from it:

      missing <= exhaustive_cap: full power-set enumeration (exact -- this
        is what was validated against Centler et al. 2006, where every
        scenario has <= 3 free candidates).

      missing > exhaustive_cap: reduced search -- test each candidate
        singly, then test the union of every candidate that individually
        passed, in ONE more check. Confirmed necessary, not theoretical:
        e_coli_core (a real, tiny 72-species BiGG network, nothing exotic)
        has 15 free candidates. Full power-set enumeration there produces
        2^15 = 32768 candidate closures PER base semi-organization, which
        measured out to 67158 total extensions, 458 seconds, and 14610
        downstream "organizations" for a 72-species network -- correct
        under the theory (each independently-addable free species really
        does double the count) but not a usable default, and not what
        Centler's own paper's scale ever exercised.

    This is a genuine completeness trade-off above the cap: a "mixed"
    subset where candidate A is only safe to add together with B (neither
    alone, only the pair) would be missed if A or B individually fails the
    singleton test. Singles + the will-conservatively-fail-if-any-single-
    conflicts union check is deliberately biased toward "report the
    common, structurally-clean cases fast" over "enumerate every
    mathematically valid but combinatorially exploding variant" -- flag
    this in reports on networks that hit the reduced regime (see
    OrganizationsResult.stats['n_free_candidates_over_cap']).
    """
    if not free_candidates:
        return set(so_masks)

    inv_idx = list(rn_data.species_to_reactions)
    extended = set(so_masks)
    for base in list(so_masks):
        missing = [c for c in free_candidates if not (base >> c) & 1]
        if not missing:
            continue

        if len(missing) <= exhaustive_cap:
            for r in range(1, len(missing) + 1):
                for combo in itertools.combinations(missing, r):
                    add_mask = 0
                    for c in combo:
                        add_mask |= 1 << c
                    candidate = base | add_mask
                    if candidate in extended:
                        continue
                    closed = closure_opt(rn_data.supp_q, rn_data.prod_q, candidate, inv_idx)
                    if closed == candidate:
                        extended.add(candidate)
        else:
            safe_singles = []
            for c in missing:
                candidate = base | (1 << c)
                if candidate in extended:
                    safe_singles.append(c)
                    continue
                closed = closure_opt(rn_data.supp_q, rn_data.prod_q, candidate, inv_idx)
                if closed == candidate:
                    extended.add(candidate)
                    safe_singles.append(c)
            if len(safe_singles) > 1:
                union_mask = base
                for c in safe_singles:
                    union_mask |= 1 << c
                if union_mask not in extended:
                    closed = closure_opt(rn_data.supp_q, rn_data.prod_q, union_mask, inv_idx)
                    if closed == union_mask:
                        extended.add(union_mask)
    return extended


def _saturate_latent_joins(so_masks: set[int], rn_data, *, max_rounds: int = 20,
                            pairwise_cap: int = 500) -> set[int]:
    """
    Close the found semi-organization set under "latent join": epm.py's
    compute_epms/compute_espm deliberately never search for unions of
    already-SSM sets that share no fundamental synergy/complementarity edge
    (see epm.latent_join's docstring -- it is provided precisely so callers
    needing full enumeration can recover these on demand instead of paying
    for exhaustive latent-join search inside the hot DFS).

    Confirmed necessary, not theoretical: on Centler et al. 2006's
    "all_sugars" scenario, the published top organization (the entire
    92-species network) is exactly such a latent join -- the DFS finds
    several smaller SOs that jointly cover the whole species set once
    unioned, but no single fundamental synergy/complementarity edge
    connects the two largest of them, so compute_espm alone stops short of
    it (found only 78/84/86; the published set is [78,84,86,92]).

    Deliberately does NOT use epm.latent_join's connectivity gate
    (_is_connected). That check exists to certify the DFS's OWN combination
    rules produced a set reachable through genuinely interacting reactions
    -- it is not part of the classical COT definition of an organization
    (Def 2.1: closed + self-maintaining, nothing about topology), and it is
    actively wrong to apply here: it always reports "disconnected" for any
    mask containing E0 species, because supp_q/prod_q have E0 stripped out
    by construction, so E0 species can never appear adjacent to anything in
    that graph regardless of the real network's structure (confirmed while
    debugging the all_sugars case above -- E0 alone made _is_connected
    return False even though the LP self-maintenance check downstream
    correctly accepts the resulting set). Two disjoint, non-interacting
    self-maintaining subsystems coexisting (e.g. an unrelated sugar pathway
    sitting inert next to another) are a legitimate organization under
    Def 2.1 regardless of whether they share a reaction edge -- Stage 5's
    LP check is the actual authority on validity, not this probe.

    Two probes per round:
      (a) grand-union probe (cheap, O(k)): closure of the union of every
          mask currently known -- catches "top of the discovered
          sub-lattice" cases like the one above in a single extra check,
          without any pairwise blow-up.
      (b) pairwise probe (O(k^2)): every pair of currently-known masks.
          Skipped once the known-SO count exceeds pairwise_cap, since it
          is quadratic and this function is called on genome-scale results
          too; the grand-union probe alone still catches the common case.

    Runs to a fixed point (no new mask found) or max_rounds, whichever
    comes first.
    """
    from .closure import build_inv_idx, closure_opt, is_ssm

    inv_idx = build_inv_idx(rn_data.supp_q, rn_data.n_species)

    def _probe(candidate: int) -> int | None:
        closed = closure_opt(rn_data.supp_q, rn_data.prod_q, candidate, inv_idx)
        if is_ssm(rn_data.supp_q, rn_data.prod_q, closed):
            return closed
        return None

    masks = set(so_masks)
    for _ in range(max_rounds):
        new_found = False

        union_all = 0
        for m in masks:
            union_all |= m
        if union_all not in masks:
            verified = _probe(union_all)
            if verified is not None and verified not in masks:
                masks.add(verified)
                new_found = True

        if len(masks) <= pairwise_cap:
            current = list(masks)
            for x, y in itertools.combinations(current, 2):
                j = _probe(x | y)
                if j is not None and j not in masks:
                    masks.add(j)
                    new_found = True

        if not new_found:
            break
    return masks


@dataclass
class SemiOrganization:
    """
    One node of the EPM/ESPM lattice, promoted with an LP verification verdict.

    species_mask  : int              — QUOTIENTED (E0-excluded) species bitmask,
                     in rn_data's index space -- pass this directly to
                     pyCOT.analysis.decomposition.decompose(species_mask,
                     rn_data, S_full) for E/F/fragile-circuit decomposition
    species_names : tuple[str, ...]  — full species set (E0 species included)
    order         : int              — 0 = EPM (elementary), >=1 = ESPM order
    is_organization : bool           — True iff LP self-maintenance verified
    flux          : dict[str, float] | None — reaction_name -> flux, if verified
    """
    species_mask: int
    species_names: tuple
    order: int
    is_organization: bool
    flux: dict | None = None


@dataclass
class OrganizationsResult:
    """
    Full output of compute_organizations.

    semiorganizations : list[SemiOrganization] — every SSM found (EPMs + ESPMs)
    organizations     : list[SemiOrganization] — the subset that is LP-verified
                         self-maintaining (empty if verify_organizations=False)
    rn_data, ercs, hier, syn_result, comp_result, epm_result, espm_result
                      : intermediate pipeline artifacts, kept for reuse
                        (e.g. by the ERC-hierarchy visualization or by
                        projects/Decomposition_Theorem's decomposition step)
    stats             : dict — pipeline timing/counts
    """
    semiorganizations: list
    organizations: list
    rn_data: object
    ercs: list
    hier: object
    syn_result: object
    comp_result: object
    epm_result: object
    espm_result: object
    stats: dict = field(default_factory=dict)


def _species_objects_for_mask(rn, rn_data, quotiented_mask: int):
    """
    Map a quotiented (E0-stripped) species bitmask back to a list of pyCOT
    Species objects for the FULL set (E0 species included -- E0 is always
    part of every closed set, it is just stripped from cot_gen's masks for
    combinatorial efficiency; see cot_types.RNData's docstring).
    """
    full_mask = quotiented_mask | rn_data.E0_mask
    names = set(rn_data.bitset_to_names(full_mask))
    return [sp for sp in rn.species() if sp.name in names]


def compute_organizations(
    rn,
    *,
    network_id: str = "",
    max_espm_order: int = 10,
    verify_organizations: bool = True,
    verbose: bool = False,
    counters=None,
) -> OrganizationsResult:
    """
    Run the full ERC -> fundamental relations -> hierarchy -> EPM/ESPM ->
    LP-verified-organizations pipeline on a pyCOT ReactionNetwork.

    Parameters
    ----------
    rn : pyCOT ReactionNetwork (e.g. from pyCOT.io.functions.read_txt)
    network_id : label carried through metrics/instrumentation
    max_espm_order : safety cap on ESPM BFS depth (epm.compute_espm)
    verify_organizations : if False, skip step 5 (LP check) entirely and
        return every SSM with is_organization=False -- useful when you only
        need the fast semi-organization lattice (e.g. for
        Decomposition_Theorem's own downstream per-circuit LP checks, which
        do their own verification at a finer grain).
    verbose : print per-stage progress (forwarded to compute_epms/compute_espm)
    counters : optional metrics.Counters for instrumentation

    Returns
    -------
    OrganizationsResult
    """
    stats: dict = {}

    if verbose:
        print(f"[organizations] Stage 1: building RNData + ERCs for '{network_id}'...")
    rn_data = build_rndata(rn, network_id=network_id)
    ercs = compute_ercs(rn_data, counters=counters, verify=False)
    stats['n_species'] = rn_data.n_species
    stats['n_reactions'] = rn_data.n_reactions
    stats['n_ercs'] = len(ercs)

    if verbose:
        print(f"[organizations] Stage 2: fundamental hierarchy + relations "
              f"over {len(ercs)} ERCs...")
    hier = build_hierarchy(ercs)
    syn_result = compute_synergies_basis_first(ercs, hier, counters=counters)
    comp_result = compute_complementarities(ercs, hier, syn_result, counters=counters)
    stats['n_fundamental_synergies'] = len(syn_result.fundamental)
    stats['n_fundamental_complementarities'] = len(comp_result.fundamental)

    if verbose:
        print("[organizations] Stage 3: EPM exploration...")
    epm_result = compute_epms(
        rn_data, ercs, hier, syn_result, comp_result,
        counters=counters, verbose=verbose,
    )
    if verbose:
        print("[organizations] Stage 4: ESPM exploration...")
    espm_result = compute_espm(
        rn_data, ercs, hier, syn_result, comp_result, epm_result,
        max_order=max_espm_order, counters=counters, verbose=verbose,
    )
    stats['n_epms'] = len(epm_result.all_epm_masks)
    stats['n_espms'] = espm_result.total_espm()
    stats['max_espm_order_reached'] = espm_result.max_order()

    # Normalize every SO mask by stripping E0 bits before doing anything
    # else. compute_ercs injects E0 itself as a P-ERC (species_mask ==
    # E0_mask, see erc.py's module docstring on "E0 as a P-ERC"), so a DFS
    # state that happens to pull the E0-ERC into its erc_set ends up with
    # E0's bits baked directly into `sp`, while an otherwise-identical
    # state that never touched the E0-ERC does not -- two different
    # integers for the exact same logical organization once E0 is added
    # back for display. Confirmed as a real duplicate-organization bug on
    # Centler's lactose/glycerol/all_sugars scenarios (same species set,
    # counted twice under two different raw masks) before this
    # normalization was added. Stripping here makes every mask this
    # function touches from this point on E0-free and canonical; E0 is
    # added back exactly once, at the point species names are produced.
    e0 = rn_data.E0_mask
    raw_so_masks = {m & ~e0 for m in espm_result.all_so_masks}

    so_order = {}
    for sp in epm_result.all_epm_masks:
        so_order[sp & ~e0] = 0
    for k, masks in espm_result.espm_by_order.items():
        for sp in masks:
            so_order[sp & ~e0] = k

    # ── Stage 4.5: free-species extension (see module docstring) ──────────
    free_candidates = _free_candidate_species(rn_data)
    all_masks = _extend_with_free_species(
        raw_so_masks, rn_data, free_candidates,
    )
    n_extended = len(all_masks) - len(raw_so_masks)
    stats['n_free_candidate_species'] = len(free_candidates)
    stats['n_free_species_extensions'] = n_extended
    stats['n_free_candidates_over_cap'] = len(free_candidates) > 10
    if verbose and n_extended:
        print(f"[organizations] Stage 4.5: free-species extension added "
              f"{n_extended} candidate sets ({len(free_candidates)} "
              f"never-solely-consumed species considered).")

    # ── Stage 4.6: latent-join saturation (see _saturate_latent_joins) ────
    n_before_joins = len(all_masks)
    all_masks = _saturate_latent_joins(all_masks, rn_data)
    n_joined = len(all_masks) - n_before_joins
    stats['n_latent_join_extensions'] = n_joined
    if verbose and n_joined:
        print(f"[organizations] Stage 4.6: latent-join saturation added "
              f"{n_joined} candidate sets not reachable via synergy/"
              f"complementarity DFS alone.")

    # Assign an order to every newly-added mask: 1 + max(order of any proper
    # sub-mask already known), matching epm.py's own order-assignment rule.
    for m in sorted(all_masks - set(so_order), key=lambda x: bin(x).count('1')):
        max_sub = -1
        for sub, sub_ord in so_order.items():
            if sub != m and (sub & m) == sub and sub_ord > max_sub:
                max_sub = sub_ord
        so_order[m] = max_sub + 1

    # ── Stage 5: LP self-maintenance verification ─────────────────────────
    if verbose:
        print(f"[organizations] Stage 5: LP-verifying "
              f"{len(all_masks)} semi-organizations...")

    semiorganizations: list[SemiOrganization] = []
    organizations: list[SemiOrganization] = []
    n_verified_true = 0
    for sp_mask in sorted(all_masks):
        order = so_order.get(sp_mask, 0)
        full_names = tuple(sorted(
            rn_data.bitset_to_names(sp_mask | rn_data.E0_mask)
        ))
        is_org = False
        flux_dict = None
        if verify_organizations:
            species_objs = _species_objects_for_mask(rn, rn_data, sp_mask)
            ok, flux_vec, _prod = check_self_maintenance(species_objs, rn)
            if ok:
                is_org = True
                n_verified_true += 1
                sub_rn = rn.sub_reaction_network(species_objs)
                rxn_names = [r.node.name for r in
                             sorted(sub_rn.reactions(), key=lambda r: r.node.index)]
                flux_dict = {name: float(v) for name, v in zip(rxn_names, flux_vec)}
        so = SemiOrganization(
            species_mask=sp_mask, species_names=full_names, order=order,
            is_organization=is_org, flux=flux_dict,
        )
        semiorganizations.append(so)
        if is_org:
            organizations.append(so)

    stats['n_semiorganizations'] = len(semiorganizations)
    stats['n_organizations'] = n_verified_true

    if verbose and verify_organizations:
        print(f"[organizations] Done: {len(semiorganizations)} semi-organizations, "
              f"{n_verified_true} verified as true organizations.")

    return OrganizationsResult(
        semiorganizations=semiorganizations,
        organizations=organizations,
        rn_data=rn_data, ercs=ercs, hier=hier,
        syn_result=syn_result, comp_result=comp_result,
        epm_result=epm_result, espm_result=espm_result,
        stats=stats,
    )
