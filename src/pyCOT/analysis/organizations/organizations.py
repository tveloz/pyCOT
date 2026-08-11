"""
organizations.py — Fundamental Organizations from the EPM/ESPM lattice.

Terminology fixed by Veloz & Bassi, "Synergy and Complementarity: The
Generative Basis of Chemical Organizations" (the paper this whole engine
implements — see src/pyCOT/analysis/organizations/, ported from
projects/COT_Fundamental_Generators_Exploration/cot_gen):

  FUNDAMENTAL organizations : verified organizations reachable through the
    fundamental-generator DFS alone (epm.compute_epms / compute_espm),
    which combines ERCs ONLY via fundamental synergy and fundamental
    complementarity. This is Theorem thm:fund_gen's object: every
    "productively novel" semi-organization has a fundamental generator,
    and nothing outside this set does. THIS IS THE DEFAULT AND ONLY THING
    compute_organizations() COMPUTES UNLESS include_spurious=True IS
    PASSED EXPLICITLY.

  SPURIOUS organizations : every other closed+self-maintaining species set
    -- reachable only by operations the paper explicitly identifies as
    generatively irrelevant: adding a species that triggers no new
    reaction (Sec. 3.2: "every time we add a species s to a persistent
    module X so that it does not participate in any reaction, we obtain a
    different persistent module because neither closure or self-
    maintainance is affected... building non-reactive closed sets is
    irrelevant... exponential"), or joining two already-persistent modules
    that share no synergy or complementarity edge (Table 2's case (F,F):
    "not synergetic, not complementary... a sequence of connected sets can
    be irrelevant"). There may be several distinct mechanisms that produce
    spurious organizations (free-species addition and latent join are the
    two implemented here) -- no attempt is made to enumerate or classify
    every such mechanism.

Why compute spurious organizations AT ALL, ever
-------------------------------------------------
Only for cross-validation against a full, pre-productive-novelty-theory
enumeration -- e.g. Centler et al. 2006's published organization tables,
computed by brute force before this decomposition existed. Reproducing
those tables exactly (species-for-species, not just by count) is a strong
end-to-end correctness check on the whole pipeline (ERC/closure computation
through the LP self-maintenance kernel). It is NOT what "compute
organizations efficiently" should mean, and it does not scale: the
spurious set can be exponentially larger than the fundamental one (see
_extend_with_free_species's docstring for a measured 2^15 case on a
72-species real network). Pass include_spurious=True only when you
specifically need the full naive enumeration for comparison; leave it off
(the default) for everything else, including all genome-scale work.

Pipeline (fundamental path, always run)
-----------------------------------------
  1. ERC computation                    erc.compute_ercs
  2. Fundamental relations              synergy.compute_synergies_basis_first,
                                         complementarity.compute_complementarities
  3. Fundamental hierarchy              hierarchy.build_hierarchy
  4. EPM/ESPM exploration               epm.compute_epms, epm.compute_espm
  5. LP self-maintenance verification   self_maintenance.check_self_maintenance

Steps 1-4 find every semi-organization reachable via fundamental synergy/
complementarity combination: a closed species set X with req(X) = 0. This
is necessary but not sufficient for the classical COT definition of self-
maintenance (Def 2.4 / self-maintainance in the paper): SSM says nothing
about whether a non-negative flux vector actually exists that keeps every
triggered reaction genuinely active while balancing net production against
consumption for every species. Step 5 promotes each SSM found into a
verified Organization by running the same LP check
(self_maintenance.check_self_maintenance) that
projects/Decomposition_Theorem/decomp/circuits.py already relies on when
decomposing organizations into E/F/fragile-circuit structure.

Bounding (not enumerating) the spurious tail
-----------------------------------------------
For each fundamental organization, it is cheap (one closure check per
candidate, no enumeration) to count how many "free" species -- never the
sole reactant of any reaction anywhere in the network -- could indepen-
dently be added to it. That count k gives 2^k - 1 as a lower bound on how
many additional spurious organizations exist in that fundamental
organization's free-species family alone (lower bound: distinct fundamental
organizations can generate overlapping families, and latent joins are a
separate, unbounded-here mechanism). This is computed unconditionally
(cheap) and reported per-organization and in aggregate, independent of
whether include_spurious actually enumerates them.

Public API
----------
compute_organizations(rn, *, network_id="", max_espm_order=10,
                       verify_organizations=True, include_spurious=False,
                       verbose=False, counters=None) -> OrganizationsResult
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
    can ever be "freely" (spuriously) addable to a semi-organization
    without a per-organization closure check ruling them out first.

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


def _free_species_bound(mask: int, rn_data, free_candidates: list[int]) -> tuple[int, int]:
    """
    Cheap, non-enumerating bound: for `mask` (a fundamental organization),
    count free candidates not already present whose SINGLE addition
    triggers no new reaction (one closure check each, no power-set).

    Returns (n_addable, spurious_lower_bound) where spurious_lower_bound =
    2**n_addable - 1 (every non-empty subset of the individually-safe
    candidates is at least a structural candidate; not all are guaranteed
    to independently pass a joint closure check, so this is a lower bound
    on the true count within this organization's free-species family, not
    an exact count -- see module docstring).
    """
    if not free_candidates:
        return 0, 0
    inv_idx = list(rn_data.species_to_reactions)
    n_addable = 0
    for c in free_candidates:
        if (mask >> c) & 1:
            continue
        candidate = mask | (1 << c)
        closed = closure_opt(rn_data.supp_q, rn_data.prod_q, candidate, inv_idx)
        if closed == candidate:
            n_addable += 1
    return n_addable, (2 ** n_addable - 1)


def _extend_with_free_species(so_masks: set[int], rn_data, free_candidates: list[int],
                               *, exhaustive_cap: int = 10) -> set[int]:
    """
    SPURIOUS-organization generator (only called when include_spurious=True).

    For every semi-organization mask in so_masks, try adding subsets of
    free_candidates not already present; keep those that don't trigger any
    new reaction (structural closure check only -- LP verification happens
    later, for every mask this returns). This is deliberately the
    generatively-irrelevant construction the paper's Sec. 3.2 describes
    (see module docstring) -- it exists only to reach a full naive
    enumeration for cross-validation, never as part of the recommended path.

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
        does double the count) but exactly the exponential blow-up the
        paper's whole methodology exists to make unnecessary. This
        function is now opt-in precisely because of that measurement.

    This is a genuine completeness trade-off above the cap: a "mixed"
    subset where candidate A is only safe to add together with B (neither
    alone, only the pair) would be missed if A or B individually fails the
    singleton test.
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


def _is_connected_raw(rn_data, X: int) -> bool:
    """
    Connectivity check per the paper's actual Def. conn_set: species
    s1, s2 directly connected iff SOME REACTION r (in the real, un-
    quotiented network) has {s1, s2} subset-or-equal supp(r) union prod(r).
    X is connected iff every pair of its species is connected (via a
    chain) within X.

    Deliberately uses rn_data.supp_raw / prod_raw (the real reaction data),
    NOT supp_q / prod_q (E0-quotiented). Using the quotiented version was a
    real bug in an earlier version of this module: E0 species never appear
    in supp_q/prod_q by construction (that is what quotienting means), so
    a connectivity check built on supp_q/prod_q structurally could never
    see an E0 species as connected to anything, regardless of the real
    network's topology, and an earlier version of this code responded by
    bypassing the check entirely rather than fixing which data it read.
    The paper's own Def. conn_set is stated over the real reaction network,
    so that is what this checks.

    Not currently used to gate anything in this module (see
    _saturate_latent_joins's docstring for why) -- kept as a correct,
    available utility, e.g. for reporting which spurious organizations are
    disconnected (Sos-external for that reason specifically) rather than
    filtering them out.
    """
    sp_list = [s for s in range(rn_data.n_species) if (X >> s) & 1]
    if len(sp_list) <= 1:
        return True

    adj: dict[int, set[int]] = {s: set() for s in sp_list}
    for r in range(rn_data.n_reactions):
        sq = rn_data.supp_raw[r]
        if not sq or (sq & X) != sq:
            continue  # reaction not active within X
        touched = sq | rn_data.prod_raw[r]
        rxn_sp = [s for s in sp_list if (touched >> s) & 1]
        for a_i, a in enumerate(rxn_sp):
            for b in rxn_sp[a_i + 1:]:
                adj[a].add(b)
                adj[b].add(a)

    start = sp_list[0]
    visited = {start}
    queue = [start]
    while queue:
        curr = queue.pop()
        for nb in adj[curr]:
            if nb not in visited:
                visited.add(nb)
                queue.append(nb)
    return len(visited) == len(sp_list)


def _saturate_latent_joins(so_masks: set[int], rn_data, *, max_rounds: int = 20,
                            pairwise_cap: int = 500) -> set[int]:
    """
    SPURIOUS-organization generator (only called when include_spurious=True).

    Closes the found semi-organization set under "latent join": epm.py's
    compute_epms/compute_espm deliberately never search for unions of
    already-SSM sets that share no fundamental synergy/complementarity edge
    -- Table 2's case (F,F) in the paper: "not synergetic, not
    complementary... a sequence of connected sets can be irrelevant."
    Provided only so a full naive enumeration can be recovered on demand
    for cross-validation (e.g. against Centler et al. 2006's published
    tables, computed before this productive-novelty theory existed).

    Confirmed necessary for that cross-validation, not theoretical: on
    Centler et al. 2006's "all_sugars" scenario, the published top
    organization (the entire 92-species network) is exactly such a latent
    join -- the DFS finds several smaller SOs that jointly cover the whole
    species set once unioned, but no single fundamental synergy/
    complementarity edge connects the two largest of them.

    Deliberately does NOT gate on _is_connected_raw. Two rounds of getting
    this wrong, worth recording precisely: a first version checked
    connectivity on the E0-quotiented data (supp_q/prod_q), which
    structurally can never see an E0 species as connected to anything
    (quotienting is what strips them out), and silently rejected the
    all_sugars=92 case; that version was fixed to check the REAL reaction
    data (supp_raw/prod_raw) instead, which is correct for what the paper
    calls connectivity (Def. conn_set) -- but then found ADP and AMP
    genuinely, correctly disconnected in Centler's own model (each has
    only an unconditional inflow and an independent decay reaction, never
    appearing in a shared reaction with any other species), which made
    the whole 92-species set fail the gate -- yet Centler's own 2006 paper
    publishes that exact 92-species set as an organization. The
    resolution: connectivity (Sos = CONNECTED semi-organizations, Lemma
    Sos_subset) is THIS PAPER's refinement of what counts as a
    *fundamental*/relevant organization, not a universal validity
    requirement on closed+SSM sets in the classical (Dittrich/Centler)
    sense this function is reconstructing. Fundamental organizations get
    connectivity for free, structurally, from the DFS's own combination
    rules (see epm.py: "Connectivity is NOT checked: it is structurally
    guaranteed by both combination rules") -- nothing here needs to
    re-check it. Spurious organizations are, by definition, outside that
    guarantee, and reconstructing Centler's pre-refinement published
    results specifically requires NOT imposing it.

    Two probes per round:
      (a) grand-union probe (cheap, O(k)): closure of the union of every
          mask currently known.
      (b) pairwise probe (O(k^2)): every pair of currently-known masks.
          Skipped once the known-SO count exceeds pairwise_cap.

    Runs to a fixed point (no new mask found) or max_rounds, whichever
    comes first.
    """
    from .closure import build_inv_idx, closure_opt, is_ssm

    inv_idx = build_inv_idx(rn_data.supp_q, rn_data.n_species)

    def _probe(candidate: int) -> int | None:
        closed = closure_opt(rn_data.supp_q, rn_data.prod_q, candidate, inv_idx)
        if not is_ssm(rn_data.supp_q, rn_data.prod_q, closed):
            return None
        return closed

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
    is_fundamental : bool            — True iff reached via the fundamental-
                     generator DFS alone (epm.compute_epms/compute_espm);
                     False iff only reachable via free-species-addition or
                     latent-join (see module docstring) -- always True
                     unless include_spurious=True was passed.
    is_organization : bool           — True iff LP self-maintenance verified
    flux          : dict[str, float] | None — reaction_name -> flux, if verified
    free_species_addable : int       — (fundamental orgs only) count of free
                     candidates whose single addition triggers no new
                     reaction; see _free_species_bound.
    spurious_lower_bound : int       — 2**free_species_addable - 1; a lower
                     bound on additional spurious organizations in this
                     organization's free-species family alone.
    """
    species_mask: int
    species_names: tuple
    order: int
    is_fundamental: bool
    is_organization: bool
    flux: dict | None = None
    free_species_addable: int = 0
    spurious_lower_bound: int = 0


@dataclass
class OrganizationsResult:
    """
    Full output of compute_organizations.

    semiorganizations : list[SemiOrganization] — every SSM found. Fundamental
                         only, unless include_spurious=True was passed.
    organizations     : list[SemiOrganization] — the subset that is LP-verified
                         self-maintaining (empty if verify_organizations=False)
    rn_data, ercs, hier, syn_result, comp_result, epm_result, espm_result
                      : intermediate pipeline artifacts, kept for reuse
                        (e.g. by the ERC-hierarchy visualization or by
                        projects/Decomposition_Theorem's decomposition step)
    stats             : dict — pipeline timing/counts. Notably
                        stats['spurious_lower_bound'] is always computed
                        (cheap) regardless of include_spurious -- see
                        module docstring's "Bounding" section.
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
    include_spurious: bool = False,
    verbose: bool = False,
    counters=None,
) -> OrganizationsResult:
    """
    Run the fundamental ERC -> fundamental relations -> hierarchy ->
    EPM/ESPM -> LP-verified-organizations pipeline on a pyCOT
    ReactionNetwork.

    By default this computes ONLY fundamental organizations -- those
    reachable by combining ERCs via fundamental synergy/complementarity
    alone (Theorem thm:fund_gen). This is the recommended, efficient path
    for all networks, including genome-scale ones.

    Parameters
    ----------
    rn : pyCOT ReactionNetwork (e.g. from pyCOT.io.functions.read_txt)
    network_id : label carried through metrics/instrumentation
    max_espm_order : safety cap on ESPM BFS depth (epm.compute_espm)
    verify_organizations : if False, skip the LP check entirely and
        return every SSM with is_organization=False -- useful when you only
        need the fast semi-organization lattice (e.g. for
        Decomposition_Theorem's own downstream per-circuit LP checks, which
        do their own verification at a finer grain).
    include_spurious : if True, additionally enumerate spurious
        organizations (free-species addition + latent join -- see module
        docstring) to reach a full naive enumeration. NOT the recommended
        default: this is exponential in the number of "free" species and
        exists only for cross-validation against pre-productive-novelty-
        theory published results (e.g. Centler et al. 2006). Leave False
        for genome-scale work.
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
    # back for display. Stripping here makes every mask this function
    # touches from this point on E0-free and canonical; E0 is added back
    # exactly once, at the point species names are produced.
    e0 = rn_data.E0_mask
    fundamental_masks = {m & ~e0 for m in espm_result.all_so_masks}

    so_order = {}
    for sp in epm_result.all_epm_masks:
        so_order[sp & ~e0] = 0
    for k, masks in espm_result.espm_by_order.items():
        for sp in masks:
            so_order[sp & ~e0] = k

    # ── Cheap free-species bound (always computed, never enumerated) ──────
    free_candidates = _free_candidate_species(rn_data)
    stats['n_free_candidate_species'] = len(free_candidates)
    bound_by_mask: dict[int, tuple[int, int]] = {
        m: _free_species_bound(m, rn_data, free_candidates) for m in fundamental_masks
    }
    stats['spurious_lower_bound'] = sum(b for _, b in bound_by_mask.values())

    all_masks = set(fundamental_masks)

    if include_spurious:
        # ── Stage 4.5: free-species extension ──────────────────────────
        all_masks = _extend_with_free_species(fundamental_masks, rn_data, free_candidates)
        n_extended = len(all_masks) - len(fundamental_masks)
        stats['n_free_species_extensions'] = n_extended
        stats['n_free_candidates_over_cap'] = len(free_candidates) > 10
        if verbose and n_extended:
            print(f"[organizations] Stage 4.5: free-species extension added "
                  f"{n_extended} spurious candidate sets ({len(free_candidates)} "
                  f"never-solely-consumed species considered).")

        # ── Stage 4.6: latent-join saturation ──────────────────────────
        n_before_joins = len(all_masks)
        all_masks = _saturate_latent_joins(all_masks, rn_data)
        n_joined = len(all_masks) - n_before_joins
        stats['n_latent_join_extensions'] = n_joined
        if verbose and n_joined:
            print(f"[organizations] Stage 4.6: latent-join saturation added "
                  f"{n_joined} spurious candidate sets not reachable via "
                  f"synergy/complementarity DFS alone.")

        # Assign an order to every newly-added mask: 1 + max(order of any
        # proper sub-mask already known), matching epm.py's own rule.
        for m in sorted(all_masks - set(so_order), key=lambda x: bin(x).count('1')):
            max_sub = -1
            for sub, sub_ord in so_order.items():
                if sub != m and (sub & m) == sub and sub_ord > max_sub:
                    max_sub = sub_ord
            so_order[m] = max_sub + 1

    # ── LP self-maintenance verification (fundamental + spurious alike) ───
    if verbose:
        print(f"[organizations] Stage 5: LP-verifying "
              f"{len(all_masks)} semi-organizations "
              f"({len(fundamental_masks)} fundamental"
              f"{f', {len(all_masks) - len(fundamental_masks)} spurious' if include_spurious else ''})...")

    semiorganizations: list[SemiOrganization] = []
    organizations: list[SemiOrganization] = []
    n_verified_true = 0
    n_fund_verified_true = 0
    for sp_mask in sorted(all_masks):
        order = so_order.get(sp_mask, 0)
        is_fund = sp_mask in fundamental_masks
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
                if is_fund:
                    n_fund_verified_true += 1
                sub_rn = rn.sub_reaction_network(species_objs)
                rxn_names = [r.node.name for r in
                             sorted(sub_rn.reactions(), key=lambda r: r.node.index)]
                flux_dict = {name: float(v) for name, v in zip(rxn_names, flux_vec)}
        n_addable, spurious_bound = bound_by_mask.get(sp_mask, (0, 0))
        so = SemiOrganization(
            species_mask=sp_mask, species_names=full_names, order=order,
            is_fundamental=is_fund, is_organization=is_org, flux=flux_dict,
            free_species_addable=n_addable, spurious_lower_bound=spurious_bound,
        )
        semiorganizations.append(so)
        if is_org:
            organizations.append(so)

    stats['n_semiorganizations'] = len(semiorganizations)
    stats['n_organizations'] = n_verified_true
    stats['n_fundamental_organizations'] = n_fund_verified_true
    if include_spurious:
        stats['n_spurious_organizations'] = n_verified_true - n_fund_verified_true

    if verbose and verify_organizations:
        print(f"[organizations] Done: {len(semiorganizations)} semi-organizations, "
              f"{n_verified_true} verified as true organizations "
              f"({n_fund_verified_true} fundamental).")

    return OrganizationsResult(
        semiorganizations=semiorganizations,
        organizations=organizations,
        rn_data=rn_data, ercs=ercs, hier=hier,
        syn_result=syn_result, comp_result=comp_result,
        epm_result=epm_result, espm_result=espm_result,
        stats=stats,
    )
