"""
self_maintenance.py — LP verification of true (Def 2.4) self-maintenance.

This is the validated LP kernel originally developed in
Persistent_Modules_Generator.py (moved here verbatim, logic and comments
unchanged) and already relied on directly by
projects/Decomposition_Theorem/decomp/circuits.py to decide, per fragile
circuit, whether an organization-candidate genuinely self-maintains.

Why this module exists separately from SO0/SOi search (so_search.py)
-----------------------------------------------------------------------
so_search.py's "SSM" (semi-self-maintaining, req_mask == 0) is a *necessary but
not sufficient* condition for the classical COT notion of self-maintenance.
SSM only says every species some active reaction needs is produced by
*some* active reaction in the set -- it says nothing about whether a
non-negative flux vector actually exists that keeps every triggered
reaction genuinely active while balancing production against consumption.
That stronger, LP-verified condition is what this module checks; see
organizations.py for how the two stages compose into the actual
"compute organizations" pipeline.

Public API
----------
minimize_sv(S, epsilon=1e-6, method='highs', margin=1e-9) -> [bool, flux|None]
check_self_maintenance(species_set, RN, epsilon=1e-6) -> (bool, flux, production)
diagnose_self_maintenance(species_set, RN, epsilon=1e-6) -> dict
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linprog


def minimize_sv(S, epsilon=1e-6, method='highs', margin=1e-9):
    """
    Check self-maintenance: is there v ≥ epsilon (componentwise), sum(v) = n_reactions,
    with S·v ≥ -margin?

    Mathematical condition (COT, Def 2.4 / flux cone V(X)): X is self-maintaining
    iff there exists v with v_r > 0 STRICTLY for every r in R_X (every reaction
    triggered by X — i.e. every reaction whose reactants lie in X — and S·v ≥ 0.
    A reaction in R_X is not optional: if its reactants are present, it is part
    of the sub-network and must be assigned a genuinely positive rate, not "off".

    Key design choices:
      - v ≥ epsilon, NOT v ≥ 0: every reaction in R_X gets a strictly positive
        floor. Letting a triggered reaction sit at v_r = 0 (an earlier version
        of this function did, via bounds=(0, None)) lets the LP silently
        "switch off" reactions whose reactants are present but whose own
        consumption would break feasibility — e.g. a pure decay/outflow
        reaction `s -> ∅` for a species s that is otherwise only regenerated
        catalytically (net zero) elsewhere: with v_r allowed to be 0, the
        solver happily sets the decay's flux to zero and never notices s is
        actually being drained, wrongly reporting self-maintenance. Caught via
        networks/FarmVariants/Farm.txt: 'infr' has a decay reaction (R16) whose
        only counter-source (R17) needs 'farmer'; without this floor, a species
        set containing infr but not farmer was (wrongly) reported as an
        elementary organization, with R16 sitting at flux 0 in the witness.
      - epsilon default (1e-6) is deliberately >> margin default (1e-9): if the
        two are close, a forced-positive drain of size ~epsilon can be absorbed
        by the -margin slack on the Sv >= -margin constraint, silently
        reintroducing the same bug at the numerical level. Verified empirically
        on the Farm.txt case above: epsilon=1e-9 or 1e-7 (too close to
        margin=1e-9) still wrongly reports self-maintaining; epsilon=1e-6 or
        larger correctly reports infeasible.
      - sum(v) = n_reactions: normalization prevents comparing across
        different overall flux scales; combined with v >= epsilon this still
        leaves a nonempty box whenever n_reactions * epsilon <= n_reactions,
        i.e. always (epsilon <= 1), so it does not by itself create infeasibility.
      - S·v ≥ -margin: small negative tolerance (not S·v ≥ +margin, which
        wrongly rejects exactly-balanced net-zero exchange pairs, e.g. a
        catalyst appearing on both sides of a reaction with equal coefficient).
    """
    S_array = np.asarray(S)
    n_species, n_reactions = S_array.shape
    if n_reactions == 0:
        return [False, None]

    # Objective: minimise total flux (LP feasibility — any feasible point works)
    c = np.ones(n_reactions)

    # S·v ≥ -margin  ↔  -S·v ≤ +margin
    A_ub = -S_array
    b_ub = margin * np.ones(n_species)

    # Normalization: sum(v) = n_reactions  (prevents trivial v = 0)
    A_eq = np.ones((1, n_reactions))
    b_eq = [float(n_reactions)]

    # v ≥ epsilon — every reaction of R_X is triggered and MUST fire at a
    # genuinely positive rate (Def 2.4); it is not optional just because its
    # own consumption might otherwise break feasibility.
    bounds = [(epsilon, None) for _ in range(n_reactions)]

    try:
        result = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                         bounds=bounds, method=method)
        if result.success:
            # Trust the LP solver: HiGHS reports success iff the primal
            # constraints are satisfied within its internal feasibility
            # tolerance (~1e-7).  Re-checking with our tighter margin=1e-9
            # would spuriously reject boundary solutions due to floating-point
            # rounding when recomputing S @ v.
            return [True, result.x]
        else:
            return [False, None]
    except Exception as e:
        print(f"Linear programming error: {e}")
        return [False, None]


def check_self_maintenance(species_set, RN, epsilon=1e-6):
    """
    Check if a set of species is self-maintaining using linear programming.

    Parameters
    ----------
    species_set : list of Species
        Set of species to check
    RN : ReactionNetwork
        The reaction network
    epsilon : float
        Minimum flux for active (triggered) reactions — must stay >> minimize_sv's
        margin (1e-9) or the strict-positivity floor can be numerically absorbed
        by the Sv >= -margin slack; see minimize_sv's docstring.

    Returns
    -------
    tuple
        (is_self_maintaining: bool, flux_vector: np.array or None, production_vector: np.array or None)
    """
    try:
        if not species_set:
            return (False, None, None)

        # Create sub-reaction network with species and their activated reactions
        sub_RN = RN.sub_reaction_network(species_set)

        # Get stoichiometric matrix
        S = sub_RN.stoichiometry_matrix()

        # Check if there are any reactions
        if S.shape[1] == 0:
            return (False, None, None)

        # Solve for self-maintaining vector
        res = minimize_sv(S, epsilon=epsilon, method='highs')

        if res[0]:  # If successful
            flux_vector = res[1]
            production_vector = np.asarray(S) @ flux_vector
            return (True, flux_vector, production_vector)
        else:
            return (False, None, None)

    except Exception as e:
        print(f"Error in check_self_maintenance: {e}")
        return (False, None, None)


def diagnose_self_maintenance(species_set, RN, epsilon=1e-6):
    """
    Return a human-readable diagnosis of the self-maintenance status of a
    species set.

    Parameters
    ----------
    species_set : list of Species
        Species forming the semi-organisation to test.
    RN : ReactionNetwork
        The full reaction network.
    epsilon : float
        Minimum flux threshold (same as used by compute_all_organizations).

    Returns
    -------
    dict with keys:
        'is_org'    : bool
        'reactions' : {rxn_name: flux}    -- present only when is_org=True
        'production': {sp_name: net_prod} -- present only when is_org=True
        'blocking'  : {sp_name: slack}    -- present only when is_org=False;
                      slack > 0 means that species has a production deficit.
                      Value 'no_producer' means no reaction in the sub-network
                      can produce it at all.
    """
    if not species_set:
        return {'is_org': False, 'blocking': {}}

    sub_RN = RN.sub_reaction_network(species_set)
    S_obj  = sub_RN.stoichiometry_matrix()
    S      = np.asarray(S_obj, dtype=float)
    n_sp, n_rx = S.shape
    sp_names = list(S_obj.species)
    rx_names = list(S_obj.reactions)

    if n_rx == 0:
        return {'is_org': False,
                'blocking': {sp.name: 'no_reactions' for sp in species_set}}

    # ── Pass 1: standard SM check ──────────────────────────────────────────
    res = minimize_sv(S, epsilon=epsilon)
    if res[0]:
        flux   = res[1]
        prod   = S @ flux
        return {
            'is_org':     True,
            'reactions':  {rx_names[j]: float(flux[j]) for j in range(n_rx)},
            'production': {sp_names[i]: float(prod[i]) for i in range(n_sp)},
        }

    # ── Pass 2: structural no-producer check ──────────────────────────────
    # Exclude all-zero rows (pure catalysts in this sub-network): their LP
    # constraint is 0·v ≥ 0, trivially satisfied, so they are not blockers.
    blocking = {}
    for i, sp in enumerate(sp_names):
        if np.all(S[i, :] <= 0) and np.any(S[i, :] < 0):
            blocking[sp] = 'no_producer'

    # ── Pass 3: slack LP — find species with production deficit ───────────
    # min  sum(s)
    # s.t. S·v + s >= -margin  (i.e. -S·v - s <= margin)
    #      sum(v) = n_rx  (normalization — avoid trivial v = 0)
    #      v >= 0,  s >= 0
    # slack > 0 for species i means even the best v ≥ 0 leaves species i
    # with a production deficit that must be covered by the slack variable.
    margin_diag = 1e-9
    c    = np.concatenate([np.zeros(n_rx), np.ones(n_sp)])
    A_ub = np.hstack([-S, -np.eye(n_sp)])
    b_ub = margin_diag * np.ones(n_sp)
    A_eq_d = np.concatenate([np.ones(n_rx), np.zeros(n_sp)])[np.newaxis, :]
    b_eq_d = [float(n_rx)]
    bds  = [(0, None)] * n_rx + [(0, None)] * n_sp

    r2 = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq_d, b_eq=b_eq_d,
                 bounds=bds, method='highs')
    if r2.success:
        for i, sp in enumerate(sp_names):
            slack = float(r2.x[n_rx + i])
            if slack > 1e-7 and sp not in blocking:
                blocking[sp] = round(slack, 6)

    return {'is_org': False, 'blocking': blocking}
