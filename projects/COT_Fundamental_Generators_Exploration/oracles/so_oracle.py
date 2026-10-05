"""
oracles/so_oracle.py — Brute-force SO0/SOi oracle (paper Defs 30-31).

Theory
------
A semi-organization (SO) is a species set X that is:
  - Closed      : prod(R_X) ⊆ X  (no new species are produced outside X)
  - SSM         : req(X) = ∅     (X produces everything its reactions need)
                  i.e. supp(R_X) ⊆ prod(R_X)
  - Reactive    : R_X ≠ ∅        (at least one reaction in X)
  - Connected   : X cannot be partitioned into two disconnected reactive subsets

Every reactive closed set is a union of ERCs, so we enumerate ERC-subsets.

Elementary SO (SO0, Def 30): a SO X such that the only SO strictly contained
  in X is E_∅ (the inflow-only set, represented as the empty set in the
  quotiented RN). In the quotiented network: SO0 = minimal non-empty SO.

Higher-order SO (SOi, i >= 1, Def 31): a SO that is not elementary (it
  contains at least one elementary SO). Generative index i: the max
  generative index among its proper sub-SOs is i-1.

Public API
----------
elementary_so_oracle(rn, ercs) -> list[int]
    Returns species bitmasks of all elementary SOs.  Only feasible for
    |ercs| <= ~20.

so_hierarchy_oracle(rn, ercs) -> dict  with keys "elementary_masks",
    "so_by_order", "all_sos", "so_order"
    Full SO structure including elementary and all-order higher-order SOs.
"""
from __future__ import annotations

from itertools import combinations


# ---------------------------------------------------------------------------
# Reaction-network helpers (mirror of pyCOT.analysis.organizations.so_search
# for oracle independence)
# ---------------------------------------------------------------------------

def _closure(supp_q, prod_q, n_rxn: int, seed: int) -> int:
    """Horn fixpoint closure of a species bitmask."""
    X = seed
    changed = True
    while changed:
        changed = False
        for r in range(n_rxn):
            sq = supp_q[r]
            if sq and (sq & X) == sq:
                new = prod_q[r] & ~X
                if new:
                    X |= new
                    changed = True
    return X


def _is_ssm(supp_q, prod_q, n_rxn: int, X: int) -> bool:
    """True if X is semi-self-maintaining: supp(R_X) ⊆ prod(R_X)."""
    supp_tot = prod_tot = 0
    for r in range(n_rxn):
        sq = supp_q[r]
        if sq and (sq & X) == sq:
            supp_tot |= sq
            prod_tot |= prod_q[r]
    return (supp_tot & ~prod_tot) == 0


def _is_reactive(supp_q, n_rxn: int, X: int) -> bool:
    return any(supp_q[r] and (supp_q[r] & X) == supp_q[r] for r in range(n_rxn))


def _is_connected(supp_q, prod_q, n_rxn: int, n_species: int, X: int) -> bool:
    """
    True if X is connected as a reaction hypergraph.
    Two species are in the same component if they co-appear in any active reaction
    (either in its support or in its products).
    """
    sp_list = [s for s in range(n_species) if (X >> s) & 1]
    if len(sp_list) <= 1:
        return True

    adj: dict[int, set[int]] = {s: set() for s in sp_list}
    sp_set = set(sp_list)

    for r in range(n_rxn):
        sq = supp_q[r]
        if not sq or (sq & X) != sq:
            continue
        pq = prod_q[r]
        # species in X that participate in this reaction
        rxn_sp = [s for s in sp_list if ((sq >> s) & 1) or ((pq >> s) & 1)]
        for idx, a in enumerate(rxn_sp):
            for b in rxn_sp[idx + 1:]:
                adj[a].add(b)
                adj[b].add(a)

    # BFS
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


# ---------------------------------------------------------------------------
# Oracle entry points
# ---------------------------------------------------------------------------

def _compute_all_sos(rn, ercs) -> list[int]:
    """
    Enumerate all SOs by iterating over ERC subsets.

    Returns sorted list of species bitmasks for all SOs.
    Raises RuntimeError if |ercs| > 25 (infeasibly many subsets).
    """
    n = len(ercs)
    if n > 25:
        raise RuntimeError(
            f"SO oracle: {n} ERCs → {2**n:,} subsets to check. "
            "Too many for brute force (limit: 25 ERCs). "
            "Use compare_oracles.py with MIN_REACTIONS / MAX_REACTIONS to select smaller networks."
        )

    masks  = [e.species_mask for e in ercs]
    sq     = rn.supp_q
    pq     = rn.prod_q
    n_rxn  = rn.n_reactions
    n_sp   = rn.n_species

    so_set: set[int] = set()

    for bits in range(1, 1 << n):
        # Union of ERC species masks for this subset
        seed = 0
        for i in range(n):
            if (bits >> i) & 1:
                seed |= masks[i]

        # Closure (handles synergetic reactions that activate new ERCs)
        X = _closure(sq, pq, n_rxn, seed)

        if not _is_reactive(sq, n_rxn, X):
            continue
        if not _is_ssm(sq, pq, n_rxn, X):
            continue
        if not _is_connected(sq, pq, n_rxn, n_sp, X):
            continue

        so_set.add(X)

    return sorted(so_set)


def elementary_so_oracle(rn, ercs) -> list[int]:
    """
    Brute-force elementary-SO (SO0) oracle.

    An elementary SO is a minimal non-empty SO in the quotiented network
    (E_∅ = ∅ since inflow species are stripped).

    Returns sorted list of species bitmasks for all elementary SOs.
    """
    all_sos = _compute_all_sos(rn, ercs)

    elementary = []
    for X in all_sos:
        # X is elementary if no proper non-empty SO is strictly contained in X
        is_elementary = not any(
            Y != X and Y != 0 and (Y & X) == Y
            for Y in all_sos
        )
        if is_elementary:
            elementary.append(X)

    return elementary


def so_hierarchy_oracle(rn, ercs) -> dict:
    """
    Full SO structure: elementary SOs (SO0) and higher-order SOs (SOi) by order.

    Returns
    -------
    dict with:
      "all_sos"          : list[int]           — all SO bitmasks
      "elementary_masks" : list[int]           — elementary (SO0) bitmasks
      "so_by_order"       : dict[int, list[int]] — {order: [SO bitmasks]}, order >= 1
      "so_order"          : dict[int, int]      — {SO bitmask: order}
    """
    all_sos = _compute_all_sos(rn, ercs)

    # Compute generative index for each SO.
    # Order 0 = elementary SO (minimal non-empty SO)
    # Order k = SO whose maximum proper sub-SO order is k-1
    so_order: dict[int, int] = {}
    elementary: list[int] = []

    # Process SOs from smallest to largest
    sorted_sos = sorted(all_sos, key=lambda x: bin(x).count('1'))

    for X in sorted_sos:
        # Find proper sub-SOs
        sub_sos = [Y for Y in all_sos if Y != X and Y != 0 and (Y & X) == Y]
        if not sub_sos:
            # No proper sub-SO → elementary (order 0)
            so_order[X] = 0
            elementary.append(X)
        else:
            max_sub_order = max(so_order.get(Y, 0) for Y in sub_sos)
            so_order[X] = max_sub_order + 1

    so_by_order: dict[int, list[int]] = {}
    for X, ord_k in so_order.items():
        if ord_k > 0:   # higher-order SOs have order >= 1
            so_by_order.setdefault(ord_k, []).append(X)

    return {
        "all_sos":          all_sos,
        "elementary_masks": elementary,
        "so_by_order":      so_by_order,
        "so_order":         so_order,
    }


def elementary_so_oracle_set(rn, ercs) -> set[int]:
    """Return elementary-SO bitmasks as a set (for comparison with efficient algorithm)."""
    return set(elementary_so_oracle(rn, ercs))
