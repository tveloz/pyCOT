"""
oracles/complementarity_oracle.py — Brute-force complementarity oracle (paper Defs 22–26).

comp_basic_set(ercs) -> set[tuple[int,int,int,int,bool]]
    All basic complementary pairs as (i, j, fwd_supply, bwd_supply, chain) with i < j.
    fwd_supply = prod_i & req_j, bwd_supply = prod_j & req_i.
    chain = True iff i, j are hierarchy-comparable (intra-chain pair). Comparable
    pairs ARE included: the subset-to-superset direction is provably always 0,
    but the superset-to-subset direction can be nonempty (see complementarity.py
    module docstring for the worked example).

comp_fund_set(ercs) -> set[tuple[int,int,int]]
    Fundamental complementarities as (prod_idx, cons_idx, species_bit).
    prod_idx ∈ minprod(s), cons_idx ∈ mincons(s).

No hierarchy object needed — containment is derived directly from species_masks.
"""
from __future__ import annotations

from itertools import combinations


def comp_basic_set(ercs) -> set[tuple]:
    """
    Brute-force basic complementarity oracle.

    (E_i, E_j) are basic complementary if prod_i & req_j ≠ 0 or prod_j & req_i ≠ 0.
    Comparable pairs are included (see module docstring): only the subset-to-
    superset direction is provably always 0, not the reverse.
    Returns set of (i, j, fwd_supply, bwd_supply, chain) with i < j.
    """
    masks = [e.species_mask for e in ercs]
    reqs  = [e.req_mask     for e in ercs]
    prods = [e.prod_mask    for e in ercs]
    n     = len(ercs)

    result: set[tuple] = set()
    for i, j in combinations(range(n), 2):
        mi, mj = masks[i], masks[j]
        chain = (mi & mj) == mi or (mi & mj) == mj
        fwd = prods[i] & reqs[j]   # E_i supplies to E_j
        bwd = prods[j] & reqs[i]   # E_j supplies to E_i
        if fwd or bwd:
            result.add((i, j, fwd, bwd, chain))

    return result


def comp_fund_set(ercs) -> set[tuple]:
    """
    Brute-force fundamental complementarity oracle.

    For each species s:
      minprod(s) = ⊆-minimal ERCs in {E : s ∈ prod(R_E)}
      mincons(s) = ⊆-minimal ERCs in {E : s ∈ req(E)}

    Returns {(prod_idx, cons_idx, species_bit)} for all valid combinations.
    """
    masks = [e.species_mask for e in ercs]
    reqs  = [e.req_mask     for e in ercs]
    prods = [e.prod_mask    for e in ercs]
    n     = len(ercs)

    all_masks = masks + reqs + prods
    max_bit = max((m.bit_length() for m in all_masks), default=0)

    result: set[tuple] = set()

    for s_bit in range(max_bit):
        s_mask = 1 << s_bit

        producers = [i for i in range(n) if prods[i] & s_mask]
        consumers = [i for i in range(n) if reqs[i] & s_mask]

        if not producers or not consumers:
            continue

        # minprod(s): no strictly contained ERC also produces s
        min_prods = [
            i for i in producers
            if not any(
                j != i
                and (masks[j] & masks[i]) == masks[j]
                and masks[j] != masks[i]
                for j in producers
            )
        ]

        # mincons(s): no strictly contained ERC also requires s
        min_cons = [
            i for i in consumers
            if not any(
                j != i
                and (masks[j] & masks[i]) == masks[j]
                and masks[j] != masks[i]
                for j in consumers
            )
        ]

        for pi in min_prods:
            for ci in min_cons:
                if pi != ci:
                    result.add((pi, ci, s_bit))

    return result
