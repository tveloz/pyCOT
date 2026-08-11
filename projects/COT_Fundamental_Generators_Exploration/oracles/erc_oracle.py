"""
oracles/erc_oracle.py — Brute-force ERC discovery (correctness anchor).

erc_oracle(supp, prod) -> dict[int, list[int]]
    For each reaction r with supp[r] ≠ 0, compute
        erc_mask = closure_oracle(supp, prod, supp[r])
    Group reactions by their erc_mask.
    Returns dict: erc_mask → list[reaction_indices].

minbas_oracle(supp_list) -> list[int]
    Return ⊆-minimal elements of supp_list.

req_prod_oracle(supp, prod, mask) -> (int, int)
    Compute (req_mask, prod_mask) for a closed set mask.

Every function here is deliberately simple and obviously correct.
The optimized equivalents in cot_gen/ must return identical results.
"""
from __future__ import annotations

from .closure_oracle import closure_oracle


def erc_oracle(
    supp: list[int] | tuple[int, ...],
    prod: list[int] | tuple[int, ...],
) -> dict[int, list[int]]:
    """
    Compute ERC groups by brute-force closure.

    Parameters
    ----------
    supp : support bitmasks (quotiented — E0 stripped).
    prod : product bitmasks (quotiented).

    Returns
    -------
    dict mapping erc_mask → sorted list of reaction indices.
    """
    groups: dict[int, list[int]] = {}
    for r_idx, s in enumerate(supp):
        if s == 0:
            continue
        erc_mask = closure_oracle(supp, prod, s)
        if erc_mask not in groups:
            groups[erc_mask] = []
        groups[erc_mask].append(r_idx)
    return groups


def minbas_oracle(masks: list[int]) -> list[int]:
    """
    Return the ⊆-antichain of masks: remove any mask that is a strict
    superset of another.

    Duplicates are removed before the antichain computation.
    """
    unique = list(dict.fromkeys(masks))
    minimal = []
    for i, m_i in enumerate(unique):
        dominated = any(
            m_j != m_i and (m_j & m_i) == m_j
            for j, m_j in enumerate(unique)
            if j != i
        )
        if not dominated:
            minimal.append(m_i)
    return minimal


def req_prod_oracle(
    supp: list[int] | tuple[int, ...],
    prod: list[int] | tuple[int, ...],
    mask: int,
) -> tuple[int, int]:
    """
    Compute (req_mask, prod_mask) for the closed set `mask`.

    req(mask)  = ( ⋃ supp[r] : supp[r] ⊆ mask ) \\ ( ⋃ prod[r] : supp[r] ⊆ mask )
    prod(mask) =   ⋃ prod[r] : supp[r] ⊆ mask
    """
    agg_supp = 0
    agg_prod = 0
    for s, p in zip(supp, prod):
        if (s & ~mask) == 0:  # supp[r] ⊆ mask
            agg_supp |= s
            agg_prod |= p
    return agg_supp & ~agg_prod, agg_prod


def compute_ercs_oracle(
    supp: list[int] | tuple[int, ...],
    prod: list[int] | tuple[int, ...],
) -> list[dict]:
    """
    Full brute-force ERC computation (oracle interface).

    Returns a list of dicts (one per ERC), each with keys:
        species_mask, reaction_indices, min_bases, req_mask, prod_mask

    List is sorted by species_mask for determinism.
    """
    groups = erc_oracle(supp, prod)
    result = []
    for erc_mask, r_indices in sorted(groups.items()):
        supports = [supp[r] for r in r_indices]
        min_bases = minbas_oracle(supports)
        req, prod_m = req_prod_oracle(supp, prod, erc_mask)
        result.append({
            "species_mask": erc_mask,
            "reaction_indices": sorted(r_indices),
            "min_bases": min_bases,
            "req_mask": req,
            "prod_mask": prod_m,
        })
    return result
