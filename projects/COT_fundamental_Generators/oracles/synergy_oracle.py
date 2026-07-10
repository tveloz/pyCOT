"""
oracles/synergy_oracle.py — Brute-force synergy detection (correctness anchor).

synergy_oracle(ercs) -> list[dict]
    For each ordered pair (i, j, k) compute all basic synergies by definition:
        ∃ b ∈ min_bases[k]: b ⊆ (mask_i | mask_j) ∧ b ⊄ mask_i ∧ b ⊄ mask_j
    No hierarchy needed — just containment checks from masks.
    Returns list of {"i", "j", "k"} dicts.

This is O(|E|³ × max|MinBas|) but uses no caching — obviously correct.
"""
from __future__ import annotations

from itertools import combinations


def synergy_oracle(ercs) -> list[dict]:
    """
    Brute-force basic synergy detection.

    Parameters
    ----------
    ercs : list[ERCData]  (species_mask and min_bases must be set)

    Returns
    -------
    list[dict] with keys "i", "j", "k", sorted for determinism.
    """
    masks     = [e.species_mask for e in ercs]
    min_bases = [e.min_bases    for e in ercs]
    n         = len(ercs)

    result = []
    for i, j in combinations(range(n), 2):
        mi, mj = masks[i], masks[j]
        # Skip comparable pairs
        if (mi & mj) == mi or (mi & mj) == mj:
            continue
        joint = mi | mj
        for k in range(n):
            if k == i or k == j:
                continue
            mk = masks[k]
            # Skip if target contained by either base
            if (mk & mi) == mk or (mk & mj) == mk:
                continue
            # Check joint cover condition
            for b in min_bases[k]:
                if (b & joint) == b and (b & mi) != b and (b & mj) != b:
                    result.append({"i": i, "j": j, "k": k})
                    break

    return sorted(result, key=lambda x: (x["i"], x["j"], x["k"]))


def synergy_set(ercs) -> set[tuple[int, int, int]]:
    """Return basic synergies as a set of (i, j, k) triples."""
    return {(r["i"], r["j"], r["k"]) for r in synergy_oracle(ercs)}


def fundamental_synergy_set(ercs) -> set[tuple[int, int, int]]:
    """
    Brute-force fundamental synergy oracle (paper Def 21).

    Steps
    -----
    1. Find all basic synergies.
    2. Filter to maximal synergies: for each pair (i,j), keep only targets k
       not strictly contained in any other target k' for the same pair.
    3. Filter maximal to fundamental: (i,j,k) is fundamental if
         - no i' ⊊ i (strict subset) such that {i',j} → k is a maximal synergy
         - no j' ⊊ j (strict subset) such that {i,j'} → k is a maximal synergy
       where {a,b} denotes the unordered pair (stored as (min,max,k) in the set).
    """
    basic = synergy_set(ercs)
    masks = [e.species_mask for e in ercs]

    # Step 1→2: maximal synergies per pair (i,j)
    pair_targets: dict[tuple[int, int], list[int]] = {}
    for (i, j, k) in basic:
        pair_targets.setdefault((i, j), []).append(k)

    maximal: set[tuple[int, int, int]] = set()
    for (i, j), targets in pair_targets.items():
        for k in targets:
            mk = masks[k]
            # k is maximal if no other target k' strictly contains it
            if not any(
                k2 != k and (mk & masks[k2]) == mk   # mk ⊆ masks[k2]
                for k2 in targets
            ):
                maximal.add((i, j, k))

    # Step 2→3: fundamental = maximal with no reducible slot
    # Lookup: unordered pair key → set of maximal k's
    pair_key_to_k: dict[frozenset, set[int]] = {}
    for (i, j, k) in maximal:
        pair_key_to_k.setdefault(frozenset([i, j]), set()).add(k)

    fundamental: set[tuple[int, int, int]] = set()
    for (i, j, k) in maximal:
        mi, mj = masks[i], masks[j]

        # Check i-slot: is there i' ⊊ i with {i',j} →^max k?
        i_reducible = any(
            i2 != i
            and (masks[i2] & mi) == masks[i2]   # masks[i2] ⊆ masks[i]
            and masks[i2] != mi                  # strict
            and k in pair_key_to_k.get(frozenset([i2, j]), set())
            for i2 in range(len(ercs))
        )
        if i_reducible:
            continue

        # Check j-slot: is there j' ⊊ j with {i,j'} →^max k?
        j_reducible = any(
            j2 != j
            and (masks[j2] & mj) == masks[j2]
            and masks[j2] != mj
            and k in pair_key_to_k.get(frozenset([i, j2]), set())
            for j2 in range(len(ercs))
        )
        if j_reducible:
            continue

        fundamental.add((i, j, k))

    return fundamental
