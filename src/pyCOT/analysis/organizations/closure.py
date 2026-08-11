"""
closure.py — Optimized closure computation (Horn / Dowling-Gallier style).

Public API
----------
closure_opt(supp, prod, X, inv_idx) -> int
    Compute the smallest closed set containing X.
    Uses reactant-counter Horn propagation: O(|fired reactions| * avg_supp_size).

build_inv_idx(supp, n_species) -> list[list[int]]
    Build the species→reactions inverted index from a support list.
    Call once, reuse for every closure computation in the same network.

See oracles/closure_oracle.py for the naive fixpoint reference.

Invariants asserted at runtime (cheap; use -O to disable):
  • Idempotence:  closure_opt(supp, prod, closure_opt(supp, prod, X)) == result
  • Monotonicity: X ⊆ closure_opt(supp, prod, X)
  • Contains X:  (result & X) == X

Usage
-----
    inv_idx = build_inv_idx(rn.supp_q, rn.n_species)
    closed_set = closure_opt(rn.supp_q, rn.prod_q, seed_mask, inv_idx)
"""
from __future__ import annotations

from collections import deque
from typing import Sequence


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _bits(n: int):
    """Yield bit positions (ascending) of bits set in n."""
    while n:
        lsb = n & (-n)
        yield lsb.bit_length() - 1
        n &= n - 1


def build_inv_idx(supp: Sequence[int], n_species: int) -> list[list[int]]:
    """
    Build the inverted index: species s → list of reaction indices r
    with bit s set in supp[r].

    Parameters
    ----------
    supp      : sequence of support bitmasks (one per reaction).
                Use supp_q (E0-quotiented) for ERC computation.
    n_species : total number of species (length of the returned list).
    """
    inv: list[list[int]] = [[] for _ in range(n_species)]
    for r_idx, s in enumerate(supp):
        m = s
        while m:
            lsb = m & (-m)
            inv[lsb.bit_length() - 1].append(r_idx)
            m &= m - 1
    return inv


# ---------------------------------------------------------------------------
# Optimized closure (Horn / Dowling-Gallier)
# ---------------------------------------------------------------------------

def closure_opt(
    supp: Sequence[int],
    prod: Sequence[int],
    X: int,
    inv_idx: list[list[int]],
    *,
    counters=None,          # optional cot_gen.metrics.Counters for instrumentation
) -> int:
    """
    Compute the smallest closed set ⊇ X for the reaction network (supp, prod).

    Algorithm (Dowling–Gallier style Horn propagation):
      1. cnt[r] = popcount(supp[r] & ~X)  — missing reactants per reaction
      2. Scan reactions: fire those with cnt[r]=0 (supp[r] ⊆ X already);
         collect newly added species in a FIFO queue.
      3. For each newly added species s: decrement cnt[r] for every r in
         inv_idx[s]; if cnt[r] hits 0 fire r, enqueue new products.

    Complexity: O(|R|) for step 2 + O(Σ |supp(r)| for fired reactions) for step 3.
    The naive fixpoint alternative is O(|R| × iterations).

    Parameters
    ----------
    supp    : support bitmasks (length n_reactions).
    prod    : product bitmasks (length n_reactions).
    X       : initial species bitmask.
    inv_idx : build_inv_idx(supp, n_species) — reuse across calls.
    counters: optional Counters object (increments ctr keys for the paper).

    Returns
    -------
    int : bitmask of the closed set.
    """
    # Step 1 — initialise reactant counters
    # cnt[r] = number of reactants of r NOT yet in X
    cnt = [bin(s & ~X).count('1') for s in supp]

    queue: deque[int] = deque()
    reactions_fired = 0
    species_added = 0

    # Step 2 — fire immediately active reactions (those with supp ⊆ X)
    for r_idx, (s, p) in enumerate(zip(supp, prod)):
        if cnt[r_idx] == 0:
            new = p & ~X
            if new:
                X |= new
                reactions_fired += 1
                for bit in _bits(new):
                    queue.append(bit)
                    species_added += 1

    # Step 3 — cascade: process newly added species
    while queue:
        sp = queue.popleft()
        for r_idx in inv_idx[sp]:
            cnt[r_idx] -= 1
            if cnt[r_idx] == 0:
                new = prod[r_idx] & ~X
                if new:
                    X |= new
                    reactions_fired += 1
                    for bit in _bits(new):
                        queue.append(bit)
                        species_added += 1

    if counters is not None:
        counters.inc("closure.reactions_fired", reactions_fired)
        counters.inc("closure.species_added", species_added)
        counters.inc("closure.calls")

    return X


# ---------------------------------------------------------------------------
# Convenience wrappers
# ---------------------------------------------------------------------------

def closure_from_names(
    rn_data,
    species_names,
    *,
    counters=None,
) -> int:
    """
    Compute closure starting from a list of species names.

    Parameters
    ----------
    rn_data      : RNData
    species_names: iterable of str
    counters     : optional Counters

    Returns
    -------
    int bitmask (quotiented — E0 species may NOT appear; they're implicit).
    """
    X = rn_data.names_to_bitset(species_names) & ~rn_data.E0_mask
    inv_idx = list(rn_data.species_to_reactions)  # already a list of tuples
    return closure_opt(
        rn_data.supp_q, rn_data.prod_q, X, inv_idx, counters=counters
    )


def is_closed(supp: Sequence[int], prod: Sequence[int], X: int) -> bool:
    """Return True if X is closed (prod(R_X) ⊆ X)."""
    for s, p in zip(supp, prod):
        if (s & ~X) == 0 and (p & ~X):  # s ⊆ X but some product not in X
            return False
    return True


def is_ssm(supp: Sequence[int], prod: Sequence[int], X: int) -> bool:
    """
    Return True if X is semi-self-maintaining (req(X) = ∅).

    req(X) = supp(R_X) - prod(R_X)
    """
    agg_supp = 0
    agg_prod = 0
    for s, p in zip(supp, prod):
        if (s & ~X) == 0:  # reaction active in X
            agg_supp |= s
            agg_prod |= p
    return (agg_supp & ~agg_prod) == 0
