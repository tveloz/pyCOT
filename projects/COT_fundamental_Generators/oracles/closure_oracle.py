"""
oracles/closure_oracle.py — Naive fixpoint closure (correctness anchor).

closure_oracle(supp, prod, X) → int

Algorithm: repeat "add products of all active reactions" until stable.
This is the textbook definition; it is O(|R| × iterations) and obviously
correct.  It is the reference implementation that the optimized closure in
cot_gen/closure.py must match on all inputs.
"""
from __future__ import annotations

from typing import Sequence


def closure_oracle(
    supp: Sequence[int],
    prod: Sequence[int],
    X: int,
) -> int:
    """
    Naive fixpoint: add products of all active reactions until stable.

    Parameters
    ----------
    supp : support bitmasks (one per reaction).
    prod : product bitmasks (one per reaction).
    X    : initial species bitmask.

    Returns
    -------
    int : bitmask of the smallest closed set ⊇ X.
    """
    prev = -1
    while prev != X:
        prev = X
        for s, p in zip(supp, prod):
            if (s & X) == s:   # supp(r) ⊆ X  →  reaction is active
                X |= p
    return X


def is_closed_oracle(supp: Sequence[int], prod: Sequence[int], X: int) -> bool:
    """Return True iff X is a fixed point of one fixpoint iteration."""
    for s, p in zip(supp, prod):
        if (s & X) == s and (p & ~X):
            return False
    return True


def is_ssm_oracle(supp: Sequence[int], prod: Sequence[int], X: int) -> bool:
    """
    Return True iff X is semi-self-maintaining.

    req(X) = supp(R_X) - prod(R_X)
    """
    agg_s = 0
    agg_p = 0
    for s, p in zip(supp, prod):
        if (s & X) == s:
            agg_s |= s
            agg_p |= p
    return (agg_s & ~agg_p) == 0
