"""
overproduction.py — Compute F, the maximal simultaneously overproducible set.

Definition (Sec. 2, decomposing_RAF_v2 / systems-05-00030-v2): species s is
overproducible in X if there is v in the flux cone V(X) = {v >= 0 : v_r > 0
for r in R_X, v_r = 0 otherwise, S v >= 0} with (S v)_s > 0.

Proposition 2.7: F, the set of ALL overproducible species, is itself
SIMULTANEOUSLY overproducible — one witness works for all of them at once.
The proof is constructive and is exactly what this module computes: if
v_1, ..., v_k individually witness s_1, ..., s_k (each S v_i >= 0 with
(S v_i)_{s_i} > 0), then v = v_1 + ... + v_k satisfies S v >= 0 (sum of
nonnegatives) and (S v)_{s_i} >= (S v_i)_{s_i} > 0 for every i (the other
terms only add nonnegative slack). So instead of one big joint LP, we solve
one small LP per candidate species and sum the witnesses — cheaper, and it
doubles as a certificate matching the paper's own proof.

Native food (E0) species are added directly, no LP needed: they are net-
produced by their own inflow reaction by construction (Prop 5.1a: F0 ⊆ F).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linprog

from .bridge import SODomain

_EPS = 1e-7
_MARGIN = 1e-9


def test_overproducible(S: np.ndarray, row: int, n_reactions: int):
    """
    LP: maximize (S v)_row  s.t.  S v >= -margin,  v >= 0,  sum(v) = n_reactions.

    Returns (is_overproducible, v or None, value).
    """
    if n_reactions == 0:
        return False, None, 0.0
    c = -S[row, :]                       # linprog minimizes -> maximize (Sv)_row
    A_ub = -S                            # -S v <= margin  <=>  S v >= -margin
    b_ub = _MARGIN * np.ones(S.shape[0])
    A_eq = np.ones((1, n_reactions))
    b_eq = [float(n_reactions)]
    bounds = [(0, None)] * n_reactions
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                   bounds=bounds, method='highs')
    if not res.success:
        return False, None, 0.0
    value = float(S[row, :] @ res.x)
    return value > _EPS, (res.x if value > _EPS else None), value


def compute_overproduced(domain: SODomain, e0_mask: int, e_mask: int):
    """
    Return (F_mask, witnesses) where witnesses is a dict
    {species_index -> flux vector (over domain.R_X)} for every
    individually-tested species that turned out overproducible
    (native food species are not included — they need no witness LP).
    """
    f_mask = e0_mask & domain.X_full_mask
    witnesses: dict[int, np.ndarray] = {}
    n_r = len(domain.R_X)
    skip = e0_mask | e_mask
    for i, sp in enumerate(domain.sp_indices):
        if (skip >> sp) & 1:
            continue
        ok, v, _value = test_overproducible(domain.S, i, n_r)
        if ok:
            f_mask |= 1 << sp
            witnesses[sp] = v
    return f_mask, witnesses


def joint_witness(domain: SODomain, witnesses: dict[int, np.ndarray]) -> np.ndarray | None:
    """
    Sum the individual witnesses into one joint flux vector demonstrating
    simultaneous overproducibility of every species in `witnesses` at once
    (the constructive content of Proposition 2.7). Returns None if there
    are no witnesses to sum.
    """
    if not witnesses:
        return None
    total = np.zeros(len(domain.R_X))
    for v in witnesses.values():
        total += v
    return total
