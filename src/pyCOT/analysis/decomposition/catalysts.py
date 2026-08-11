"""
catalysts.py — Compute E, the catalyst set.

A species is a catalyst iff it has zero net stoichiometric effect in EVERY
reaction of R_X it appears in (Def. 2.3-ish of the Decomposition Theorem:
"the reaction leaves the species' amount unchanged"). Equivalently: its row
in S (restricted to R_X) is the all-zero row.

Native food (E0) species can never satisfy this: they are produced net-
positive by at least one inflow reaction (that is how E0 is defined/seeded
in cot_gen — see io_pyCOT._inflow_product_mask), so their row is never
all-zero. E and F0 are therefore automatically disjoint; no priority rule
is needed between catalysts.py and overproduction.py.
"""
from __future__ import annotations

import numpy as np

from .bridge import SODomain

_TOL = 1e-9


def compute_catalysts(domain: SODomain) -> int:
    """Return the bitmask of catalyst species within domain.X_full_mask."""
    e_mask = 0
    S = domain.S
    for i, sp in enumerate(domain.sp_indices):
        row = S[i, :]
        if row.size == 0 or np.all(np.abs(row) <= _TOL):
            e_mask |= 1 << sp
    return e_mask
