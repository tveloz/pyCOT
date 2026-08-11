"""
core.py — Standalone decomposition of one semi-organization.

decompose(sp_mask, rn_data, S_full) -> DecompositionResult

    X = (E ∪ F) ∪ D1 ∪ ... ∪ Dm,  is_organization iff every Di self-maintains.

This is the single-SO entry point. For walking an entire EPM/ESPM lattice
incrementally (reusing monotonicity of F and unchanged circuits across the
hierarchy), see hierarchy.py.
"""
from __future__ import annotations

import numpy as np

from .bridge import SODomain, so_domain
from .catalysts import compute_catalysts
from .overproduction import compute_overproduced
from .circuits import compute_fragile_circuits
from .types import DecompositionResult


def decompose(sp_mask: int, rn_data, S_full: np.ndarray) -> DecompositionResult:
    """
    Decompose semi-organization `sp_mask` (cot_gen quotiented species-index
    space) using the real stoichiometry in `S_full` (see bridge.build_full_stoich).
    """
    domain = so_domain(sp_mask, rn_data, S_full)
    e_mask = compute_catalysts(domain)
    f_mask, _witnesses = compute_overproduced(domain, rn_data.E0_mask, e_mask)
    circuits = compute_fragile_circuits(domain, rn_data, e_mask, f_mask)
    is_org = all(c.is_self_maintaining for c in circuits)
    return DecompositionResult(
        sp_mask=sp_mask, E_mask=e_mask, F_mask=f_mask,
        circuits=circuits, R_X=tuple(domain.R_X), is_organization=is_org,
    )
