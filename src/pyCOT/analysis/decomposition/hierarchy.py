"""
hierarchy.py — Incremental decomposition over an entire EPM/ESPM lattice.

Walks a cot_gen SOLattice (cot_gen/deep_report.py: nodes = SO species
masks, order_of = EPM/ESPM order, parents_of = immediate sub-SOs in the
Hasse diagram) bottom-up, reusing two structural guarantees instead of
re-deriving everything from scratch at every node:

  1. F monotonicity (Prop 2.19): if s in F(parent) and no reaction newly
     activated in the child (R_X(child) \\ R_X(parent)) net-consumes s,
     then s in F(child) — no LP needed, just a sign check on the new
     reactions' column(s) for row s.

  2. Fragile-circuit invariance: dynamical connection is recomputed fresh
     for the child (it's cheap — union-find, no LP), but if a resulting
     circuit Di is IDENTICAL (same species AND same R*_i) to a circuit
     already verified at some parent, its self-maintenance verdict is
     reused directly rather than re-solving the LP — this is exactly
     Theorem 2.16's point: R*_i-local self-maintenance only depends on
     Di's own species/reactions, so if those didn't change, the verdict
     can't have changed either.

Species/circuits that don't satisfy either shortcut fall back to the
regular per-node LP work from catalysts.py / overproduction.py / circuits.py
— this module never trades correctness for speed, only skips work that is
PROVABLY redundant.
"""
from __future__ import annotations

from .bridge import so_domain
from .catalysts import compute_catalysts
from .overproduction import compute_overproduced
from .circuits import dynamical_components, build_circuit
from .types import DecompositionResult, FragileCircuit

_TOL = 1e-9


def _mask_to_indices(mask: int) -> list[int]:
    out = []
    m = mask
    while m:
        lsb = m & (-m)
        out.append(lsb.bit_length() - 1)
        m &= m - 1
    return out


def _carry_forward_f(domain, e_mask: int, parent_results: list[DecompositionResult]) -> int:
    """
    Species carried into F(child) for free via Prop 2.19, checked against
    EVERY parent independently (a species only needs ONE valid parent).
    """
    carried = 0
    for parent in parent_results:
        parent_R_X = set(parent.R_X)
        new_r = [r for r in domain.R_X if r not in parent_R_X]
        for sp in _mask_to_indices(parent.F_mask):
            if (carried >> sp) & 1 or (e_mask >> sp) & 1:
                continue
            row = domain.row_of(sp)
            if row < 0:
                continue
            if all(domain.S[row, domain.col_of(r)] >= -_TOL for r in new_r):
                carried |= 1 << sp
    return carried


def _reused_or_fresh_circuit(domain, rn_data, di_mask: int,
                              parent_results: list[DecompositionResult]) -> FragileCircuit:
    r_star = tuple(sorted(
        r for r in domain.R_X if (rn_data.supp_raw[r] & di_mask) != 0
    ))
    for parent in parent_results:
        for c in parent.circuits:
            if c.species_mask == di_mask and c.reaction_ids == r_star:
                return c  # identical circuit already verified — reuse verdict
    return build_circuit(domain, rn_data, di_mask)


def decompose_hierarchy(so_lattice, rn_data, S_full, *, verbose: bool = False):
    """
    Decompose every SO in `so_lattice`, bottom-up.

    Returns dict[sp_mask -> DecompositionResult].
    """
    results: dict[int, DecompositionResult] = {}
    nodes_sorted = sorted(so_lattice.nodes,
                           key=lambda sp: (so_lattice.order_of.get(sp, 0), sp))

    n_lp_skipped_f = 0
    n_lp_skipped_circuit = 0

    for sp in nodes_sorted:
        domain = so_domain(sp, rn_data, S_full)
        e_mask = compute_catalysts(domain)

        parent_masks = so_lattice.parents_of.get(sp, [])
        parent_results = [results[p] for p in parent_masks if p in results]

        carried_f = _carry_forward_f(domain, e_mask, parent_results) if parent_results else 0
        remaining_skip = e_mask | carried_f
        f_mask_new, _w = compute_overproduced(domain, rn_data.E0_mask, remaining_skip)
        f_mask = f_mask_new | carried_f
        n_lp_skipped_f += bin(carried_f).count('1')

        c_mask = domain.X_full_mask & ~(e_mask | f_mask)
        components = dynamical_components(domain, rn_data, c_mask)
        circuits = []
        for comp in components:
            circ = _reused_or_fresh_circuit(domain, rn_data, comp, parent_results)
            if parent_results and any(
                circ is c for p in parent_results for c in p.circuits
            ):
                n_lp_skipped_circuit += 1
            circuits.append(circ)

        is_org = all(c.is_self_maintaining for c in circuits)
        results[sp] = DecompositionResult(
            sp_mask=sp, E_mask=e_mask, F_mask=f_mask, circuits=circuits,
            R_X=tuple(domain.R_X), is_organization=is_org,
        )

        if verbose:
            print(f"  order {so_lattice.order_of.get(sp, 0):>2}  "
                  f"{results[sp].summary()}")

    if verbose:
        print(f"[decompose_hierarchy] carried-forward F species: {n_lp_skipped_f}, "
              f"reused circuit verdicts: {n_lp_skipped_circuit} / {len(nodes_sorted)} nodes")

    return results
