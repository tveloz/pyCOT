"""
decomp — Computational implementation of Veloz's Decomposition Theorem.

Given a semi-organization X (closed, req(X)=0), decomposes it as

    X = (E ∪ F) ∪ D1 ∪ ... ∪ Dm

  E  : catalysts — zero net stoichiometric effect in every reaction of R_X.
  F  : the maximal simultaneously overproducible set (flux-cone witness).
  Di : fragile circuits — equivalence classes of C = X \\ (E ∪ F) under
       dynamical connection (species chained through reactions that draw
       at least one reactant from C, never routing through F).

Theorem 2.16 (Veloz & Razeto-Barry 2017 / decomposing_RAF_v2): X is
self-maintaining iff every Di is self-maintaining w.r.t. its own path
R*_i (reactions of R_X that consume a species of Di).

Public API
----------
decomp.core.decompose(sp_mask, rn, rn_data) -> DecompositionResult
    Standalone decomposition of one semi-organization.
decomp.hierarchy.decompose_hierarchy(so_lattice, rn, rn_data) -> dict[int, DecompositionResult]
    Incremental decomposition of every SO in an EPM/ESPM lattice, reusing
    monotonicity of F (Prop 2.19) and unchanged fragile circuits to skip
    recomputation where provably safe.
"""
