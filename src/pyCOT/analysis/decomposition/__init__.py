"""
pyCOT.analysis.decomposition — Veloz's Decomposition Theorem, computational form.

Given a semi-organization X (closed, req(X)=0 -- e.g. any node of the
EPM/ESPM lattice produced by pyCOT.analysis.organizations.compute_organizations),
decomposes it as

    X = (E ∪ F) ∪ D1 ∪ ... ∪ Dm

  E  : catalysts — zero net stoichiometric effect in every reaction of R_X.
  F  : the maximal simultaneously overproducible set (flux-cone witness).
  Di : fragile circuits — equivalence classes of C = X \\ (E ∪ F) under
       dynamical connection (species chained through reactions that draw
       at least one reactant from C, never routing through F).

Theorem 2.16 (Veloz & Razeto-Barry 2017 / decomposing_RAF_v2): X is
self-maintaining iff every Di is self-maintaining w.r.t. its own path
R*_i (reactions of R_X that consume a species of Di). Each Di's check
reuses pyCOT.analysis.organizations.self_maintenance.minimize_sv -- the
same LP kernel pyCOT.analysis.organizations.compute_organizations uses to
verify whole organizations, applied here at circuit granularity.

Ported from projects/Decomposition_Theorem/decomp/ (relative imports
unchanged; this package was already decoupled from any cot_gen-specific
representation -- it takes a plain species bitmask + duck-typed rn/rn_data,
so nothing in its own logic needed to change on the move).

Public API
----------
decompose(sp_mask, rn_data, S_full) -> DecompositionResult
    Standalone decomposition of one semi-organization.
    (core.decompose; also available as decomposition.core.decompose)
decompose_hierarchy(so_lattice, rn_data, S_full, *, verbose=False)
        -> dict[int, DecompositionResult]
    Incremental decomposition of every SO in an EPM/ESPM lattice, reusing
    monotonicity of F (Prop 2.19) and unchanged fragile circuits to skip
    recomputation where provably safe.
    (hierarchy.decompose_hierarchy)
build_org_graph(org_results, so_lattice=None) -> OrgGraph
    Containment graph over decomposed organizations, for the Hasse-style
    E/F/circuit visualizations in pyCOT.visualization.decomposition_viz.

Visualization: see pyCOT.visualization.decomposition_viz for
plot_organization_hasse / plot_organization_chains / plot_decomposition_evolution.
"""
from .types import DecompositionResult, FragileCircuit
from .bridge import SODomain, so_domain, build_full_stoich
from .core import decompose
from .hierarchy import decompose_hierarchy
from .catalysts import compute_catalysts
from .overproduction import compute_overproduced, test_overproducible
from .circuits import compute_fragile_circuits, build_circuit, dynamical_components
from .witness import compute_witness, OrganizationWitness
from .org_graph import OrgGraph, build_org_graph, enumerate_root_to_leaf_chains, select_divergent_chains
from .dependency import DependencyDAG, build_dependency_dag, explain_circuit_failure

__all__ = [
    "DecompositionResult", "FragileCircuit",
    "SODomain", "so_domain", "build_full_stoich",
    "decompose", "decompose_hierarchy",
    "compute_catalysts",
    "compute_overproduced", "test_overproducible",
    "compute_fragile_circuits", "build_circuit", "dynamical_components",
    "compute_witness", "OrganizationWitness",
    "OrgGraph", "build_org_graph", "enumerate_root_to_leaf_chains", "select_divergent_chains",
    "DependencyDAG", "build_dependency_dag", "explain_circuit_failure",
]
