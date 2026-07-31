"""
maxraf_decomp.py -- decompose ONLY the maxRAF-generated set X = gen(F0,
maxRAF) under each inflow scenario, instead of enumerating the whole
semi-organization lattice (cot_gen's EPM/ESPM machinery, which does not
scale to genome-size networks or even to a sparse-food e_coli_core
scenario -- see prior conversation).

Everything here is polynomial: maxRAF is polynomial (Hordijk & Steel),
and decompose_species_set() is one LP per fragile circuit, not a
combinatorial lattice search -- this is the SAME pipeline
scripts/run_raf_comparison.py already runs for one scenario; this module
just makes it reusable across many (network, inflow) pairs, and adds a
second piece: building a NESTED CHAIN of organizations directly from the
dependency DAG's depth order, with no additional search at all.

Why the chain is free: Theorem 2.16 says X is self-maintaining iff every
fragile circuit self-maintains on its OWN local path -- a per-circuit,
local LP already solved once for the full maxRAF decomposition. So any
prefix X_t = E ∪ F ∪ (circuits with depth <= t) is self-maintaining
automatically, using flags already computed -- PROVIDED X_t is actually
closed (checked explicitly below, not assumed: growing circuits are not
guaranteed to remain closed as a strict subset, since a reaction of a
later circuit's path could happen to already have all its reactants
available earlier).
"""
from __future__ import annotations

from .crs import gen
from .cot_bridge import cot_translate, closure as cot_closure
from .decomp_shim import decompose_species_set, build_shim
from .dependency import build_dependency_dag
from .raf_algo import compute_maxRAF


def analyze_maxraf_decomposition(crs, food_tokens) -> dict:
    """`crs` is built ONCE per network (either induction mode); `food_tokens`
    overrides crs.food for this one scenario -- cheap, since CRS.food is a
    free-standing field independent of the reaction list itself."""
    crs.food = frozenset(food_tokens) & crs.species
    maxraf = compute_maxRAF(crs)
    X_raf = gen(crs.food, [crs.reactions[n] for n in maxraf]) if maxraf else frozenset(crs.food)

    net_translated = cot_translate(crs)
    result = decompose_species_set(net_translated, X_raf)
    shim, S_full = build_shim(net_translated)

    sp_index = dict(shim.species_index)
    F0_mask = 0
    for s in crs.food:
        if s in sp_index:
            F0_mask |= 1 << sp_index[s]

    dag = build_dependency_dag(result, net_translated, shim, S_full, F0_mask)

    return dict(
        food=set(crs.food), maxraf=maxraf, X_raf=X_raf, result=result,
        shim=shim, net_translated=net_translated, S_full=S_full,
        F0_mask=F0_mask, dag=dag,
    )


def organization_chain(analysis: dict) -> dict:
    """Nested chain X_0 subset X_1 subset ... built by adding one
    dependency-depth's worth of circuits at a time. Each entry reports
    whether it is genuinely closed (checked, not assumed) and hence a
    genuine organization (closed + every included circuit self-maintains)."""
    result = analysis["result"]
    dag = analysis["dag"]
    shim = analysis["shim"]
    net = analysis["net_translated"]

    depths_present = sorted({d for d in dag.depth if d is not None})
    unreachable = [i for i, d in enumerate(dag.depth) if d is None]

    base_mask = result.E_mask | result.F_mask
    chain = []
    included: set[int] = set()
    for t in depths_present:
        included |= {i for i, d in enumerate(dag.depth) if d == t}
        sp_mask = base_mask
        for i in included:
            sp_mask |= result.circuits[i].species_mask
        names = set(shim.bitset_to_names(sp_mask))
        closed_names = set(cot_closure(net, names))
        is_closed = closed_names == names
        all_sm = all(result.circuits[i].is_self_maintaining for i in included)
        chain.append(dict(
            depth=t, n_circuits=len(included), n_species=len(names),
            is_closed=is_closed, extra_if_not_closed=sorted(closed_names - names),
            is_organization=is_closed and all_sm,
            species=names,
        ))
    return dict(chain=chain, unreachable_circuits=unreachable, n_circuits_total=len(result.circuits))
