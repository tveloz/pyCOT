"""
dependency.py — Food-dependency structure among a decomposed SO's fragile
circuits, adapting decomposing_RAF_v3.tex's (subsec:dag) interface / depth /
dependency-DAG machinery directly to COT's own E/F/circuit decomposition
(Veloz & Razeto-Barry 2017b) -- no RAF layer involved, per the paper's own
Definitions 10-12 but read against Theorem 1 (Complexity 2022) instead of
against a CRS translation.

Native food = rn_data.E0_mask (species with an inflow reaction -- "natively
overproducible" in the paper's terms). A fragile circuit Di's interface
I(Di) is the subset of E ∪ F that its own path R*_i draws on as a reactant.
Some interface species are native food; the rest are "contextually
overproduced" -- surplus that exists only because ANOTHER circuit's path
(or the shared food reactions) produces it as a byproduct. The dependency
DAG records edges Di ≻ Dj: "Di's path produces a contextually-overproduced
species that Dj's interface needs."

What this DOES explain: the qualitative precondition structure between
circuits -- which one's activity enables which other's inputs -- and it
catches a failure mode Theorem 2.16's LOCAL check cannot see by
construction: two circuits mutually needing each other's surplus before
either can run (the paper's activatability caveat, Remark rmk:reading's
"mutual supply" case). That shows up here as a genuine graph cycle among
circuits and is flagged explicitly as non-activatable, rather than silently
assigned an arbitrary order.

What this does NOT explain: why a circuit fails Theorem 2.16's own local
self-maintenance LP (circuits.py's minimize_sv on S_i). That check is
already fully local to Di's own stoichiometry and does not care about food
availability at all -- Theorem 1 (Complexity 2022) decomposes X's
self-maintenance into independent per-circuit checks, each one free to
assume E ∪ F is available in whatever amount needed. A circuit can have a
perfectly satisfiable interface (reachable at depth 1, straight from native
food) and still fail self-maintenance outright, because its OWN internal
stoichiometry cannot be balanced by any nonnegative combination of its own
path reactions (see explain_circuit_failure below for that separate
question).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .types import DecompositionResult, FragileCircuit


def _mask_to_indices(mask: int) -> list[int]:
    out = []
    m = mask
    while m:
        lsb = m & (-m)
        out.append(lsb.bit_length() - 1)
        m &= m - 1
    return out


def interface(circuit: FragileCircuit, rn_data, food_like_mask: int) -> int:
    """I(Di): species of E ∪ F that Di's own path R*_i consumes as a reactant."""
    mask = 0
    for r in circuit.reaction_ids:
        mask |= rn_data.supp_raw[r]
    return mask & food_like_mask


def _produces(rn_data, reaction_ids, species_mask: int) -> int:
    """Species of `species_mask` that appear as a PRODUCT of some reaction
    in `reaction_ids` (relational "can supply", not a net-positivity claim)."""
    mask = 0
    for r in reaction_ids:
        mask |= rn_data.prod_raw[r]
    return mask & species_mask


@dataclass
class DependencyDAG:
    """
    Dependency structure among the fragile circuits of one decomposed SO.

    circuits            : same order as DecompositionResult.circuits.
    native_food_mask     : E0 ∩ (E ∪ F) -- natively overproducible interface
                            (has its own inflow reaction).
    food_reaction_mask   : (E ∪ F) \\ E0 species produced by a "food
                            reaction" -- one whose support touches NO
                            circuit species at all, so it runs independent
                            of any circuit's own activity (R_F in the
                            paper's notation, e.g. the grass/grain
                            amplifiers here). Available on the same footing
                            as native food for dependency purposes, even
                            though it isn't literally an inflow.
    contextual_food_mask : (E ∪ F) \\ E0 -- overproduced only as a byproduct
                            of SOME reaction (food or circuit).
    interfaces           : circuit index -> I(Di) mask.
    edges                : (i, j) meaning Di ≻ Dj (Di's path supplies part
                            of Dj's interface that neither native food nor
                            a food reaction covers).
    depth                : circuit index -> 1 (interface fully covered by
                            native food + food reactions) or 1 + max(depth
                            of its producer circuits); circuits inside a
                            cycle (see `cycles`) get depth None.
    cycles               : list of circuit-index tuples, each a maximal set
                            of mutually-dependent circuits (size > 1) --
                            i.e. NOT activatable in isolation from each
                            other (Assumption 4 of the RAF paper fails for
                            this SO). Empty in the ordinary/acyclic case.
    unmet                : interface species covered by NEITHER food (native
                            or food-reaction) NOR any circuit's path --
                            should be empty for any SO that is genuinely
                            semi-self-maintaining (Def. 3); a non-empty
                            result signals a bug upstream, not a modeling
                            fact, and is surfaced rather than silently
                            ignored.
    """
    circuits: list[FragileCircuit]
    native_food_mask: int
    contextual_food_mask: int
    food_reaction_mask: int = 0
    interfaces: dict[int, int] = field(default_factory=dict)
    edges: list[tuple[int, int]] = field(default_factory=list)
    depth: dict[int, int | None] = field(default_factory=dict)
    cycles: list[tuple[int, ...]] = field(default_factory=list)
    unmet: dict[int, int] = field(default_factory=dict)

    def summary(self, rn_data) -> str:
        lines = []
        for i, c in enumerate(self.circuits):
            iface = self.interfaces.get(i, 0)
            native = iface & self.native_food_mask
            from_food_rxn = iface & self.food_reaction_mask
            ctx = iface & self.contextual_food_mask & ~self.food_reaction_mask
            d = self.depth.get(i)
            d_str = str(d) if d is not None else "CYCLIC (non-activatable alone)"
            lines.append(
                f"  D{i+1} {{{', '.join(rn_data.bitset_to_names(c.species_mask))}}} "
                f"depth={d_str}"
            )
            if native:
                lines.append(f"       native food used:       {rn_data.bitset_to_names(native)}")
            if from_food_rxn:
                lines.append(f"       food-reaction supplied: {rn_data.bitset_to_names(from_food_rxn)}")
            if ctx:
                lines.append(f"       circuit-supplied:       {rn_data.bitset_to_names(ctx)}")
            unmet_i = self.unmet.get(i, 0)
            if unmet_i:
                lines.append(f"       ** UNMET interface **: {rn_data.bitset_to_names(unmet_i)}")
        if self.edges:
            lines.append("  dependency edges (source → target, source's path supplies target's interface):")
            for i, j in self.edges:
                lines.append(f"    D{i+1} → D{j+1}")
        if self.cycles:
            for cyc in self.cycles:
                names = ", ".join(f"D{i+1}" for i in cyc)
                lines.append(f"  ** UNSOLVABLE DEPENDENCY **: {{{names}}} mutually require "
                              f"each other's surplus -- none can activate first")
        return "\n".join(lines) if lines else "  (no fragile circuits)"


def build_dependency_dag(result: DecompositionResult, rn_data) -> DependencyDAG:
    food_like = result.E_mask | result.F_mask
    native = rn_data.E0_mask & food_like
    contextual = food_like & ~rn_data.E0_mask
    circuit_species = 0
    for c in result.circuits:
        circuit_species |= c.species_mask

    dag = DependencyDAG(circuits=result.circuits, native_food_mask=native,
                         contextual_food_mask=contextual)

    ifaces = {i: interface(c, rn_data, food_like) for i, c in enumerate(result.circuits)}
    dag.interfaces = ifaces

    # "Food reactions" (R_F in the paper's notation): reactions of R_X that
    # don't touch ANY circuit species at all. These run independent of any
    # circuit's own activity, so whatever contextual-food species they
    # produce is available at depth 1, on the same footing as native food
    # -- exactly like R6 (grass+water+fertilizer -> 2grass) or R8 (the
    # grain amplifier) here: neither reaction depends on a fragile circuit,
    # so a circuit needing grass or grain does NOT thereby depend on
    # another circuit, even though grass/grain are themselves "contextual"
    # (non-native) food.
    food_reactions = [r for r in result.R_X if (rn_data.supp_raw[r] & circuit_species) == 0]
    food_produced = 0
    for r in food_reactions:
        food_produced |= rn_data.prod_raw[r]
    food_reaction_mask = food_produced & contextual
    food_available = native | food_reaction_mask
    dag.food_reaction_mask = food_reaction_mask

    # Which circuits can supply a given (still-unresolved) contextual-food species.
    producers_of: dict[int, list[int]] = {}
    for i, c in enumerate(result.circuits):
        supplied = _produces(rn_data, c.reaction_ids, contextual)
        for sp in _mask_to_indices(supplied):
            producers_of.setdefault(sp, []).append(i)

    edges: set[tuple[int, int]] = set()
    unmet: dict[int, int] = {}
    for j, iface in ifaces.items():
        ctx_needed = iface & contextual & ~food_available
        for sp in _mask_to_indices(ctx_needed):
            srcs = producers_of.get(sp, [])
            if not srcs:
                unmet[j] = unmet.get(j, 0) | (1 << sp)
                continue
            for i in srcs:
                if i != j:
                    edges.add((i, j))
    dag.edges = sorted(edges)
    dag.unmet = unmet

    # Depth via fixed-point / cycle detection (Tarjan-lite): a circuit whose
    # interface is entirely food-available (native OR produced by a pure
    # food reaction) is depth 1; otherwise 1 + max(depth of its circuit
    # producers), once ALL producers have a settled depth. Circuits that
    # never settle (because they only depend on each other) form the
    # cycles list instead.
    children: dict[int, set[int]] = {i: set() for i in ifaces}   # i -> circuits i depends on
    for i, j in edges:
        children[j].add(i)

    depth: dict[int, int] = {}
    remaining = set(ifaces.keys())
    for i in list(remaining):
        if (ifaces[i] & contextual & ~food_available) == 0 or not children[i]:
            depth[i] = 1
            remaining.discard(i)

    changed = True
    while changed and remaining:
        changed = False
        for i in list(remaining):
            deps = children[i]
            if deps <= depth.keys():
                depth[i] = 1 + max((depth[d] for d in deps), default=0)
                remaining.discard(i)
                changed = True

    dag.depth = {i: depth.get(i) for i in ifaces}

    if remaining:
        # Strongly-connected components among the leftover (mutually
        # dependent) circuits -- report each nontrivial SCC as one
        # unsolvable-dependency block.
        dag.cycles = _tarjan_sccs({i: children[i] & remaining for i in remaining})

    return dag


def _tarjan_sccs(graph: dict[int, set[int]]) -> list[tuple[int, ...]]:
    index_counter = [0]
    stack: list[int] = []
    lowlink: dict[int, int] = {}
    index: dict[int, int] = {}
    on_stack: dict[int, bool] = {}
    result: list[tuple[int, ...]] = []

    def strongconnect(v):
        index[v] = index_counter[0]
        lowlink[v] = index_counter[0]
        index_counter[0] += 1
        stack.append(v)
        on_stack[v] = True
        for w in graph.get(v, ()):
            if w not in index:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif on_stack.get(w):
                lowlink[v] = min(lowlink[v], index[w])
        if lowlink[v] == index[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            if len(comp) > 1:
                result.append(tuple(sorted(comp)))

    for v in graph:
        if v not in index:
            strongconnect(v)
    return result


def explain_circuit_failure(circuit: FragileCircuit, domain, rn_data) -> str:
    """
    For a circuit that FAILED Theorem 2.16's local self-maintenance LP:
    print its own stoichiometric submatrix (species x path reactions) and
    flag any species with no strictly-positive entry in its row at all --
    an immediate, LP-independent proof that no nonnegative combination of
    Di's own path reactions can ever net-produce it (a "structurally
    starved" species, the same phenomenon as Example ex:gap's 2x->y,y->x
    but caught relationally instead of numerically).

    This is purely diagnostic (the LP result is already authoritative);
    it exists to make an infeasible circuit's OWN arithmetic legible rather
    than leaving "is_self_maintaining=False" unexplained.
    """
    if circuit.is_self_maintaining:
        return "  (self-maintaining -- no failure to explain)"

    di_species = _mask_to_indices(circuit.species_mask)
    r_star = circuit.reaction_ids
    if not r_star:
        names = ", ".join(rn_data.bitset_to_names(circuit.species_mask))
        return (f"  no reaction in R_X consumes any of {{{names}}} -- "
                f"this circuit has no path at all and can never be produced.")

    rows = [domain.row_of(sp) for sp in di_species]
    cols = [domain.col_of(r) for r in r_star]
    S_i = domain.S[np.ix_(rows, cols)]

    sp_names = [rn_data.species_name(sp) for sp in di_species]
    r_names = [rn_data.reaction_name(r) for r in r_star]

    lines = [f"  path R*_i = {{{', '.join(r_names)}}}"]
    starved = []
    for k, sp in enumerate(di_species):
        row = S_i[k, :]
        if (row > 1e-9).sum() == 0:
            starved.append(sp_names[k])
    if starved:
        lines.append(f"  ** structurally starved **: {{{', '.join(starved)}}} "
                      f"never appears with a positive coefficient in any "
                      f"reaction of its own path -- no flux on R*_i can net-produce "
                      f"it, so self-maintenance is impossible regardless of rates.")
    else:
        lines.append("  every species has at least one net-positive producer within "
                      "the path, but no nonnegative rate assignment balances all of "
                      "them simultaneously (an infeasible ratio constraint, as in "
                      "Example ex:gap: 2x→y, y→x forces v2≥2v1 and v1≥v2).")
    lines.append(f"  stoichiometry (rows=species, cols=reactions):")
    lines.append(f"    {'':12s} " + " ".join(f"{n:>10s}" for n in r_names))
    for k, name in enumerate(sp_names):
        lines.append(f"    {name:12s} " + " ".join(f"{S_i[k, c]:10.2f}" for c in range(len(r_names))))
    return "\n".join(lines)
