"""
generators.py — Synergy reachability and primitive ERC identification (Stage 3).

Theory (Veloz & Bassi)
-----------------------
Given a set S of ERCs and a collection of fundamental synergies:

  Synergy closure / reachability:
    Starting from S, repeatedly add any ERC k such that
    some (i,j)→k in fundamental_synergies with i∈S and j∈S.
    The fixpoint is "all ERCs reachable from S via synergy".

  Primitive ERC:
    An ERC k is primitive if it is NOT the target of any fundamental synergy.
    Primitive ERCs cannot be "generated" from smaller ERCs via synergy.

  Generative basis:
    The set of all primitive ERCs.  If all non-primitive ERCs are reachable
    from the basis (via iterated synergy), the basis "generates" the whole
    ERC poset.

Public API
----------
compute_generators(ercs, hier, syn_result, *, counters=None) -> GeneratorResult

GeneratorResult:
  .primitive_indices  — sorted list of primitive ERC indices
  .basis_reach        — frozenset of ERC indices reachable from all primitives
  .coverage           — fraction of ERCs reachable from the basis (1.0 = full)
  .unreachable        — frozenset of ERC indices not reachable from the basis

Helpers:
  reachable_from(seed_indices, syn_tuples) -> frozenset[int]
  primitive_ercs(ercs, syn_result) -> list[int]
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

@dataclass
class GeneratorResult:
    """Output of compute_generators."""
    primitive_indices: list[int]          # sorted ERC indices of primitives
    basis_reach: frozenset[int]           # ERCs reachable from all primitives
    coverage: float                       # fraction of ERCs reachable
    unreachable: frozenset[int]           # ERC indices not reachable from basis

    def is_complete(self) -> bool:
        """True if the generative basis reaches all ERCs."""
        return self.coverage == 1.0


# ---------------------------------------------------------------------------
# Core algorithms
# ---------------------------------------------------------------------------

def reachable_from(
    seed_indices: Iterable[int],
    syn_tuples,                    # list[SynergyTuple] or [(i,j,k,...)]
) -> frozenset[int]:
    """
    Compute the synergy-closure of seed_indices under fundamental synergies.

    Iterates until no new ERC can be added:
      if i∈reach AND j∈reach AND (i,j)→k is a synergy THEN k∈reach.

    Parameters
    ----------
    seed_indices : iterable of int (ERC indices)
    syn_tuples   : iterable of objects with .i, .j, .k attributes
                   (typically SynergyTuple from synergy.py)

    Returns
    -------
    frozenset[int] of all reachable ERC indices (including seed).
    """
    reach: set[int] = set(seed_indices)
    # Pre-extract (i,j,k) triples for speed
    triples = [(s.i, s.j, s.k) for s in syn_tuples]

    changed = True
    while changed:
        changed = False
        for (i, j, k) in triples:
            if i in reach and j in reach and k not in reach:
                reach.add(k)
                changed = True

    return frozenset(reach)


def primitive_ercs(ercs, syn_result) -> list[int]:
    """
    Return sorted list of primitive ERC indices.

    An ERC k is primitive if it is NOT the target (k) of any fundamental synergy.
    Primitive ERCs form the generative basis: they must be "given" externally
    since they cannot be derived by combining other ERCs via synergy.

    Parameters
    ----------
    ercs       : list[ERCData]
    syn_result : SynergyResult (from compute_synergies with level='fundamental')

    Returns
    -------
    Sorted list[int] of primitive ERC indices.
    """
    synergy_targets = {s.k for s in syn_result.fundamental}
    return sorted(i for i in range(len(ercs)) if i not in synergy_targets)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def compute_generators(
    ercs,
    hier,
    syn_result,
    *,
    counters=None,
) -> GeneratorResult:
    """
    Compute primitive ERCs and their generative reach.

    Steps
    -----
    1. Find primitive ERCs (not targets of any fundamental synergy).
    2. Compute reachable_from(primitives, fundamental_synergies).
    3. Measure coverage and identify unreachable ERCs.

    Parameters
    ----------
    ercs       : list[ERCData]
    hier       : HierarchyData   (unused currently; reserved for future filters)
    syn_result : SynergyResult   (must include .fundamental)
    counters   : optional Counters

    Returns
    -------
    GeneratorResult
    """
    n = len(ercs)
    prims = primitive_ercs(ercs, syn_result)

    if n == 0:
        result = GeneratorResult(
            primitive_indices=[],
            basis_reach=frozenset(),
            coverage=1.0,
            unreachable=frozenset(),
        )
    else:
        basis_reach = reachable_from(prims, syn_result.fundamental)
        all_indices = frozenset(range(n))
        unreachable = all_indices - basis_reach
        coverage = len(basis_reach) / n

        result = GeneratorResult(
            primitive_indices=prims,
            basis_reach=basis_reach,
            coverage=coverage,
            unreachable=unreachable,
        )

    if counters is not None:
        counters.inc("gen.n_primitives", len(result.primitive_indices))
        counters.inc("gen.n_reachable",  len(result.basis_reach))
        counters.inc("gen.n_unreachable", len(result.unreachable))

    return result
