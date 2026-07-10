"""
synergy.py — Bitset-based ERC synergy detection (Stage 2).

Theory (from Veloz & Bassi §3–4)
---------------------------------
Given ERCs E1, E2, ET (with ET not contained in E1 or E2):

  Basic synergy (E1, E2) → ET:
    ∃ b ∈ MinBas(ET) such that
        b ⊆ (E1.mask | E2.mask)   (joint closure covers the basis)
      ∧ b ⊄ E1.mask               (E1 alone doesn't cover it)
      ∧ b ⊄ E2.mask               (E2 alone doesn't cover it)

  Maximal synergy: among all basic synergies (E1,E2)→ET,
    keep only those ET that are not contained in another target ET'
    of a basic synergy with the same pair (E1, E2).

  Fundamental synergy: a maximal synergy (E1,E2)→ET is fundamental if
    NO pair (E1', E2') with E1' ⊆ E1, E2' ⊆ E2 (E1'≠E1 or E2'≠E2)
    has a maximal synergy to ET.

Public API
----------
compute_synergies(ercs, hierarchy, *, level="maximal", counters=None)
    -> SynergyResult (aggregate)

SynergyTuple — named record for one synergy relation

Indices are positions in the `ercs` list (same as HierarchyData).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Sequence


# ---------------------------------------------------------------------------
# SynergyTuple
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SynergyTuple:
    """One synergistic relationship: (i, j) → k."""
    i: int       # index of first ERC  (ercs[i].species_mask)
    j: int       # index of second ERC (ercs[j].species_mask)
    k: int       # index of target ERC (ercs[k].species_mask)
    level: str   # "basic" | "maximal" | "fundamental"


@dataclass
class SynergyResult:
    """Aggregate output of compute_synergies."""
    basic:       list[SynergyTuple] = field(default_factory=list)
    maximal:     list[SynergyTuple] = field(default_factory=list)
    fundamental: list[SynergyTuple] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Core bitset predicate
# ---------------------------------------------------------------------------

def _has_joint_cover(m1: int, m2: int, min_bases: list[int]) -> bool:
    """
    Return True if joint mask (m1 | m2) covers at least one min_base
    that neither m1 nor m2 covers alone.
    """
    joint = m1 | m2
    for b in min_bases:
        if (b & joint) == b:          # joint covers this basis
            if (b & m1) != b and (b & m2) != b:  # neither alone covers it
                return True
    return False


def _covers_any(m: int, min_bases: list[int]) -> bool:
    """Return True if mask m covers at least one element of min_bases."""
    return any((b & m) == b for b in min_bases)


# ---------------------------------------------------------------------------
# Basic synergies
# ---------------------------------------------------------------------------

def _basic_synergies_pair(
    i: int,
    j: int,
    masks: list[int],
    min_bases: list[list[int]],
    hier,
) -> list[int]:
    """
    Return list of target indices k where (i, j) has a basic synergy.

    Filters:
      • i and j must be incomparable (neither contains the other)
      • k must not be contained in i or j (target not below either base)
      • joint mask covers a min_base of k that neither i nor j covers alone
    """
    if not hier.can_interact(i, j):
        return []

    m1, m2 = masks[i], masks[j]
    targets = []
    for k, (mk, bk) in enumerate(zip(masks, min_bases)):
        if k == i or k == j:
            continue
        # Skip if target is contained by either ERC
        if k in hier.descendants[i] or k in hier.descendants[j]:
            continue
        if _has_joint_cover(m1, m2, bk):
            targets.append(k)
    return targets


def compute_basic_synergies(
    ercs,
    hier,
    *,
    counters=None,
) -> list[SynergyTuple]:
    """
    Compute all basic synergies.

    O(|E|³ × max|MinBas|) in the worst case.
    In practice much faster because most pairs are comparable (filtered early).
    """
    masks     = [e.species_mask for e in ercs]
    min_bases = [e.min_bases    for e in ercs]
    n = len(ercs)

    basic = []
    pairs_checked = 0
    pairs_comparable = 0

    for i, j in combinations(range(n), 2):
        pairs_checked += 1
        targets = _basic_synergies_pair(i, j, masks, min_bases, hier)
        if not targets and not hier.can_interact(i, j):
            pairs_comparable += 1
        for k in targets:
            basic.append(SynergyTuple(i=i, j=j, k=k, level="basic"))

    if counters is not None:
        counters.inc("syn.pairs_checked",    pairs_checked)
        counters.inc("syn.pairs_comparable", pairs_comparable)
        counters.inc("syn.basic_count",      len(basic))

    return basic


# ---------------------------------------------------------------------------
# Maximal synergies
# ---------------------------------------------------------------------------

def compute_maximal_synergies(
    basic: list[SynergyTuple],
    ercs,
    hier,
    *,
    counters=None,
) -> list[SynergyTuple]:
    """
    Filter basic synergies to keep only maximal targets per (i, j) pair.

    Maximal: target k is maximal for pair (i,j) iff no other target k'
    for the same pair satisfies ercs[k].mask ⊊ ercs[k'].mask.

    Returns new SynergyTuple objects with level="maximal".
    """
    masks = [e.species_mask for e in ercs]

    # Group by pair
    from collections import defaultdict
    pair_targets: dict[tuple[int,int], list[int]] = defaultdict(list)
    for s in basic:
        pair_targets[(s.i, s.j)].append(s.k)

    maximal = []
    for (i, j), ks in pair_targets.items():
        for k in ks:
            # k is maximal iff no other k' in ks has masks[k] ⊊ masks[k']
            mk = masks[k]
            dominated = any(
                k2 != k and (mk & masks[k2]) == mk and mk != masks[k2]
                for k2 in ks
            )
            if not dominated:
                maximal.append(SynergyTuple(i=i, j=j, k=k, level="maximal"))

    if counters is not None:
        counters.inc("syn.maximal_count", len(maximal))

    return maximal


# ---------------------------------------------------------------------------
# Fundamental synergies
# ---------------------------------------------------------------------------

def compute_fundamental_synergies(
    maximal: list[SynergyTuple],
    ercs,
    hier,
    *,
    counters=None,
) -> list[SynergyTuple]:
    """
    Filter maximal synergies to keep only fundamental ones.

    A maximal synergy (i, j) → k is fundamental iff no pair (i', j')
    with i' in descendants[i]∪{i}, j' in descendants[j]∪{j}, (i',j')≠(i,j)
    has a maximal synergy to the same target k.

    This is O(|maximal| × (avg|desc|)²) which can be expensive for large
    hierarchies.  The bitset approach is much faster than the name-based one
    because comparing pair memberships uses set lookups not closure calls.
    """
    # Index maximal by target k  →  set of (i, j) pairs that synergize to k
    from collections import defaultdict
    target_to_pairs: dict[int, set[tuple[int, int]]] = defaultdict(set)
    for s in maximal:
        target_to_pairs[s.k].add((s.i, s.j))

    fundamental = []
    for s in maximal:
        i, j, k = s.i, s.j, s.k
        pairs_to_k = target_to_pairs[k]

        # Sub-pairs: i' ∈ desc[i]∪{i}, j' ∈ desc[j]∪{j}, (i',j') ≠ (i,j)
        i_set = hier.descendants[i] | {i}
        j_set = hier.descendants[j] | {j}

        is_fundamental = True
        for i2 in i_set:
            for j2 in j_set:
                if i2 == i and j2 == j:
                    continue
                # canonical pair (smaller index first)
                pair = (min(i2, j2), max(i2, j2))
                if pair in pairs_to_k:
                    is_fundamental = False
                    break
            if not is_fundamental:
                break

        if is_fundamental:
            fundamental.append(SynergyTuple(i=i, j=j, k=k, level="fundamental"))

    if counters is not None:
        counters.inc("syn.fundamental_count", len(fundamental))

    return fundamental


# ---------------------------------------------------------------------------
# Fused direct computation of fundamental synergies (no intermediate storage)
# ---------------------------------------------------------------------------

def compute_fundamental_direct(
    ercs,
    hier,
    *,
    counters=None,
) -> SynergyResult:
    """
    Compute basic + maximal + fundamental synergies in a single pair-at-a-time
    pass, without materialising the full basic list globally.

    Complexity: O(|E|^2 * |E| * B) — identical to the staged version — but
    the inner work per pair is bounded by |E|*B + |desc_i|*|desc_j|.
    Peak memory: O(|E| + max_targets_per_pair) instead of O(basic_total).

    Algorithm per incomparable pair (i, j):
      1. Collect basic targets k  (joint mask covers a min_base of k).
      2. Filter to maximal k      (no other basic target strictly contains k).
      3. For each maximal k:
           Scan all sub-pairs (i', j') with i'∈desc(i)∪{i}, j'∈desc(j)∪{j}.
           If ANY sub-pair achieves a basic synergy to k → not fundamental.
      4. Emit: basic (all), maximal (filtered), fundamental (filtered).
    """
    masks     = [e.species_mask for e in ercs]
    min_bases = [e.min_bases    for e in ercs]
    n = len(ercs)

    all_basic:       list[SynergyTuple] = []
    all_maximal:     list[SynergyTuple] = []
    all_fundamental: list[SynergyTuple] = []

    pairs_checked = 0
    pairs_comparable = 0

    for i, j in combinations(range(n), 2):
        if not hier.can_interact(i, j):
            pairs_comparable += 1
            continue
        pairs_checked += 1

        m1, m2 = masks[i], masks[j]
        desc_i = hier.descendants[i]
        desc_j = hier.descendants[j]

        # Step 1: collect basic targets
        basic_ks: list[int] = []
        for k in range(n):
            if k == i or k == j:
                continue
            if k in desc_i or k in desc_j:
                continue
            if _has_joint_cover(m1, m2, min_bases[k]):
                basic_ks.append(k)

        for k in basic_ks:
            all_basic.append(SynergyTuple(i=i, j=j, k=k, level="basic"))

        if not basic_ks:
            continue

        # Step 2: maximal targets — k not strictly dominated by another basic_k
        mk_map = {k: masks[k] for k in basic_ks}
        maximal_ks: list[int] = []
        for k in basic_ks:
            mk = mk_map[k]
            if not any(
                k2 != k and (mk & mk_map[k2]) == mk and mk != mk_map[k2]
                for k2 in basic_ks
            ):
                maximal_ks.append(k)

        for k in maximal_ks:
            all_maximal.append(SynergyTuple(i=i, j=j, k=k, level="maximal"))

        # Step 3: fundamental — no sub-pair of (i,j) also achieves k
        i_set = desc_i | {i}
        j_set = desc_j | {j}

        for k in maximal_ks:
            bk = min_bases[k]
            is_fund = True
            for i2 in i_set:
                if not is_fund:
                    break
                m_i2 = masks[i2]
                for j2 in j_set:
                    if i2 == i and j2 == j:
                        continue
                    if _has_joint_cover(m_i2, masks[j2], bk):
                        is_fund = False
                        break
            if is_fund:
                all_fundamental.append(SynergyTuple(i=i, j=j, k=k, level="fundamental"))

    if counters is not None:
        counters.inc("syn.pairs_checked",      pairs_checked)
        counters.inc("syn.pairs_comparable",   pairs_comparable)
        counters.inc("syn.basic_count",        len(all_basic))
        counters.inc("syn.maximal_count",      len(all_maximal))
        counters.inc("syn.fundamental_count",  len(all_fundamental))

    return SynergyResult(
        basic=all_basic,
        maximal=all_maximal,
        fundamental=all_fundamental,
    )


# ---------------------------------------------------------------------------
# Basis-and-target-first algorithm (alternative to pair-first)
# ---------------------------------------------------------------------------

def compute_synergies_basis_first(
    ercs,
    hier,
    *,
    counters=None,
) -> SynergyResult:
    """
    Alternative fundamental-synergy algorithm: iterate targets k and min-bases b,
    collecting basic synergies only from ERC pairs that are *relevant* to each
    basis (i.e., each ERC covers some but not all species of b).

    Rationale
    ---------
    The pair-first algorithm checks every O(|E|²) incomparable pair against
    every target k's min-bases — O(|E|³·B) total.  For sparse networks, most
    ERCs share few species with any given basis, so |relevant_k_b| << |E|
    and the inner pair scan over relevant ERCs is much smaller in practice.

    Correctness
    -----------
    Equivalence to the pair-first algorithm is exact (proved by De Morgan on
    the joint-cover predicate): (i,j)→k is detected iff ∃ b∈MinBas(k) s.t.
      bi = b & masks[i] ≠ ∅, ≠ b  (i covers some but not all of b)
      bj = b & masks[j] ≠ ∅, ≠ b  (j covers some but not all of b)
      bi | bj == b                  (joint coverage is complete)
    which is precisely the basic-synergy predicate.

    The maximal and fundamental filters are reused from the pair-first pipeline
    unchanged.

    Parameters
    ----------
    ercs     : list[ERCData]   (output of compute_ercs)
    hier     : HierarchyData   (output of build_hierarchy)
    counters : optional Counters for instrumentation

    Returns
    -------
    SynergyResult with .basic, .maximal, .fundamental populated.
    """
    masks         = [e.species_mask for e in ercs]
    min_bases_all = [e.min_bases    for e in ercs]
    n             = len(ercs)

    # pair (i,j) → set of basic target indices k
    pair_to_targets: dict[tuple[int, int], set[int]] = {}

    for k in range(n):
        bk = min_bases_all[k]
        if not bk:
            continue

        for b in bk:
            # ERCs that cover SOME but NOT ALL species of b.
            # (ERCs that cover all of b have b ⊆ mask → mask contains the closure
            #  of b = k, so k ⊆ mask, meaning i is an ancestor of k; any basis of
            #  k is fully inside such an ERC → it can trigger k alone → not a
            #  synergy source.  ERCs that cover none of b contribute nothing.)
            relevant: list[int] = []
            for i in range(n):
                bi = b & masks[i]
                if bi == 0 or bi == b:
                    continue
                relevant.append(i)

            # Check all incomparable pairs from relevant
            r = len(relevant)
            for a in range(r):
                i  = relevant[a]
                bi = b & masks[i]
                # Species of b that i misses — j must cover all of them
                rest = b ^ bi    # = b & ~bi, since bi ⊆ b
                for c in range(a + 1, r):
                    j = relevant[c]
                    if not hier.can_interact(i, j):
                        continue
                    if (rest & masks[j]) != rest:
                        continue   # j doesn't cover remaining species
                    # (i, j) jointly covers b; neither alone does → basic synergy to k
                    p = (min(i, j), max(i, j))
                    if p not in pair_to_targets:
                        pair_to_targets[p] = set()
                    pair_to_targets[p].add(k)

    # Build basic synergy list (sorted for determinism)
    all_basic: list[SynergyTuple] = [
        SynergyTuple(i=p[0], j=p[1], k=k, level="basic")
        for p, ks in sorted(pair_to_targets.items())
        for k in sorted(ks)
    ]

    # Reuse existing maximal and fundamental filters
    maximal     = compute_maximal_synergies(all_basic, ercs, hier, counters=counters)
    fundamental = compute_fundamental_synergies(maximal, ercs, hier, counters=counters)

    if counters is not None:
        counters.inc("syn.bf.basic_count",       len(all_basic))
        counters.inc("syn.bf.maximal_count",     len(maximal))
        counters.inc("syn.bf.fundamental_count", len(fundamental))

    return SynergyResult(basic=all_basic, maximal=maximal, fundamental=fundamental)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_synergies(
    ercs,
    hier,
    *,
    level: str = "maximal",
    direct: bool = False,
    basis_first: bool = False,
    counters=None,
) -> SynergyResult:
    """
    Compute ERC synergies up to the requested level.

    Parameters
    ----------
    ercs        : list[ERCData]  (output of compute_ercs)
    hier        : HierarchyData  (output of build_hierarchy)
    level       : "basic" | "maximal" | "fundamental"
    direct      : if True and level="fundamental", use the fused single-pass
                  algorithm that avoids materialising the full basic list.
                  Identical results, lower peak memory.
    basis_first : if True, use the target-and-basis-first algorithm which
                  iterates targets k then min-bases b and only examines ERC
                  pairs relevant to each basis.  Produces identical results;
                  typically faster on large sparse networks.
                  (Implies level="fundamental"; all three levels are returned.)
    counters    : optional Counters for instrumentation

    Returns
    -------
    SynergyResult with .basic, .maximal, .fundamental lists populated
    up to the requested level.
    """
    if basis_first:
        return compute_synergies_basis_first(ercs, hier, counters=counters)

    if level == "fundamental" and direct:
        return compute_fundamental_direct(ercs, hier, counters=counters)

    result = SynergyResult()

    result.basic = compute_basic_synergies(ercs, hier, counters=counters)

    if level in ("maximal", "fundamental"):
        result.maximal = compute_maximal_synergies(
            result.basic, ercs, hier, counters=counters
        )

    if level == "fundamental":
        result.fundamental = compute_fundamental_synergies(
            result.maximal, ercs, hier, counters=counters
        )

    return result
