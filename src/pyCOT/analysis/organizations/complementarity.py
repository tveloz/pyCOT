"""
complementarity.py — ERC complementarity per Veloz & Bassi (paper Defs 22–26).

Theory
------
Supply from E to E' (Def 22):
  supl(E ⇀ E') = prod(R_E) ∩ req(E')
  In bitset: prod_mask_i & req_mask_j

Basic complementarity (Def 23):
  E and E' are complementary if supl(E⇀E') ∪ supl(E'⇀E) ≠ ∅.
  One ERC produces a species the other needs externally.

  This is NOT restricted to incomparable pairs. For comparable E ⊊ E':
  prod(R_E) ⊆ prod(R_E') (every reaction triggered within closure(E) is
  also triggered within closure(E'), since R_E ⊆ R_E'), so
  supl(E⇀E') = prod(R_E) ∩ req(E') = ∅ — the smaller ERC can never supply
  the larger one. But the converse is NOT guaranteed empty: supl(E'⇀E) =
  prod(R_E') ∩ req(E) can be nonempty, because req(E) is computed from
  E's own (smaller) reaction set R_E and is not generally a subset of
  req(E') or disjoint from prod(R_E'). Concretely: r0: a=>2a, r1: a+b=>c,
  r2: d=>b, r3: d=>a+d gives E=closure({a,b})={a,b,c} with req(E)={b},
  and E'=closure({d})={a,b,c,d} ⊋ E with prod(R_E')∋b — so E' supplies b
  to E even though E ⊊ E'. This is "intra-chain" complementarity, as
  opposed to the "inter-chain" case where E and E' are incomparable.
  Each CompPair/FundComp carries a `chain` flag recording which case it is.

Pure complementarity (Def 24):
  Basic and NOT synergetic (pass syn_result to compute this).

Fundamental complementarity (Def 26):
  For species s ∈ supl(E ⇀ E'), the supply is fundamental if:
    E  ∈ minprod(s)  — no ERC E'' ⊊ E also produces s
    E' ∈ mincons(s)  — no ERC E'' ⊊ E' also requires s
  where:
    minprod(s) = ⊆-minimal {E ∈ ε : s ∈ prod(R_E)}
    mincons(s) = ⊆-minimal {E ∈ ε : s ∈ req(E)}

Public API
----------
compute_complementarities(ercs, hier, syn_result=None, *, counters=None) -> CompResult

CompPair(i, j, fwd_supply, bwd_supply, is_pure, chain)
  i, j         : ERC indices (i < j)
  fwd_supply   : prod_i & req_j  (species E_i supplies to E_j, as bitmask)
  bwd_supply   : prod_j & req_i  (species E_j supplies to E_i, as bitmask)
  is_pure      : True if pair is basic-complementary and NOT synergetic
  chain        : True if i and j are hierarchy-comparable ("intra-chain");
                 False if incomparable ("inter-chain"). For chain pairs,
                 exactly one of fwd_supply/bwd_supply is guaranteed 0 (the
                 subset-ERC-to-superset-ERC direction).

FundComp(prod_idx, cons_idx, species, chain)
  prod_idx : ERC index of the ⊆-minimal producer of `species`
  cons_idx : ERC index of the ⊆-minimal consumer of `species`
  species  : bit index (0-based) of the supplied species
  chain    : True if prod_idx and cons_idx are hierarchy-comparable

CompResult(.basic, .pure, .fundamental)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CompPair:
    """One basic complementary pair of ERCs (i < j)."""
    i: int           # ERC index (smaller)
    j: int           # ERC index (larger)
    fwd_supply: int  # prod_i & req_j — species E_i supplies to E_j (bitmask)
    bwd_supply: int  # prod_j & req_i — species E_j supplies to E_i (bitmask)
    is_pure: bool = False  # True if pair is not synergetic
    chain: bool = False    # True if i, j are hierarchy-comparable (intra-chain)


@dataclass(frozen=True)
class FundComp:
    """One fundamental complementarity: minimal producer supplies s to minimal consumer."""
    prod_idx: int   # index of the ⊆-minimal ERC that produces `species`
    cons_idx: int   # index of the ⊆-minimal ERC that requires `species`
    species:  int   # species bit index (0-based position)
    chain: bool = False  # True if prod_idx, cons_idx are hierarchy-comparable


@dataclass
class CompResult:
    """Aggregate output of compute_complementarities."""
    basic:       list[CompPair] = field(default_factory=list)
    pure:        list[CompPair] = field(default_factory=list)
    fundamental: list[FundComp] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.basic)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _bits(mask: int):
    """Yield set bit positions in ascending order."""
    m = mask
    while m:
        lsb = m & (-m)
        yield lsb.bit_length() - 1
        m &= m - 1


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def compute_complementarities(
    ercs,
    hier,
    syn_result=None,
    *,
    counters=None,
) -> CompResult:
    """
    Compute basic, pure, and fundamental complementarities.

    Parameters
    ----------
    ercs       : list[ERCData]  (from compute_ercs)
    hier       : HierarchyData  (from build_hierarchy)
    syn_result : SynergyResult | None
                 If provided, pure complementarity is also computed
                 (basic pairs that are NOT synergetic).
    counters   : optional Counters

    Returns
    -------
    CompResult with .basic, .pure, .fundamental lists.
    """
    masks = [e.species_mask for e in ercs]
    reqs  = [e.req_mask     for e in ercs]
    prods = [e.prod_mask    for e in ercs]
    n     = len(ercs)

    # ── Build synergetic pairs set (for pure complementarity) ─────────────────
    syn_pairs: set[tuple[int, int]] = set()
    if syn_result is not None:
        for s in syn_result.basic:
            syn_pairs.add((min(s.i, s.j), max(s.i, s.j)))

    # ── Basic + pure complementarity ──────────────────────────────────────────
    basic_pairs: list[CompPair] = []

    for i, j in combinations(range(n), 2):
        # Comparable ("chain") pairs are included: the subset-ERC-to-superset-ERC
        # direction is always 0 by construction (see module docstring), but the
        # reverse direction can be nonempty — that is intra-chain complementarity.
        chain = hier.is_comparable(i, j)
        fwd = prods[i] & reqs[j]   # E_i produces what E_j needs
        bwd = prods[j] & reqs[i]   # E_j produces what E_i needs
        if not (fwd or bwd):
            continue
        # Chain pairs can never be synergetic (synergy requires incomparable
        # ERCs), so every intra-chain complementary pair is trivially pure.
        is_pure = chain or ((syn_result is not None) and ((i, j) not in syn_pairs))
        basic_pairs.append(CompPair(i=i, j=j, fwd_supply=fwd, bwd_supply=bwd, is_pure=is_pure, chain=chain))

    pure_pairs = [cp for cp in basic_pairs if cp.is_pure] if syn_result is not None else []

    # ── Fundamental complementarity ───────────────────────────────────────────
    all_masks = masks + reqs + prods
    max_bit = max((m.bit_length() for m in all_masks), default=0)

    fundamental: list[FundComp] = []

    for s_bit in range(max_bit):
        s_mask = 1 << s_bit

        # ERCs that produce species s
        producers = [i for i in range(n) if prods[i] & s_mask]
        # ERCs that require species s externally
        consumers = [i for i in range(n) if reqs[i] & s_mask]

        if not producers or not consumers:
            continue

        # minprod(s): producers with no strictly contained producer
        min_prods = [
            i for i in producers
            if not any(
                j != i
                and (masks[j] & masks[i]) == masks[j]   # masks[j] ⊆ masks[i]
                and masks[j] != masks[i]                  # strict subset
                for j in producers
            )
        ]

        # mincons(s): consumers with no strictly contained consumer
        min_cons = [
            i for i in consumers
            if not any(
                j != i
                and (masks[j] & masks[i]) == masks[j]
                and masks[j] != masks[i]
                for j in consumers
            )
        ]

        for pi in min_prods:
            for ci in min_cons:
                if pi != ci:
                    fundamental.append(FundComp(
                        prod_idx=pi, cons_idx=ci, species=s_bit,
                        chain=hier.is_comparable(pi, ci),
                    ))

    result = CompResult(basic=basic_pairs, pure=pure_pairs, fundamental=fundamental)

    if counters is not None:
        counters.inc("comp.n_basic",       len(basic_pairs))
        counters.inc("comp.n_pure",        len(pure_pairs))
        counters.inc("comp.n_fundamental", len(fundamental))
        counters.inc("comp.n_basic_intra_chain", sum(1 for cp in basic_pairs if cp.chain))
        counters.inc("comp.n_basic_inter_chain", sum(1 for cp in basic_pairs if not cp.chain))
        counters.inc("comp.n_fundamental_intra_chain", sum(1 for fc in fundamental if fc.chain))
        counters.inc("comp.n_fundamental_inter_chain", sum(1 for fc in fundamental if not fc.chain))

    return result
