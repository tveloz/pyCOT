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
  Only incomparable pairs are relevant (for comparable E ⊊ E': prod(R_E) ⊆ prod(R_E'),
  so prod(R_E) ∩ req(E') = ∅ since req(E') = supp(R_E') \\ prod(R_E')).

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

CompPair(i, j, fwd_supply, bwd_supply, is_pure)
  i, j         : ERC indices (i < j)
  fwd_supply   : prod_i & req_j  (species E_i supplies to E_j, as bitmask)
  bwd_supply   : prod_j & req_i  (species E_j supplies to E_i, as bitmask)
  is_pure      : True if pair is basic-complementary and NOT synergetic

FundComp(prod_idx, cons_idx, species)
  prod_idx : ERC index of the ⊆-minimal producer of `species`
  cons_idx : ERC index of the ⊆-minimal consumer of `species`
  species  : bit index (0-based) of the supplied species

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


@dataclass(frozen=True)
class FundComp:
    """One fundamental complementarity: minimal producer supplies s to minimal consumer."""
    prod_idx: int   # index of the ⊆-minimal ERC that produces `species`
    cons_idx: int   # index of the ⊆-minimal ERC that requires `species`
    species:  int   # species bit index (0-based position)


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
        if not hier.can_interact(i, j):   # skip comparable pairs
            continue
        fwd = prods[i] & reqs[j]   # E_i produces what E_j needs
        bwd = prods[j] & reqs[i]   # E_j produces what E_i needs
        if not (fwd or bwd):
            continue
        is_pure = (syn_result is not None) and ((i, j) not in syn_pairs)
        basic_pairs.append(CompPair(i=i, j=j, fwd_supply=fwd, bwd_supply=bwd, is_pure=is_pure))

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
                    fundamental.append(FundComp(prod_idx=pi, cons_idx=ci, species=s_bit))

    result = CompResult(basic=basic_pairs, pure=pure_pairs, fundamental=fundamental)

    if counters is not None:
        counters.inc("comp.n_basic",       len(basic_pairs))
        counters.inc("comp.n_pure",        len(pure_pairs))
        counters.inc("comp.n_fundamental", len(fundamental))

    return result
