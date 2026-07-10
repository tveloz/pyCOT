"""
gold_networks.py — Hand-crafted small networks with known ERC structure.

Each network is defined as (supp, prod, expected_ercs).
supp / prod are lists of int bitmasks.  Species indices match the comments.

These are used in test_erc.py for exact correctness checks.
"""
from __future__ import annotations
from dataclasses import dataclass


@dataclass
class GoldNet:
    name: str
    n_species: int
    species: list[str]           # name at index i
    supp: list[int]              # supp bitmask per reaction
    prod: list[int]              # prod bitmask per reaction
    E0_mask: int                 # expected E0 bitmask
    # expected_ercs: list of (erc_mask, n_min_bases, is_persistent)
    expected_ercs: list[tuple[int, int, bool]]
    note: str = ""


# ---------------------------------------------------------------------------
# Gold 1: Simple loop  s0 → s1 → s2 → s0
# ---------------------------------------------------------------------------
#
# s0=bit0, s1=bit1, s2=bit2
# r0: {s0} → {s1}
# r1: {s1} → {s2}
# r2: {s2} → {s0}
#
# closure({s0}) = closure({s1}) = closure({s2}) = {s0, s1, s2}
# One ERC: mask=0b111=7, 3 reactions each firing from their single reactant.
# MinBas = {1, 2, 4} (three singletons, all are minimal).
# req = {}: s0 produced by r2, s1 by r0, s2 by r1 → persistent.
#
GOLD1_LOOP = GoldNet(
    name="loop",
    n_species=3,
    species=["s0", "s1", "s2"],
    supp=[0b001, 0b010, 0b100],  # r0 needs s0, r1 needs s1, r2 needs s2
    prod=[0b010, 0b100, 0b001],  # r0 gives s1, r1 gives s2, r2 gives s0
    E0_mask=0,
    expected_ercs=[
        (0b111, 3, True),   # one ERC, 3 min-bases, persistent
    ],
    note="Minimal persistent network; ERC = full set",
)


# ---------------------------------------------------------------------------
# Gold 2: Three-ERC hierarchy
# ---------------------------------------------------------------------------
#
# s0=bit0, s1=bit1, s2=bit2, s3=bit3
# r0: {s0}     → {s1}       # s0 produces s1
# r1: {s1}     → {s0}       # s1 produces s0  → E1 = {s0, s1}
# r2: {s0, s2} → {s3}       # together produce s3  → E3 = full set
# r3: {s3}     → {s2}       # s3 regenerates s2  → E2 = {s2, s3}
#
# closure({s0}) = closure({s1}) = {s0,s1} = 0b0011 = 3
# closure({s3})                 = {s2,s3} = 0b1100 = 12  ← non-persistent (needs s3)
# closure({s0,s2})              = {s0,s1,s2,s3} = 0b1111 = 15
#
# ERC1 (mask=3):  r0,r1  MinBas=[{s0},{s1}]  req=∅  → persistent
# ERC2 (mask=12): r3     MinBas=[{s3}]        req={s3} → NOT persistent
# ERC3 (mask=15): r2     MinBas=[{s0,s2}]     req=∅  → persistent
#
GOLD2_HIERARCHY = GoldNet(
    name="hierarchy",
    n_species=4,
    species=["s0", "s1", "s2", "s3"],
    supp=[0b0001, 0b0010, 0b0101, 0b1000],
    prod=[0b0010, 0b0001, 0b1000, 0b0100],
    E0_mask=0,
    expected_ercs=[
        (0b0011, 2, True),   # E1 = {s0, s1}: two min-bases, persistent
        (0b1100, 1, False),  # E2 = {s2, s3}: needs external s3, non-persistent
        (0b1111, 1, True),   # E3 = full set: single min-basis {s0,s2}, persistent
    ],
    note="Three-ERC hierarchy; E2={s2,s3} is non-persistent; E3 requires s0 AND s2",
)


# ---------------------------------------------------------------------------
# Gold 3: MinBas example  a→b, a+b→2b
# ---------------------------------------------------------------------------
#
# a=bit0, b=bit1
# r0: {a}    → {b}    # a alone produces b
# r1: {a, b} → {2b}   # a+b together produce b (stoich 2, but bitmask same as b)
#
# closure({a}) = {a, b}: r0 fires (supp={a}⊆{a}) → adds b
# closure({b}) = {b}: no reaction has supp⊆{b} alone (r0 needs a, r1 needs a)
# Actually wait: r0 has supp={a}. closure({a}) = {a,b}. Only one ERC.
# MinBas: supp(r0)={a}=0b01, supp(r1)={a,b}=0b11
# Minimal element: {a} is ⊊ {a,b} → MinBas = [{a}]
#
# ERC: mask=0b11=3, MinBas=[0b01], NOT persistent (a consumed but never produced)
# req = {a} (a consumed by r0 and r1 but no reaction produces a)
#
GOLD3_MINBAS = GoldNet(
    name="minbas",
    n_species=2,
    species=["a", "b"],
    supp=[0b01, 0b11],
    prod=[0b10, 0b10],
    E0_mask=0,
    expected_ercs=[
        (0b11, 1, False),  # one ERC, MinBas=[{a}], NOT persistent (a never produced)
    ],
    note="r1's support {a,b} is dominated by r0's {a} → MinBas has one element; a never produced → not persistent",
)


# ---------------------------------------------------------------------------
# Gold 4: Inflow / E0  ∅→a, a→b
# ---------------------------------------------------------------------------
#
# a=bit0, b=bit1
# r0: {} → {a}       # inflow: a is always available
# r1: {a} → {b}      # a produces b
#
# E0 = closure({a}) = {a, b} (r1 fires immediately)
# supp_q: r0 → 0 (inflow), r1 → {a}&~E0 = 0 (a is in E0)
# All reactions have supp_q=0 → no ERCs.
#
GOLD4_INFLOW = GoldNet(
    name="inflow_E0",
    n_species=2,
    species=["a", "b"],
    supp=[0b00, 0b01],
    prod=[0b01, 0b10],
    E0_mask=0b11,              # both a and b in E0
    expected_ercs=[],          # no non-trivial ERCs
    note="Everything in E0 — no ERCs",
)


# ---------------------------------------------------------------------------
# Gold 5: Non-persistent ERC  (requires external species)
# ---------------------------------------------------------------------------
#
# a=bit0, b=bit1, c=bit2
# r0: {a, b} → {c}   # consumes a and b, produces c
# r1: {c}    → {b}   # c regenerates b
#
# No inflow → E0 = 0.
# closure({a,b}) = {a,b,c}: r0 fires → c added; r1 fires → b added (already there)
# closure({c}) = {b,c}: r1 fires → b added
#
# ERC1 = closure({c}) = {b,c} = 0b110 = 6
#   reactions: [r1]
#   MinBas: [{c}=0b100]
#   req: supp(R_{b,c}) = {c}; prod(R_{b,c}) = {b}; req = {c}&~{b} = {c} ≠ ∅
#   → NOT persistent (needs external c)
#
# ERC2 = closure({a,b}) = {a,b,c} = 0b111 = 7
#   reactions: [r0, r1]
#   MinBas: [{a,b}=0b011]; also {c} gives only {b,c} not {a,b,c}, so it's
#           not a basis for ERC2.  Wait: supp_q for r0={a,b}=0b011, r1={c}=0b100.
#           closure(0b011) = 0b111 ✓, closure(0b100) = 0b110 ≠ 0b111 → r1 belongs to ERC1.
#   So ERC2 contains only r0 (with supp {a,b}), and r0's basis is {a,b}.
#   req: supp=({a,b}) ∪ {} = {a,b}... wait, R_{a,b,c} = {r : supp⊆{a,b,c}} = {r0, r1}
#        prod = {c} ∪ {b} = {b,c}
#        agg_supp = {a,b} ∪ {c} = {a,b,c}
#        req = {a,b,c} & ~{b,c} = {a}  → NOT persistent (needs external a)
#
GOLD5_NONPERSISTENT = GoldNet(
    name="nonpersistent",
    n_species=3,
    species=["a", "b", "c"],
    supp=[0b011, 0b100],
    prod=[0b100, 0b010],
    E0_mask=0,
    expected_ercs=[
        (0b110, 1, False),   # ERC {b,c}: not persistent (needs c)
        (0b111, 1, False),   # ERC {a,b,c}: not persistent (needs a)
    ],
    note="Two non-persistent ERCs; second requires external a",
)


ALL_GOLD = [GOLD1_LOOP, GOLD2_HIERARCHY, GOLD3_MINBAS, GOLD4_INFLOW, GOLD5_NONPERSISTENT]
