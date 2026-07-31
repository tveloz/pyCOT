"""
cot_bridge.py — Outflow-free CRS -> COT translation (Def 4.1) and the
relational bridge (Theorem 4.3).

Translation (Def 4.1): given a CRS Q = (M, R, C, F0), build Qtilde = (M, Rtilde):
  (i)  an inflow reaction (empty) -> s for every s in F0;
  (ii) for each r in R with reactants A, products B, and EACH catalyst
       c in C(r), the reaction A + c -> B + c (catalyst on both sides, net
       zero effect). The paper writes "a chosen catalyst c in C(r)"
       (singular); we instead emit one translated reaction PER catalyst in
       C(r), which is strictly more general (any single available catalyst
       should be able to trigger r) and reduces to the paper's own
       construction exactly when |C(r)| == 1 -- true in every worked
       example in the source paper.
No outflow reactions are added (this is the whole point: see Remark 4.2 --
adding outflow degenerates the decomposition to X = F with no fragile
circuits, which is exactly the blind spot this paper is designed to avoid).

Theorem 4.3 (relational bridge): for a CLOSED RAF R', X = gen(F0, R') is a
semi-organization of Qtilde (closed + semi-self-maintaining, Def 2.1-2.3).
is_semi_organization() below checks this DIRECTLY (not by trusting the
theorem), so every use of the bridge in this project is empirically
verified, not assumed.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .crs import CRS, Reaction


@dataclass(frozen=True)
class TReaction:
    """One translated-network reaction: signed net coefficients per species."""
    name: str
    reactant_coeffs: dict[str, float]
    product_coeffs: dict[str, float]

    def support(self) -> frozenset[str]:
        return frozenset(s for s, c in self.reactant_coeffs.items() if c > 0)

    def products(self) -> frozenset[str]:
        return frozenset(s for s, c in self.product_coeffs.items() if c > 0)

    def net(self, species: str) -> float:
        return self.product_coeffs.get(species, 0.0) - self.reactant_coeffs.get(species, 0.0)


@dataclass
class TranslatedNet:
    species: frozenset[str]
    reactions: list[TReaction]


def cot_translate(crs: CRS) -> TranslatedNet:
    """
    Def 4.1(ii) reads "a chosen catalyst c in C(r)" -- this presupposes
    C(r) is non-empty. A reaction with NO catalyst at all has no valid
    translated form and is simply ABSENT from Qtilde (it can never satisfy
    RAF condition (b) either, so it is correctly invisible to both sides).
    This matters for Theorem 4.3's proof, which explicitly characterizes
    every triggered non-inflow reaction of Qtilde as "the translation
    A+c->B+c of some r in R" -- there is no third, bare-reaction case.
    Silently including catalyst-less reactions as bare translations (an
    earlier version of this function did) breaks closure: such a reaction
    can be triggered by X = gen(F0,R') without belonging to any RAF, and
    can produce species outside X, so X is no longer closed in Qtilde --
    empirically caught on BIOMD0000000407 (binding/complex-formation
    reactions like Casp3+XIAP -> XIAP_Casp3 have no induced catalyst and
    were silently violating Theorem 4.3 until this fix).
    """
    reactions: list[TReaction] = []
    for s in sorted(crs.food):
        reactions.append(TReaction(f"inflow_{s}", {}, {s: 1.0}))

    for r in crs.reactions.values():
        if not r.catalysts:
            continue
        rc = {s: r.reactant_coeffs.get(s, 1.0) for s in r.reactants}
        pc = {s: r.product_coeffs.get(s, 1.0) for s in r.products}
        catalysts = sorted(r.catalysts)
        for c in catalysts:
            rc_c = dict(rc)
            pc_c = dict(pc)
            rc_c[c] = rc_c.get(c, 0.0) + 1.0
            pc_c[c] = pc_c.get(c, 0.0) + 1.0
            suffix = f"__cat_{c}" if len(catalysts) > 1 else ""
            reactions.append(TReaction(f"{r.name}{suffix}", rc_c, pc_c))

    return TranslatedNet(species=crs.species, reactions=reactions)


# ---------------------------------------------------------------------------
# Relational layer (Def 2.1-2.3) — pure set logic, no LP
# ---------------------------------------------------------------------------

def triggered(net: TranslatedNet, X) -> list[TReaction]:
    Xs = set(X)
    return [r for r in net.reactions if r.support() <= Xs]

def prod_of(reactions) -> frozenset[str]:
    out = set()
    for r in reactions:
        out |= r.products()
    return frozenset(out)

def supp_of(reactions) -> frozenset[str]:
    out = set()
    for r in reactions:
        out |= r.support()
    return frozenset(out)

def is_closed(net: TranslatedNet, X) -> bool:
    RX = triggered(net, X)
    return prod_of(RX) <= set(X)

def req(net: TranslatedNet, X) -> frozenset[str]:
    RX = triggered(net, X)
    return frozenset(supp_of(RX) - prod_of(RX))

def is_semi_self_maintaining(net: TranslatedNet, X) -> bool:
    return len(req(net, X)) == 0

def is_semi_organization(net: TranslatedNet, X) -> bool:
    return is_closed(net, X) and is_semi_self_maintaining(net, X)


def closure(net: TranslatedNet, X) -> frozenset[str]:
    """Smallest closed superset of X (Def 2.1), by iterated product addition."""
    reached = set(X)
    changed = True
    while changed:
        changed = False
        RX = triggered(net, reached)
        p = prod_of(RX)
        if not (p <= reached):
            reached |= p
            changed = True
    return frozenset(reached)
