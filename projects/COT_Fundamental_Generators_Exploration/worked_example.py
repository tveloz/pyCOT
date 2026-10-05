"""
worked_example.py — the companion-paper "Fully Worked Network".

Shared by tests/test_so_search.py (correctness + the vertical-lift regression)
and conjectures/benchmark_suite.py (Conjecture 1 vs Conjecture 2 comparison):
8 reactions over 13 species, with a sterile ERC (E_dagger, requires a
never-produced species) and a higher-order SO (M3) that is reachable ONLY via
a hierarchy vertical lift — not via synergy or minimal-producer
complementarity. This makes it the single clearest concrete illustration of
what the paper's two §6.4 conjectures disagree about.
"""
from __future__ import annotations

from pyCOT.analysis.organizations.cot_types import RNData

WORKED_NAMES = ["a", "b", "c", "d", "e", "f", "g", "h", "k", "p", "q", "w", "z"]
_WORKED_IDX = {n: i for i, n in enumerate(WORKED_NAMES)}


def wmask(*species: str) -> int:
    m = 0
    for s in species:
        m |= 1 << _WORKED_IDX[s]
    return m


def worked_rndata() -> RNData:
    """
    r_A: f+a->2a+p   r_B: p+b->2b+f   r_C: p+c->2c+f
    r_D: g+d->2d+q   r_E: q+e->2e+g   r_F: p+h->2h
    r_H: p+z->z+b    r_G: k+w->2k+p
    """
    reactions = [
        ("r_A", wmask("f", "a"), wmask("a", "p")),
        ("r_B", wmask("p", "b"), wmask("b", "f")),
        ("r_C", wmask("p", "c"), wmask("c", "f")),
        ("r_D", wmask("g", "d"), wmask("d", "q")),
        ("r_E", wmask("q", "e"), wmask("e", "g")),
        ("r_F", wmask("p", "h"), wmask("h")),
        ("r_H", wmask("p", "z"), wmask("z", "b")),
        ("r_G", wmask("k", "w"), wmask("k", "p")),
    ]
    supp = [r[1] for r in reactions]
    prod = [r[2] for r in reactions]
    names = [r[0] for r in reactions]
    n_sp = len(WORKED_NAMES)
    inv = tuple(
        tuple(r for r, s in enumerate(supp) if (s >> i) & 1)
        for i in range(n_sp)
    )
    return RNData(
        n_species=n_sp,
        species_names=tuple(WORKED_NAMES),
        species_index=tuple((n, i) for i, n in enumerate(WORKED_NAMES)),
        n_reactions=len(reactions),
        reaction_names=tuple(names),
        supp_raw=tuple(supp), prod_raw=tuple(prod),
        E0_mask=0, supp_q=tuple(supp), prod_q=tuple(prod),
        species_to_reactions=inv,
    )


def names_of(rn: RNData, mask: int) -> list[str]:
    return sorted(rn.species_name(i) for i in range(rn.n_species) if (mask >> i) & 1)
