"""
raf_algo.py — RAF, maxRAF, irrRAF, and closed RAF (Hordijk & Steel; relational
reading per Veloz, decomposing_RAF_v2.pdf, Sec. 3 / Def 3.1 / Def 3.3).

A non-empty subset R' subset R is a RAF relative to food set F0 iff, writing
W = gen(F0, R') (Sec. 3, crs.gen — a flat set union, NOT an iterative
reachability closure; see that function's docstring for why):

  (a) Food-generation (relational) : for every r in R' and every reactant
      s of r, either s in F0 or s in prod(R') (i.e. s in W).
  (b) Reflexive autocatalysis      : every reaction in R' is catalysed by
      at least one species in W.

This is the paper's deliberately WEAKER "relational" reading, not the
classical Hordijk-Steel constructive one (reactions orderable so each
reactant predates its own production). The relational reading is what
admits fragile circuits (mutually/cyclically self-sustaining species
groups) as valid RAF members -- the whole point of this project's bridge
to Chemical Organization Theory. Do not swap in a constructive/iterative
closure check; that silently reverts to the degenerate classical reading
(Remark 3.2: "a constructively food-generated RAF has an empty fragile
part... precisely the degeneracy we set out to avoid").

RAFs are closed under union, so there is a unique maxRAF, computed by the
standard polynomial-time iterative-removal algorithm: repeatedly strip any
reaction that currently fails (a) or (b) w.r.t. the shrinking gen(F0,*),
until nothing more can be removed (or the set is empty -- no RAF exists).

A RAF R' is CLOSED (Def 3.3) if it contains every reaction of R whose
reactants and at least one catalyst already lie in gen(F0, R') -- i.e. no
"already enabled" reaction has been arbitrarily left out. maxRAF is always
closed (see close_raf's docstring for the argument); arbitrary sub-RAFs
(e.g. an irrRAF found by trimming) generally are not, and must be closed
explicitly before feeding them to the CRS<->COT bridge (Theorem 4.3
requires a CLOSED RAF).

An irrRAF (irreducible RAF) is a RAF with no proper subset that is itself
a RAF. Finding ALL irrRAFs is generally exponential; this module provides
the standard practical substitute: greedy random-order "trimming" from a
starting RAF (default: maxRAF), repeated with different removal orders to
sample the space of irrRAFs. For small systems, exhaustive irrRAF search
is also provided.
"""
from __future__ import annotations

import random

from .crs import CRS, Reaction, gen


def is_raf(crs: CRS, reaction_names) -> bool:
    """
    Direct check: is `reaction_names` (iterable of names) a RAF?

    By convention a RAF is non-empty: the empty set satisfies both
    conditions vacuously (no reactions to violate them), but is not
    considered a legitimate RAF -- otherwise trimming would always be
    able to "succeed" by removing every last reaction.
    """
    names = set(reaction_names)
    if not names:
        return False
    rxns = [crs.reactions[n] for n in names]
    W = gen(crs.food, rxns)
    for r in rxns:
        if not (r.reactants <= W):
            return False
        if not (r.catalysts & W):
            return False
    return True


def compute_maxRAF(crs: CRS) -> frozenset[str]:
    """
    The standard iterative-removal algorithm. Returns the (unique) maxRAF
    as a frozenset of reaction names -- empty if no RAF exists.
    """
    current = set(crs.reactions.keys())
    while True:
        rxns = [crs.reactions[n] for n in current]
        W = gen(crs.food, rxns)
        keep = {
            n for n in current
            if crs.reactions[n].reactants <= W and (crs.reactions[n].catalysts & W)
        }
        if keep == current:
            return frozenset(current)
        if not keep:
            return frozenset()
        current = keep


def close_raf(crs: CRS, reaction_names) -> frozenset[str]:
    """
    Close a RAF per Def 3.3: repeatedly add any reaction of `crs.reactions`
    (not yet included) whose reactants and at least one catalyst already
    lie in gen(F0, current), until no more can be added.

    maxRAF is automatically closed: if some r not in maxRAF had its
    reactants/catalyst already available in gen(F0, maxRAF), then
    maxRAF u {r} would still satisfy both RAF conditions for every member
    (gen() only grows, so no existing member's conditions can break, and r
    itself now satisfies both by hypothesis) -- contradicting maximality.
    So this function matters mainly for closing smaller RAFs (e.g. an
    irrRAF found by trim_to_irr_raf) before feeding them to the CRS<->COT
    bridge, which requires a closed RAF (Theorem 4.3).
    """
    current = set(reaction_names)
    while True:
        W = gen(crs.food, [crs.reactions[n] for n in current])
        addable = {
            name for name, r in crs.reactions.items()
            if name not in current and r.reactants <= W and (r.catalysts & W)
        }
        if not addable:
            return frozenset(current)
        current |= addable


def trim_to_irr_raf(crs: CRS, reaction_subset, *, rng: random.Random | None = None) -> frozenset[str]:
    """
    Greedily remove reactions (in random order) from `reaction_subset` one
    at a time, keeping the removal iff the remainder is still a RAF, until
    no more reactions can be removed. The result is AN irrRAF contained in
    `reaction_subset` (not necessarily unique -- depends on removal order).
    `reaction_subset` must itself already be a RAF.
    """
    rng = rng or random.Random()
    current = set(reaction_subset)
    order = list(current)
    rng.shuffle(order)
    for n in order:
        if n not in current:
            continue
        trial = current - {n}
        if is_raf(crs, trial):
            current = trial
    return frozenset(current)


def sample_irr_rafs(crs: CRS, *, n_samples: int = 50, seed: int | None = None,
                     start: frozenset[str] | None = None) -> list[frozenset[str]]:
    """
    Sample up to `n_samples` irrRAFs by repeated random-order trimming from
    `start` (default: maxRAF). Returns the distinct irrRAFs found (order
    not guaranteed, duplicates removed). Not exhaustive -- see
    enumerate_irr_rafs for small systems.
    """
    base = start if start is not None else compute_maxRAF(crs)
    if not base:
        return []
    rng = random.Random(seed)
    found: set[frozenset[str]] = set()
    for _ in range(n_samples):
        found.add(trim_to_irr_raf(crs, base, rng=rng))
    return sorted(found, key=lambda s: (len(s), sorted(s)))


def enumerate_irr_rafs(crs: CRS, *, start: frozenset[str] | None = None,
                        max_reactions: int = 22) -> list[frozenset[str]]:
    """
    Exhaustive irrRAF enumeration by brute-force subset search over `start`
    (default: maxRAF) -- exponential, only for small systems. Guards
    against runaway cost via `max_reactions` (raises if exceeded).

    A RAF is irreducible iff none of its proper subsets is a RAF; since
    RAF-ness is not monotone in an obvious way, we check all subsets in
    increasing size order and keep the minimal RAFs found (any RAF that is
    a superset of an already-found irrRAF is skipped -- it cannot itself
    be irreducible).
    """
    base = sorted(start if start is not None else compute_maxRAF(crs))
    if not base:
        return []
    if len(base) > max_reactions:
        raise ValueError(
            f"enumerate_irr_rafs: {len(base)} reactions exceeds max_reactions="
            f"{max_reactions} (exponential search) -- use sample_irr_rafs instead."
        )

    from itertools import combinations

    irr: list[frozenset[str]] = []
    n = len(base)
    for size in range(1, n + 1):
        for combo in combinations(base, size):
            s = frozenset(combo)
            if any(prev <= s for prev in irr):
                continue
            if is_raf(crs, s):
                irr.append(s)
    return irr
