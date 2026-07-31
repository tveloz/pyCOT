"""
validate_raf_algo.py — Correctness checks for maxRAF / irrRAF / close_raf.

Run directly: python projects/RAF_Comparison/tests/validate_raf_algo.py

Two tiers of oracle:
  1. The paper's own worked examples (decomposing_RAF_v2.pdf, Examples 4.5,
     4.6, 5.8) -- hand-derived by the source itself, the strongest oracle
     available. Network N realizes; its perturbation N' does not.
  2. Hand-derived toy examples exercising edge cases: no-RAF-exists,
     food-catalyzed trivial RAF, self/bootstrap-catalyzed RAF, multi-round
     iterative pruning, maxRAF decomposing into two irrRAFs, and (new) a
     genuinely MUTUAL two-species mini-cycle that the relational reading
     must accept and the old (buggy) constructive-closure reading would
     have rejected -- a direct regression test for the Def-3.1 fix.
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
if _proj not in sys.path:
    sys.path.insert(0, _proj)

from raf.crs import CRS, gen
from raf.raf_algo import compute_maxRAF, is_raf, sample_irr_rafs, enumerate_irr_rafs, close_raf

failures = []


def check(name, cond, detail=""):
    status = "OK" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail and not cond else ""))
    if not cond:
        failures.append(name)


# ===========================================================================
# Tier 1: the paper's own Network N / N' (Sec 4.3, Examples 4.5-4.6)
# ===========================================================================

print("Network N (Example 4.5) — closed RAF, realization succeeds")
# Note: inflow of food species is NOT a reaction of the CRS itself -- Sec 3
# models food purely via F0 (a plain species set); "inflow reactions"
# (∅->s) only appear later, synthesized by the COT translation (Def 4.1(i)).
# So r0 from the paper's Eq. (8)-style listing is deliberately omitted here;
# R = {r1,r2,r3,r4} only.
network_N = CRS.build(
    reactions=[
        ("r1", ["f", "a"], ["b", "g"], ["b"]),
        ("r2", ["b"], ["a"], ["a"]),
        ("r3", ["g", "c"], ["d"], ["d"]),
        ("r4", ["d"], ["c", "h"], ["c"]),
    ],
    food=["f"],
)
maxraf_N = compute_maxRAF(network_N)
check("maxRAF(N) == all 4 reactions", maxraf_N == frozenset({"r1", "r2", "r3", "r4"}),
      f"got {maxraf_N}")
gen_N = gen(network_N.food, [network_N.reactions[n] for n in maxraf_N])
check("gen(F0, maxRAF) == all of M", gen_N == network_N.species, f"got {gen_N}")
check("maxRAF(N) is closed (== close_raf(maxRAF))",
      close_raf(network_N, maxraf_N) == maxraf_N)

print("\nNetwork N' (Example 4.6) — r1 perturbed to f+2a->b+g, closed RAF is unchanged")
network_Np = CRS.build(
    reactions=[
        ("r1p", ["f", "a"], ["b", "g"], ["b"], {"a": 2}, {}),   # f + 2a -> b + g
        ("r2", ["b"], ["a"], ["a"]),
        ("r3", ["g", "c"], ["d"], ["d"]),
        ("r4", ["d"], ["c", "h"], ["c"]),
    ],
    food=["f"],
)
maxraf_Np = compute_maxRAF(network_Np)
check("maxRAF(N') is still the full 4-reaction set (RAF-ness is stoichiometry-blind)",
      maxraf_Np == frozenset({"r1p", "r2", "r3", "r4"}), f"got {maxraf_Np}")
# NOTE: RAF conditions (a)/(b) never look at stoichiometric coefficients, only
# at reactant/product SET membership -- {f,a} vs {f,a,a} both just mean
# "reactants are {f,a}" under a set-based Reaction.reactants. So N' is a
# closed RAF and a semi-organization (Theorem 4.3) exactly like N; the
# failure of self-maintenance in N' (Example 4.6) is a purely STOICHIOMETRIC
# fact (the flux cone is empty), invisible at the RAF/relational level. This
# is exactly the semi-organization/organization gap the paper is about --
# verified properly once decomp/ is wired in (see validate_bridge.py).

print("\nExample 5.8 — dependency structure of Network N's maxRAF")
# D1={a,b} depth 1 (interface {f} subset F0), D2={c,d} depth 2 (interface {g}).
# Deferred to validate_dependency.py once decomp/ circuits are computed on
# the translated network -- this file only checks the RAF-algebra layer.

# ===========================================================================
# Tier 2: hand-derived toy examples
# ===========================================================================

print("\nExample 1 — no RAF exists (dangling catalyst never produced)")
crs1 = CRS.build(
    reactions=[("r1", ["a"], ["b"], ["c"])],   # catalyst c never produced, not food
    food=["a"],
)
m1 = compute_maxRAF(crs1)
check("maxRAF is empty", m1 == frozenset(), f"got {m1}")

print("\nExample 2 — trivial single-reaction RAF, catalyzed by food")
crs2 = CRS.build(
    reactions=[("r1", ["a", "b"], ["c"], ["a"])],   # catalyzed by food species a
    food=["a", "b"],
)
m2 = compute_maxRAF(crs2)
check("maxRAF == {r1}", m2 == frozenset({"r1"}), f"got {m2}")
check("is_raf({r1}) is True", is_raf(crs2, {"r1"}))

print("\nExample 3 — self-catalytic (bootstrap) single-reaction RAF")
crs3 = CRS.build(
    reactions=[("r1", ["a", "b"], ["c"], ["c"])],   # catalyzed by its own product
    food=["a", "b"],
)
m3 = compute_maxRAF(crs3)
check("maxRAF == {r1} (bootstrap catalysis via own product)", m3 == frozenset({"r1"}), f"got {m3}")

print("\nExample 4 — requires 2 rounds of iterative pruning")
crs4 = CRS.build(
    reactions=[
        ("r1", ["a"], ["b"], ["d"]),   # needs d, only produced by r2
        ("r2", ["b"], ["d"], ["b"]),   # self-catalyzed by b
        ("r3", ["d"], ["e"], ["z"]),   # z never produced -> must be pruned
    ],
    food=["a"],
)
m4 = compute_maxRAF(crs4)
check("maxRAF == {r1, r2} (r3 pruned)", m4 == frozenset({"r1", "r2"}), f"got {m4}")

print("\nExample 5 — maxRAF decomposes into two independent irrRAFs")
crs5 = CRS.build(
    reactions=[
        ("r1", ["a", "b"], ["c"], ["c"]),   # independent bootstrap RAF #1
        ("r2", ["a", "b"], ["d"], ["d"]),   # independent bootstrap RAF #2
    ],
    food=["a", "b"],
)
m5 = compute_maxRAF(crs5)
check("maxRAF == {r1, r2}", m5 == frozenset({"r1", "r2"}), f"got {m5}")
irr5_exhaustive = enumerate_irr_rafs(crs5)
check("exhaustive irrRAFs == {{r1}, {r2}}",
      sorted(irr5_exhaustive) == sorted([frozenset({"r1"}), frozenset({"r2"})]),
      f"got {irr5_exhaustive}")
irr5_sampled = sample_irr_rafs(crs5, n_samples=50, seed=0)
check("sampled irrRAFs also find both {r1} and {r2}",
      set(irr5_sampled) == {frozenset({"r1"}), frozenset({"r2"})},
      f"got {irr5_sampled}")

print("\nExample 6 (NEW, regression test for the Def-3.1 relational-reading fix)")
print("  a genuinely mutual 2-cycle: r1: x->y (cat y), r2: y->x (cat x), no food link at all")
crs6 = CRS.build(
    reactions=[
        ("r1", ["x"], ["y"], ["y"]),
        ("r2", ["y"], ["x"], ["x"]),
    ],
    food=[],
)
# Under the OLD (buggy) constructive-closure reading, cl(F0={}, {r1,r2})
# would stay empty forever (neither reaction's reactant is ever "reached"
# from an empty seed), so is_raf would wrongly return False here. Under the
# correct relational reading, gen({}, {r1,r2}) = {x,y} (flat union of
# products) directly, and both conditions hold.
check("is_raf({r1,r2}) is True under relational reading (mutual bootstrap cycle)",
      is_raf(crs6, {"r1", "r2"}))
m6 = compute_maxRAF(crs6)
check("maxRAF == {r1, r2}", m6 == frozenset({"r1", "r2"}), f"got {m6}")
check("neither {r1} nor {r2} alone is a RAF (each needs the other's product)",
      not is_raf(crs6, {"r1"}) and not is_raf(crs6, {"r2"}))

print("\n" + ("ALL CHECKS PASSED" if not failures else f"FAILURES: {failures}"))
sys.exit(0 if not failures else 1)
