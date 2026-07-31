"""
validate_dependency.py — Correctness checks for the dependency DAG and the
irrRAF / indecomposability correspondence (Theorems 5.4, 5.5, 5.6),
against decomposing_RAF_v2.pdf's own worked Example 5.8 (Network N).

Run directly: python projects/RAF_Comparison/tests/validate_dependency.py

Expected (verbatim from the paper): D1={a,b} depth 1, interface {f} subset
F0; D2={c,d} depth 2, interface {g}; single edge D1 >- D2. Unique irrRAF is
{r1,r2} (D1's own path). The maxRAF {r1,r2,r3,r4} is indecomposable (its
circuits form a single chain) yet reducible (it properly contains the
irrRAF {r1,r2}).
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
if _proj not in sys.path:
    sys.path.insert(0, _proj)

from raf.crs import CRS, gen
from raf.raf_algo import compute_maxRAF, is_raf, enumerate_irr_rafs
from raf.cot_bridge import cot_translate
from raf.decomp_shim import decompose_species_set, build_shim
from raf.dependency import build_dependency_dag, irrRAF_candidate_names

failures = []


def check(name, cond, detail=""):
    status = "OK" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail and not cond else ""))
    if not cond:
        failures.append(name)


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
X_N = gen(network_N.food, [network_N.reactions[n] for n in maxraf_N])
net_N = cot_translate(network_N)
result_N = decompose_species_set(net_N, X_N)
shim_N, S_full_N = build_shim(net_N)

F0_mask = 0
for s in network_N.food:
    F0_mask |= 1 << shim_N.species_idx(s)

dag = build_dependency_dag(result_N, net_N, shim_N, S_full_N, F0_mask)
check("DAG is activatable (every circuit eventually reachable from food)",
      dag.is_activatable(), f"activated={dag.activated}")

# Identify D1={a,b} and D2={c,d} by species regardless of list order.
def species_of(i):
    return set(shim_N.bitset_to_names(result_N.circuits[i].species_mask))

idx_ab = next(i for i in range(len(result_N.circuits)) if species_of(i) == {"a", "b"})
idx_cd = next(i for i in range(len(result_N.circuits)) if species_of(i) == {"c", "d"})

check("D1={a,b} has depth 1", dag.depth[idx_ab] == 1, f"got {dag.depth[idx_ab]}")
check("D2={c,d} has depth 2", dag.depth[idx_cd] == 2, f"got {dag.depth[idx_cd]}")

interface_ab = set(shim_N.bitset_to_names(dag.interfaces[idx_ab]))
interface_cd = set(shim_N.bitset_to_names(dag.interfaces[idx_cd]))
check("I(D1) == {f} (subset of F0)", interface_ab == {"f"}, f"got {interface_ab}")
check("I(D2) == {g}", interface_cd == {"g"}, f"got {interface_cd}")

check("Single edge D1 >- D2", dag.edges == [(idx_ab, idx_cd)], f"got {dag.edges}")

# Theorem 5.4: irrRAFs are exactly the minimal (depth-1) circuits' own paths.
minimal = dag.minimal_circuits()
check("Exactly one minimal (depth-1) circuit: D1", minimal == [idx_ab], f"got {minimal}")

candidate = irrRAF_candidate_names(result_N, shim_N, idx_ab)
check("D1's path reads back as {r1, r2}", candidate == frozenset({"r1", "r2"}), f"got {candidate}")
check("That candidate is indeed a RAF of the original CRS", is_raf(network_N, candidate))

irr_exhaustive = enumerate_irr_rafs(network_N, start=maxraf_N)
check("It is THE unique irrRAF of the network (matches independent brute-force enumeration)",
      irr_exhaustive == [frozenset({"r1", "r2"})], f"got {irr_exhaustive}")

# Theorem 5.5/5.6: the full maxRAF's circuit set {D1,D2} is indecomposable
# (single maximal element D2, since D1 >- D2) yet reducible (properly
# contains the irrRAF {r1,r2}).
check("{D1,D2} is indecomposable (single maximal element = D2)",
      dag.is_indecomposable({idx_ab, idx_cd}) and dag.maximal_elements({idx_ab, idx_cd}) == [idx_cd])
check("maxRAF is reducible: it properly contains the irrRAF {r1,r2}",
      frozenset({"r1", "r2"}) < maxraf_N)

print("\n" + ("ALL CHECKS PASSED" if not failures else f"FAILURES: {failures}"))
sys.exit(0 if not failures else 1)
