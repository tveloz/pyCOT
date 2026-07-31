"""
validate_bridge.py — Correctness checks for the outflow-free CRS<->COT
bridge (Theorem 4.3) and the stoichiometric realization / decomposition
(Examples 4.5, 4.6, 5.8 of decomposing_RAF_v2.pdf).

Run directly: python projects/RAF_Comparison/tests/validate_bridge.py

Oracle: the paper's own fully-worked Network N (realizes) and its
perturbation N' (does not realize) -- exact expected E/F/circuits/
is_organization values are stated in the source text itself, making this
about as strong an oracle as is available short of the author's own code.
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
if _proj not in sys.path:
    sys.path.insert(0, _proj)

from raf.crs import CRS, gen
from raf.raf_algo import compute_maxRAF
from raf.cot_bridge import cot_translate, is_semi_organization, is_closed, is_semi_self_maintaining
from raf.decomp_shim import decompose_species_set

failures = []


def check(name, cond, detail=""):
    status = "OK" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail and not cond else ""))
    if not cond:
        failures.append(name)


def make_network(coeff_a_in_r1=1):
    rc1 = {"a": coeff_a_in_r1} if coeff_a_in_r1 != 1 else {}
    return CRS.build(
        reactions=[
            ("r1", ["f", "a"], ["b", "g"], ["b"], rc1, {}),
            ("r2", ["b"], ["a"], ["a"]),
            ("r3", ["g", "c"], ["d"], ["d"]),
            ("r4", ["d"], ["c", "h"], ["c"]),
        ],
        food=["f"],
    )


print("Network N (Example 4.5) — Theorem 4.3 + full decomposition")
network_N = make_network()
maxraf_N = compute_maxRAF(network_N)
X_N = gen(network_N.food, [network_N.reactions[n] for n in maxraf_N])
net_N = cot_translate(network_N)

check("X = gen(F0,maxRAF) is closed in the translation", is_closed(net_N, X_N))
check("X is semi-self-maintaining in the translation", is_semi_self_maintaining(net_N, X_N))
check("Theorem 4.3: X is a semi-organization", is_semi_organization(net_N, X_N))

result_N = decompose_species_set(net_N, X_N)
check("E is empty (no COT-sense catalysts: RAF catalysts b,a,d,c are real "
      "reactants/products elsewhere, not net-zero everywhere)",
      result_N.E_mask == 0, f"E={result_N.E_mask:#b}")
import raf.decomp_shim as _shim_mod
shim_N, _ = _shim_mod.build_shim(net_N)
F_names = set(shim_N.bitset_to_names(result_N.F_mask))
check("F == {f, g, h} (native food f, contextually overproduced g,h)",
      F_names == {"f", "g", "h"}, f"got {F_names}")
circuit_species = sorted(set(shim_N.bitset_to_names(c.species_mask)) for c in result_N.circuits)
expected_circuits = sorted([{"a", "b"}, {"c", "d"}])
check("Circuits == {{a,b}, {c,d}} (matches D1, D2 of Example 5.8)",
      circuit_species == expected_circuits, f"got {circuit_species}")
check("X is a full organization (Example 4.5: realization succeeds)",
      result_N.is_organization)
check("Both circuits individually self-maintain",
      all(c.is_self_maintaining for c in result_N.circuits))

print("\nNetwork N' (Example 4.6) — r1: f+2a->b+g, realization fails")
network_Np = make_network(coeff_a_in_r1=2)
maxraf_Np = compute_maxRAF(network_Np)
X_Np = gen(network_Np.food, [network_Np.reactions[n] for n in maxraf_Np])
net_Np = cot_translate(network_Np)

check("Theorem 4.3 still holds: X' is a semi-organization (RAF-ness is "
      "stoichiometry-blind, so this doesn't change)",
      is_semi_organization(net_Np, X_Np))

result_Np = decompose_species_set(net_Np, X_Np)
check("X' is NOT a full organization (Example 4.6: flux cone is empty)",
      not result_Np.is_organization)
shim_Np, _ = _shim_mod.build_shim(net_Np)
failing = [c for c in result_Np.circuits if not c.is_self_maintaining]
failing_species = sorted(set(shim_Np.bitset_to_names(c.species_mask)) for c in failing)
# NOTE on what this does NOT check: Def 2.6/Prop 2.7's F(X) is only DEFINED
# under the standing assumption that X is already an organization (Sec 2.4
# preamble: "Throughout this subsection X is a fixed organization"). For a
# broken semi-org like X', F/circuits genuinely requires the SEPARATE Sec
# 4.4 apparatus (organizational core O(X), leaking vs. dependent circuits,
# a Farkas-certificate-driven construction) to reproduce the paper's exact
# fine-grained claim ("{a,b} is the leaking circuit, {c,d} is a dependent
# circuit that collapses because {a,b} can't supply g"). That apparatus is
# NOT implemented here (deferred, see project memory) -- what IS verified
# below is Theorem 2.16's actual headline claim (X self-maintaining iff
# every circuit self-maintains), which decomp's engine answers directly by
# running the overproducibility LP over the WHOLE of X' at once: since v1
# (r1') is forced to 0 by {a,b}'s own infeasibility (v2>=2v1 and v1>=v2
# force v1=v2=0), g's only producer is silenced, so g (and everything
# downstream: c,d,h) is ALSO not overproducible -- correctly collapsing
# F down to just {f} and merging the rest into one large non-self-
# maintaining circuit. This is a coarser but self-consistent superset
# answer to "is X self-maintaining", not a bug.
check("F collapses to just {f} (g,h are no longer overproducible once "
      "{a,b}'s own flux is forced to zero)",
      set(shim_Np.bitset_to_names(result_Np.F_mask)) == {"f"},
      f"got {set(shim_Np.bitset_to_names(result_Np.F_mask))}")
check("The remaining species merge into one non-self-maintaining circuit "
      "(a,b,c,d,g,h all downstream of the same broken flux)",
      len(failing) == 1 and failing_species[0] == {"a", "b", "c", "d", "g", "h"},
      f"failing circuits: {failing_species}")

print("\n" + ("ALL CHECKS PASSED" if not failures else f"FAILURES: {failures}"))
sys.exit(0 if not failures else 1)
