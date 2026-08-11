"""
witness.py — A concrete flux witness for a whole organization: a single
vector v over R_X with v_r > 0 for EVERY triggered reaction and S v >= 0
(Def. 4, Complexity 2022), plus the resulting net production (S v)_s per
species -- directly the quantity Sec. 4.2's "overproducible w.r.t. X" test
compares against zero, made concrete and printable instead of just a
species/F membership bit.

Reuses the SAME minimize_sv LP that circuits.py already trusts for the
per-circuit Theorem 2.16 check (same epsilon/margin, same zero-flux bug
fix -- see the pycot_core_self_maintenance_fix memory), just called once on
the WHOLE domain (X_full, R_X) instead of per-circuit submatrices. This is
sound directly from Theorem 1 (Complexity 2022): X self-maintains iff every
Cj ∪ F ∪ E does, so whenever DecompositionResult.is_organization is True,
minimize_sv on the full domain.S is guaranteed feasible -- this module just
asks for that witness explicitly rather than only trusting the boolean.

CAVEAT on `compute_witness`'s numbers: minimize_sv's LP has a FIXED flux
sum (sum(v) = n_reactions) and a constant objective, so it returns an
ARBITRARY feasible vertex, not a "typical" or maximal operating point --
in practice this usually dumps nearly the whole flux budget onto one
reaction (often the food inflow) and leaves every other triggered reaction
sitting at the epsilon floor (1e-6). That is a perfectly valid Def. 4
witness (every reaction IS strictly positive), but its net-production
numbers read as "almost everything is exactly zero" even for species that
COULD be produced in real surplus. `max_overproduction` below answers the
more useful question directly: for each species in F, what is the LARGEST
net production any feasible flux can achieve? (The same per-species LP
overproduction.py already solves to build F in the first place --
test_overproducible's `value` -- just no longer discarded.)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .bridge import SODomain
from .overproduction import test_overproducible
from pyCOT.analysis.organizations.self_maintenance import minimize_sv

_EPSILON_FLOOR = 1e-4  # display threshold: minimize_sv's epsilon is 1e-6, so
                        # anything below ~1e-4 is that floor, not real signal.


@dataclass
class OrganizationWitness:
    """
    sp_mask            : the organization's X_full species mask.
    reaction_ids       : R_X, same order as `flux`.
    flux               : witness v, v_r > 0 for every r in reaction_ids
                          (see module docstring: an arbitrary feasible
                          vertex, not a maximal/typical one).
    net_production     : species index -> (S v)_species under `flux`.
    max_overproduction : species index -> the LARGEST net production any
                          feasible flux over R_X can achieve for that
                          species, individually maximized (only computed
                          for species in F -- catalysts and circuit
                          members are excluded, their max is 0 by
                          definition/Theorem 2.16). This is the number
                          that answers "how much can this organization
                          actually produce of X", as opposed to
                          `net_production`'s "how much does one arbitrary
                          feasible flux happen to produce".
    """
    sp_mask: int
    reaction_ids: tuple[int, ...]
    flux: tuple[float, ...]
    net_production: dict[int, float]
    max_overproduction: dict[int, float] = field(default_factory=dict)

    def report(self, rn_data) -> str:
        lines = ["  witness flux (reactions with rate > epsilon-floor, sorted by rate desc):"]
        by_rate = sorted(zip(self.reaction_ids, self.flux), key=lambda p: -p[1])
        for r, rate in by_rate:
            note = "  (epsilon floor)" if rate < _EPSILON_FLOOR else ""
            lines.append(f"    {rn_data.reaction_name(r):8s} v={rate:12.6f}{note}")
        lines.append("  net production per species under this witness (S v):")
        for sp, val in sorted(self.net_production.items(), key=lambda p: -p[1]):
            v = 0.0 if abs(val) < _EPSILON_FLOOR else val
            tag = "overproduced" if v > 0 else ("balanced" if v >= 0 else "** NEGATIVE (bug) **")
            lines.append(f"    {rn_data.species_name(sp):12s} {v:12.6f}  {tag}")
        if self.max_overproduction:
            lines.append("  MAXIMUM achievable net production per F species (separately optimized):")
            for sp, val in sorted(self.max_overproduction.items(), key=lambda p: -p[1]):
                lines.append(f"    {rn_data.species_name(sp):12s} {val:12.6f}")
        return "\n".join(lines)


def compute_witness(domain: SODomain, f_mask: int = 0) -> OrganizationWitness | None:
    """None if the domain is not actually self-maintaining (should only
    happen if called on a non-organization; callers should check
    DecompositionResult.is_organization first).

    `f_mask`, if given, additionally reports max_overproduction for every
    species in it (pass DecompositionResult.F_mask)."""
    ok, v = minimize_sv(domain.S)
    if not ok:
        return None
    v = np.asarray(v, dtype=float)
    sv = domain.S @ v
    net = {sp: float(sv[i]) for i, sp in enumerate(domain.sp_indices)}

    max_op: dict[int, float] = {}
    n_r = len(domain.R_X)
    for i, sp in enumerate(domain.sp_indices):
        if not (f_mask >> sp) & 1:
            continue
        _ok, _v, value = test_overproducible(domain.S, i, n_r)
        max_op[sp] = value

    return OrganizationWitness(
        sp_mask=domain.sp_mask,
        reaction_ids=tuple(domain.R_X),
        flux=tuple(float(x) for x in v),
        net_production=net,
        max_overproduction=max_op,
    )
