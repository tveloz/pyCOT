"""
circuits.py — Fragile circuits: dynamical connection on C = X_full \\ (E ∪ F),
and per-circuit self-maintenance (Theorem 2.16).

Dynamical connection (finer than plain structural/graph connection): two
species of C are directly connected iff some reaction r of R_X has AT
LEAST ONE REACTANT in C, and both species appear in r's support or
products. Requiring a C-reactant is what keeps this "dynamical" rather
than structural — it specifically excludes routing a connection through a
reaction whose only C-touching role is producing a C species from F
reactants, which would wrongly fuse independent circuits that merely share
a common F-supplied input.

Theorem 2.16: X is self-maintaining iff every fragile circuit Di is
self-maintaining w.r.t. its own path R*_i = {r in R_X : support(r) has a
species of Di}, i.e. iff there is v_i >= 0 over R*_i with S_i @ v_i >= 0,
where S_i is the stoichiometry matrix restricted to Di's rows and R*_i's
columns. This is a LOCAL check — it does not care what F/E-supplied inputs
those reactions also draw on, only Di's own net balance.
"""
from __future__ import annotations

import numpy as np

from .bridge import SODomain
from .types import FragileCircuit
from pyCOT.analysis.organizations.self_maintenance import minimize_sv


class _UnionFind:
    def __init__(self, items):
        self.parent = {i: i for i in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def _mask_to_indices(mask: int) -> list[int]:
    out = []
    m = mask
    while m:
        lsb = m & (-m)
        out.append(lsb.bit_length() - 1)
        m &= m - 1
    return out


def dynamical_components(domain: SODomain, rn_data, c_mask: int) -> list[int]:
    """Return the list of connected-component masks partitioning c_mask."""
    c_species = _mask_to_indices(c_mask)
    if not c_species:
        return []
    uf = _UnionFind(c_species)
    for r in domain.R_X:
        supp = rn_data.supp_raw[r]
        if (supp & c_mask) == 0:
            continue  # no reactant in C -> does not connect anything in C
        touched = (supp | rn_data.prod_raw[r]) & c_mask
        idxs = _mask_to_indices(touched)
        for k in range(1, len(idxs)):
            uf.union(idxs[0], idxs[k])

    groups: dict[int, int] = {}
    for sp in c_species:
        root = uf.find(sp)
        groups[root] = groups.get(root, 0) | (1 << sp)
    return list(groups.values())


def _reactions_consuming(rn_data, R_X: list[int], di_mask: int) -> list[int]:
    return sorted(r for r in R_X if (rn_data.supp_raw[r] & di_mask) != 0)


def build_circuit(domain: SODomain, rn_data, di_mask: int) -> FragileCircuit:
    r_star = _reactions_consuming(rn_data, domain.R_X, di_mask)
    di_species = _mask_to_indices(di_mask)

    if not r_star:
        # No reaction consumes any Di species -> cannot self-maintain
        # unless Di is empty (never happens: di_mask is non-empty by
        # construction) or every Di species also never needs production
        # (impossible: a non-E, non-F species is by definition not
        # overproducible, i.e. it has no net-positive producer here).
        return FragileCircuit(species_mask=di_mask, reaction_ids=(), is_self_maintaining=False)

    rows = [domain.row_of(sp) for sp in di_species]
    cols = [domain.col_of(r) for r in r_star]
    S_i = domain.S[np.ix_(rows, cols)]

    ok, v = minimize_sv(S_i)
    flux = tuple(float(x) for x in v) if ok else None
    return FragileCircuit(species_mask=di_mask, reaction_ids=tuple(r_star),
                           is_self_maintaining=bool(ok), flux=flux)


def compute_fragile_circuits(domain: SODomain, rn_data, e_mask: int, f_mask: int) -> list[FragileCircuit]:
    c_mask = domain.X_full_mask & ~(e_mask | f_mask)
    components = dynamical_components(domain, rn_data, c_mask)
    return [build_circuit(domain, rn_data, comp) for comp in components]
