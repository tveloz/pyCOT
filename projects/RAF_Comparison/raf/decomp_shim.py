"""
decomp_shim.py — Bridge a TranslatedNet (cot_bridge.py) into the exact
minimal interface Decomposition_Theorem's decomp package needs, so that
project's already-validated E/F/fragile-circuit/self-maintenance engine
(catalysts.py, overproduction.py, circuits.py, core.py) can be reused
DIRECTLY on RAF-derived networks -- no reimplementation, no drift between
the two projects' notions of "organization".

Decomposition_Theorem's decompose(sp_mask, rn_data, S_full) only reads
rn_data.E0_mask and rn_data.supp_raw/prod_raw (see that project's
bridge.py/circuits.py) plus a signed (species x reactions) matrix aligned
to the same compact index -- exactly what we build here from a
TranslatedNet, with a stable name<->index mapping for reporting.
"""
from __future__ import annotations

import os
import sys

import numpy as np

from .cot_bridge import TranslatedNet

_here = os.path.dirname(os.path.abspath(__file__))
_decomp_proj = os.path.normpath(os.path.join(_here, "..", "..", "Decomposition_Theorem"))
if _decomp_proj not in sys.path:
    sys.path.insert(0, _decomp_proj)

from decomp.core import decompose            # noqa: E402
from decomp.types import DecompositionResult  # noqa: E402


class RNShim:
    """Minimal duck-typed stand-in for cot_gen's RNData, built from a TranslatedNet."""

    def __init__(self, species_names: list[str], supp_raw: list[int], prod_raw: list[int],
                 reaction_names: list[str], E0_mask: int = 0):
        self.species_names = tuple(species_names)
        self.species_index = tuple((n, i) for i, n in enumerate(species_names))
        self.n_species = len(species_names)
        self.supp_raw = tuple(supp_raw)
        self.prod_raw = tuple(prod_raw)
        self.reaction_names = tuple(reaction_names)
        self.n_reactions = len(reaction_names)
        self.E0_mask = E0_mask

    def species_name(self, i: int) -> str:
        return self.species_names[i]

    def species_idx(self, name: str) -> int:
        return dict(self.species_index)[name]

    def reaction_name(self, r: int) -> str:
        return self.reaction_names[r]

    def bitset_to_names(self, mask: int) -> list[str]:
        out = []
        m = mask
        i = 0
        while m:
            if m & 1:
                out.append(self.species_names[i])
            m >>= 1
            i += 1
        return sorted(out)

    def names_to_bitset(self, names) -> int:
        idx = dict(self.species_index)
        v = 0
        for name in names:
            v |= 1 << idx[name]
        return v


def build_shim(net: TranslatedNet) -> tuple[RNShim, np.ndarray]:
    """
    Build (RNShim, S_full) for `net`. Species indexed by sorted name
    (deterministic); reactions indexed in `net.reactions` order.
    """
    species_names = sorted(net.species)
    idx = {s: i for i, s in enumerate(species_names)}
    n_sp, n_rx = len(species_names), len(net.reactions)

    supp_raw = [0] * n_rx
    prod_raw = [0] * n_rx
    reaction_names = []
    S = np.zeros((n_sp, n_rx))
    for j, r in enumerate(net.reactions):
        reaction_names.append(r.name)
        for s in r.support():
            supp_raw[j] |= 1 << idx[s]
        for s in r.products():
            prod_raw[j] |= 1 << idx[s]
        for s, c in r.reactant_coeffs.items():
            S[idx[s], j] -= c
        for s, c in r.product_coeffs.items():
            S[idx[s], j] += c

    shim = RNShim(species_names, supp_raw, prod_raw, reaction_names, E0_mask=0)
    return shim, S


def decompose_species_set(net: TranslatedNet, X) -> DecompositionResult:
    """
    Decompose species set X (an iterable of species names) within the
    translated network `net`, reusing Decomposition_Theorem's engine.

    Note E0_mask=0 in the shim: unlike cot_gen's own quotiented ERC
    space, X here is expressed directly over `net`'s real species names
    (RAF's F0 is already folded into X via gen(F0,R'), and the shim's own
    inflow reactions are ordinary reactions of `net`) -- there is no
    separate "food is implicit and excluded from X" convention to
    reproduce here, so the bridge.so_domain() call's X_full = X | E0_mask
    reduces to plain X.
    """
    shim, S_full = build_shim(net)
    sp_mask = 0
    for s in X:
        sp_mask |= 1 << shim.species_idx(s)
    return decompose(sp_mask, shim, S_full)
