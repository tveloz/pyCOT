"""
bridge.py — Species-index / stoichiometric-matrix bridge between cot_gen's
bitmask-indexed ERC world and pyCOT's real ReactionNetwork.

Why this exists
----------------
cot_gen's ERC species masks are built from the QUOTIENTED support/product
masks (`RNData.supp_q`/`prod_q`, which zero out `E0_mask` — the native-food
closure, i.e. species reachable "for free" from inflow reactions). This
means no EPM/ESPM species set produced by cot_gen ever contains an E0
species: quotienting removes them from prod_q entirely, so closure_opt can
never add them to species_mask (see cot_gen/erc.py, cot_gen/io_pyCOT.py).

The Decomposition Theorem needs the REAL picture: native food species are
part of the network's actual stoichiometry and are trivially overproducible
(Prop 5.1a of the RAF-COT paper: native food F0 ⊆ F). So every SO handed to
this module is expanded to X_full = X ∪ E0 before any LP/stoichiometric
work — otherwise R_X would be missing reactions that touch food species and
the whole E/F/circuit classification would be wrong.

Index alignment (verified, not assumed)
----------------------------------------
`io_pyCOT._species_list`/`_build_masks` build cot_gen's compact index via
`sorted(rn.species(), key=lambda s: s.index)` / `sorted(rn.reactions(),
key=lambda r: r.node.index)`. `ReactionNetwork._build_matrix` (which backs
`stoichiometry_matrix()`) uses the IDENTICAL sort keys. Consequently
cot_gen species/reaction index i is, row/column-for-row/column, the same
index into `rn.stoichiometry_matrix()` — no remapping needed, PROVIDED `rn`
is the same object (or an unmutated copy) used to build `rn_data`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def build_full_stoich(rn) -> np.ndarray:
    """
    Build the full (n_species x n_reactions) signed stoichiometry matrix,
    row/column order aligned with RNData's compact index (see module
    docstring). Call once per network and reuse — this is the expensive
    step, everything downstream just slices it.
    """
    return np.asarray(rn.stoichiometry_matrix(), dtype=float)


@dataclass
class SODomain:
    """
    The real, LP-ready view of one semi-organization X.

    sp_mask       : the ORIGINAL (quotiented) species mask handed in.
    X_full_mask   : sp_mask | E0_mask — what actually gets analyzed.
    sp_indices    : sorted species indices of X_full_mask (row order of S).
    R_X           : sorted reaction indices with support ⊆ X_full_mask
                    (column order of S). Includes inflow reactions
                    (supp_raw == 0), which witness E0's overproducibility.
    S             : (len(sp_indices) x len(R_X)) submatrix of the full
                    stoichiometry matrix, restricted and reordered to
                    (sp_indices, R_X).
    """
    sp_mask: int
    X_full_mask: int
    sp_indices: list[int]
    R_X: list[int]
    S: np.ndarray

    def row_of(self, sp_idx: int) -> int:
        return self._row.get(sp_idx, -1)

    def col_of(self, r_idx: int) -> int:
        return self._col.get(r_idx, -1)

    def __post_init__(self):
        self._row = {s: i for i, s in enumerate(self.sp_indices)}
        self._col = {r: j for j, r in enumerate(self.R_X)}


def _mask_to_sorted_indices(mask: int) -> list[int]:
    out = []
    m = mask
    while m:
        lsb = m & (-m)
        out.append(lsb.bit_length() - 1)
        m &= m - 1
    return out


def so_domain(sp_mask: int, rn_data, S_full: np.ndarray) -> SODomain:
    """
    Build the LP-ready domain for semi-organization `sp_mask`.

    R_X is derived purely from bitmask arithmetic against `rn_data.supp_raw`
    (RAW, not quotiented — a reaction's real support may include E0
    species, which are now part of X_full and must be checked against).
    """
    X_full = sp_mask | rn_data.E0_mask
    sp_indices = _mask_to_sorted_indices(X_full)

    R_X = [r for r, supp in enumerate(rn_data.supp_raw) if (supp & ~X_full) == 0]

    S = S_full[np.ix_(sp_indices, R_X)] if sp_indices and R_X else np.zeros(
        (len(sp_indices), len(R_X))
    )
    return SODomain(sp_mask=sp_mask, X_full_mask=X_full, sp_indices=sp_indices, R_X=R_X, S=S)
