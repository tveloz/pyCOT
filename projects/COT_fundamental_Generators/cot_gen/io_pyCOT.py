"""
io_pyCOT.py — Bridge from pyCOT ReactionNetwork to cot_gen RNData.

Public API
----------
build_rndata(rn, network_id="") -> RNData
    Convert a pyCOT ReactionNetwork to an RNData (Stage-0 output).
    Steps:
      1. Enumerate species and assign indices 0..n_species-1.
      2. Build supp_raw, prod_raw as bitmasks from reaction objects.
      3. Compute E0 = closure_oracle(supp_raw, prod_raw, inflow_product_mask).
      4. Build supp_q, prod_q = raw & ~E0_mask.
      5. Build species_to_reactions inverted index from supp_q.
      6. Assert Stage-0 invariants.
      7. Return RNData.

load_rndata(path, *, network_id=None) -> RNData
    Convenience: read_txt(path) then build_rndata.

Note: pyCOT uses rustworkx internally; names may carry " (OriginalName)"
suffixes in some biomodels (cosmetic only — we store raw names as-is).
"""
from __future__ import annotations

import os
import sys

# Layout: pyCOT/projects/COT_fundamental_Generators/cot_gen/io_pyCOT.py
# _here  = .../cot_gen/
# _proj  = .../COT_fundamental_Generators/   (cot_gen + oracles importable from here)
# _repo  = .../pyCOT/                        (pyCOT importable from _repo/src)
_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_src  = os.path.normpath(os.path.join(_here, "..", "..", "..", "src"))
for _p in (_proj, _src):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyCOT.io.functions import read_txt           # noqa: E402

from .cot_types import RNData
from oracles.closure_oracle import closure_oracle  # noqa: E402


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _species_list(rn):
    """
    Return a list of (index, name) pairs for all species in `rn`.

    pyCOT's rn.species() returns Species objects with .index and .name.
    We sort by .index so our sequential assignment matches pyCOT's.
    """
    return sorted(rn.species(), key=lambda s: s.index)


def _build_masks(rn, sp_to_idx: dict[str, int]):
    """
    Build supp_raw and prod_raw bitmasks from pyCOT reaction objects.

    pyCOT Reaction objects expose:
      r.support_names()  -> list[str]
      r.products_names() -> list[str]
      r.node             -> ReactionNode(index=int, name=str, ...)
    Note: r.name is a bound method, not an attribute; use r.node.name.
    """
    reactions = sorted(rn.reactions(), key=lambda r: r.node.index)
    supp_raw = []
    prod_raw = []
    reaction_names = []
    for rxn in reactions:
        s_mask = 0
        for sp_name in rxn.support_names():
            if sp_name in sp_to_idx:
                s_mask |= 1 << sp_to_idx[sp_name]
        p_mask = 0
        for sp_name in rxn.products_names():
            if sp_name in sp_to_idx:
                p_mask |= 1 << sp_to_idx[sp_name]
        supp_raw.append(s_mask)
        prod_raw.append(p_mask)
        reaction_names.append(rxn.node.name)
    return supp_raw, prod_raw, reaction_names


def _inflow_product_mask(supp_raw: list[int], prod_raw: list[int]) -> int:
    """
    Collect the union of products of all inflow reactions (supp = 0).
    These species are available 'for free' and seed E0.
    """
    mask = 0
    for s, p in zip(supp_raw, prod_raw):
        if s == 0:
            mask |= p
    return mask


def _build_inv_idx(supp_q: list[int], n_species: int) -> list[list[int]]:
    inv: list[list[int]] = [[] for _ in range(n_species)]
    for r_idx, s in enumerate(supp_q):
        m = s
        while m:
            lsb = m & (-m)
            inv[lsb.bit_length() - 1].append(r_idx)
            m &= m - 1
    return inv


# ---------------------------------------------------------------------------
# Main public API
# ---------------------------------------------------------------------------

def build_rndata(rn, network_id: str = "") -> RNData:
    """
    Convert a pyCOT ReactionNetwork to an RNData (Stage-0).

    Parameters
    ----------
    rn         : pyCOT ReactionNetwork (from read_txt).
    network_id : string identifier used in metrics rows.

    Returns
    -------
    RNData with all Stage-0 fields populated and invariants verified.
    """
    # --- 1. Species universe -------------------------------------------------
    sp_list = _species_list(rn)
    n_species = len(sp_list)
    species_names = tuple(sp.name for sp in sp_list)
    # Build our own index (compact 0..n_species-1) by sorted order
    sp_to_idx: dict[str, int] = {sp.name: i for i, sp in enumerate(sp_list)}
    species_index = tuple(sorted(sp_to_idx.items(), key=lambda kv: kv[1]))

    # --- 2. Reactions --------------------------------------------------------
    supp_raw_list, prod_raw_list, reaction_names = _build_masks(rn, sp_to_idx)
    n_reactions = len(supp_raw_list)

    # --- 3. Compute E0 -------------------------------------------------------
    seed = _inflow_product_mask(supp_raw_list, prod_raw_list)
    E0_mask = closure_oracle(supp_raw_list, prod_raw_list, seed)

    # --- 4. Quotient ---------------------------------------------------------
    supp_q_list = [s & ~E0_mask for s in supp_raw_list]
    prod_q_list = [p & ~E0_mask for p in prod_raw_list]

    # --- 5. Inverted index ---------------------------------------------------
    inv_idx = _build_inv_idx(supp_q_list, n_species)
    species_to_reactions = tuple(tuple(lst) for lst in inv_idx)

    # --- 6. Assemble RNData --------------------------------------------------
    rn_data = RNData(
        n_species=n_species,
        species_names=species_names,
        species_index=species_index,
        n_reactions=n_reactions,
        reaction_names=tuple(reaction_names),
        supp_raw=tuple(supp_raw_list),
        prod_raw=tuple(prod_raw_list),
        E0_mask=E0_mask,
        supp_q=tuple(supp_q_list),
        prod_q=tuple(prod_q_list),
        species_to_reactions=species_to_reactions,
    )

    # --- 7. Assert Stage-0 invariants ----------------------------------------
    rn_data.assert_stage0_invariants()

    return rn_data


def load_rndata(path: str, *, network_id: str | None = None) -> RNData:
    """
    Load a .txt reaction network file and return RNData.

    Parameters
    ----------
    path       : path to the .txt file (pyCOT format).
    network_id : if None, derived from the file basename.
    """
    if network_id is None:
        network_id = os.path.splitext(os.path.basename(path))[0]
    rn = read_txt(path)
    return build_rndata(rn, network_id=network_id)
