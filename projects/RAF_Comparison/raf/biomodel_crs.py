"""
biomodel_crs.py — Build a CRS from a real BioModels/BiGG network already
loaded through cot_gen (RNData) + the original pyCOT ReactionNetwork.

There is no explicit catalysis annotation in these networks (COT reaction
networks carry only stoichiometry, no separate catalysis relation, Sec 2
of decomposing_RAF_v2.pdf / Def 7 of systems-05-00030-v2.pdf). We INDUCE
catalysis using COT's own definition of a catalyst (Def 2.8 in the newer
paper, Def 7 in the 2017 paper): a species with net-zero stoichiometric
effect that appears on both sides of a single reaction. This is the most
defensible reverse-engineering of "which species act as a RAF catalyst"
available from stoichiometry alone -- it will typically find catalysis to
be RARE in genome-scale flux-balance-style networks (where cofactors like
ATP/ADP are usually regenerated over a CYCLE of reactions, not within one
reaction), which is itself an informative empirical result about the gap
between RAF theory's per-reaction catalysis requirement and COT's broader
per-cycle self-maintenance -- report it, don't hide it by special-casing.

F0 (food set) is the union of products of DIRECT inflow reactions (raw
support == 0) -- matching Def 2.18/Sec 3's F0 exactly, not cot_gen's own
broader E0_mask (which is the full reachability closure from inflow,
already established as a DIFFERENT, wider notion in Decomposition_Theorem;
see that project's bridge.py docstring). Inflow reactions themselves are
not part of R (Sec 3 models food purely via F0, not via a reaction) --
see the RAF-algebra layer's own test networks for the same convention.
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_decomp_proj = os.path.normpath(os.path.join(_here, "..", "..", "Decomposition_Theorem"))
if _decomp_proj not in sys.path:
    sys.path.insert(0, _decomp_proj)

from decomp.bridge import build_full_stoich  # noqa: E402

from .crs import CRS

_TOL = 1e-9


def crs_from_biomodel(rn_pycot, rn_data) -> tuple[CRS, dict]:
    """
    Returns (crs, stats) where stats reports how much induced catalysis was
    actually found (n_reactions, n_with_catalyst, n_food) -- always inspect
    this: a CRS with very few catalyzed reactions will have a tiny or empty
    maxRAF, which is a real result, not a construction bug.
    """
    S_full = build_full_stoich(rn_pycot)
    n_species, n_reactions = rn_data.n_species, rn_data.n_reactions

    F0_mask = 0
    for r in range(n_reactions):
        if rn_data.supp_raw[r] == 0:
            F0_mask |= rn_data.prod_raw[r]
    F0_names = rn_data.bitset_to_names(F0_mask)

    reactions = []
    n_with_catalyst = 0
    n_skipped_inflow = 0
    for r in range(n_reactions):
        if rn_data.supp_raw[r] == 0:
            n_skipped_inflow += 1
            continue
        reactant_coeffs: dict[str, float] = {}
        product_coeffs: dict[str, float] = {}
        catalysts: list[str] = []
        col = S_full[:, r]
        supp_r, prod_r = rn_data.supp_raw[r], rn_data.prod_raw[r]
        for sp_idx in range(n_species):
            supp_bit = (supp_r >> sp_idx) & 1
            prod_bit = (prod_r >> sp_idx) & 1
            if not supp_bit and not prod_bit:
                continue
            name = rn_data.species_name(sp_idx)
            coeff = float(col[sp_idx])
            if supp_bit and prod_bit and abs(coeff) <= _TOL:
                catalysts.append(name)
            else:
                if coeff < -_TOL:
                    reactant_coeffs[name] = -coeff
                elif coeff > _TOL:
                    product_coeffs[name] = coeff
                elif supp_bit:
                    reactant_coeffs[name] = reactant_coeffs.get(name, 0.0) + 1.0
                elif prod_bit:
                    product_coeffs[name] = product_coeffs.get(name, 0.0) + 1.0
        if catalysts:
            n_with_catalyst += 1
        rname = rn_data.reaction_name(r)
        reactions.append((rname, list(reactant_coeffs), list(product_coeffs), catalysts,
                           reactant_coeffs, product_coeffs))

    crs = CRS.build(reactions, food=F0_names)
    stats = {
        "n_reactions_total": n_reactions,
        "n_reactions_in_crs": len(reactions),
        "n_inflow_skipped": n_skipped_inflow,
        "n_with_induced_catalyst": n_with_catalyst,
        "n_food_species": len(F0_names),
    }
    return crs, stats


def crs_from_biomodel_no_catalysis(rn_pycot, rn_data) -> tuple[CRS, dict]:
    """
    Third induction mode: NO catalyst induction at all, C(r) = empty set for
    every reaction. Not a heuristic -- for networks like e_coli_core where
    no per-reaction catalysis is genuinely recoverable from the metabolite
    stoichiometry (see module docstring's cofactor-pools discussion), this
    is the honest baseline: with every C(r) empty, Def 2.1(b) ("there is a
    catalyst c in C(r) with c in W") fails for every reaction, so no
    reaction can ever belong to any RAF -- maxRAF = empty set immediately,
    for every food set, with no search needed to see it. Useful as an
    explicit companion to the induced-catalyst modes above: it isolates
    what the COT/decomposition layer (which never reads `catalysts` at
    all) delivers completely independent of any catalyst-induction guess.
    """
    S_full = build_full_stoich(rn_pycot)
    n_species, n_reactions = rn_data.n_species, rn_data.n_reactions

    F0_mask = 0
    for r in range(n_reactions):
        if rn_data.supp_raw[r] == 0:
            F0_mask |= rn_data.prod_raw[r]
    F0_names = rn_data.bitset_to_names(F0_mask)

    reactions = []
    n_skipped_inflow = 0
    for r in range(n_reactions):
        if rn_data.supp_raw[r] == 0:
            n_skipped_inflow += 1
            continue
        reactant_coeffs: dict[str, float] = {}
        product_coeffs: dict[str, float] = {}
        col = S_full[:, r]
        supp_r, prod_r = rn_data.supp_raw[r], rn_data.prod_raw[r]
        for sp_idx in range(n_species):
            supp_bit = (supp_r >> sp_idx) & 1
            prod_bit = (prod_r >> sp_idx) & 1
            if not supp_bit and not prod_bit:
                continue
            name = rn_data.species_name(sp_idx)
            coeff = float(col[sp_idx])
            if coeff < -_TOL:
                reactant_coeffs[name] = -coeff
            elif coeff > _TOL:
                product_coeffs[name] = coeff
            elif supp_bit and prod_bit:
                reactant_coeffs[name] = reactant_coeffs.get(name, 0.0) + 1.0
                product_coeffs[name] = product_coeffs.get(name, 0.0) + 1.0
            elif supp_bit:
                reactant_coeffs[name] = reactant_coeffs.get(name, 0.0) + 1.0
            elif prod_bit:
                product_coeffs[name] = product_coeffs.get(name, 0.0) + 1.0
        rname = rn_data.reaction_name(r)
        reactions.append((rname, list(reactant_coeffs), list(product_coeffs), [],
                           reactant_coeffs, product_coeffs))

    crs = CRS.build(reactions, food=F0_names)
    stats = {
        "n_reactions_total": n_reactions,
        "n_reactions_in_crs": len(reactions),
        "n_inflow_skipped": n_skipped_inflow,
        "n_with_induced_catalyst": 0,
        "n_food_species": len(F0_names),
    }
    return crs, stats


def crs_from_biomodel_cofactor_pools(rn_pycot, rn_data) -> tuple[CRS, dict]:
    """
    Second induction mode, for BiGG-style flux-balance networks where the
    net-zero-in-one-reaction heuristic above structurally cannot find
    anything (see that function's docstring and cofactor_pools.py): a
    species is flagged as a RAF catalyst of reaction r if r performs a
    genuine cofactor-pool interconversion (nad_c -> nadh_c, atp_c -> adp_c,
    coa_c <-> any *coa_c, etc. -- see cofactor_pools.COFACTOR_POOLS).

    Unlike the net-zero mode, this does NOT remove the flagged species from
    reactant_coeffs/product_coeffs -- the real chemistry (nad_c genuinely
    consumed, nadh_c genuinely produced) is preserved exactly, and the
    catalyst annotation is added on top, using the fact that `Reaction.
    catalysts` is a free-standing field independent of stoichiometry (see
    crs.py / raf_algo.py's `is_raf`, which only ever tests `r.catalysts & W`
    and never reads the coefficient dicts).
    """
    from .cofactor_pools import induce_catalysts

    S_full = build_full_stoich(rn_pycot)
    n_species, n_reactions = rn_data.n_species, rn_data.n_reactions

    F0_mask = 0
    for r in range(n_reactions):
        if rn_data.supp_raw[r] == 0:
            F0_mask |= rn_data.prod_raw[r]
    F0_names = rn_data.bitset_to_names(F0_mask)
    all_species = {rn_data.species_name(i) for i in range(n_species)}

    reactions = []
    n_with_catalyst = 0
    n_skipped_inflow = 0
    for r in range(n_reactions):
        if rn_data.supp_raw[r] == 0:
            n_skipped_inflow += 1
            continue
        reactant_coeffs: dict[str, float] = {}
        product_coeffs: dict[str, float] = {}
        col = S_full[:, r]
        supp_r, prod_r = rn_data.supp_raw[r], rn_data.prod_raw[r]
        for sp_idx in range(n_species):
            supp_bit = (supp_r >> sp_idx) & 1
            prod_bit = (prod_r >> sp_idx) & 1
            if not supp_bit and not prod_bit:
                continue
            name = rn_data.species_name(sp_idx)
            coeff = float(col[sp_idx])
            if coeff < -_TOL:
                reactant_coeffs[name] = -coeff
            elif coeff > _TOL:
                product_coeffs[name] = coeff
            elif supp_bit and prod_bit:
                # genuinely net-zero within this one reaction (rare, but
                # still a real reactant+product pair, not a pool artifact)
                reactant_coeffs[name] = reactant_coeffs.get(name, 0.0) + 1.0
                product_coeffs[name] = product_coeffs.get(name, 0.0) + 1.0
            elif supp_bit:
                reactant_coeffs[name] = reactant_coeffs.get(name, 0.0) + 1.0
            elif prod_bit:
                product_coeffs[name] = product_coeffs.get(name, 0.0) + 1.0

        catalysts = sorted(induce_catalysts(
            set(reactant_coeffs), set(product_coeffs), all_species))
        if catalysts:
            n_with_catalyst += 1
        rname = rn_data.reaction_name(r)
        reactions.append((rname, list(reactant_coeffs), list(product_coeffs), catalysts,
                           reactant_coeffs, product_coeffs))

    crs = CRS.build(reactions, food=F0_names)
    stats = {
        "n_reactions_total": n_reactions,
        "n_reactions_in_crs": len(reactions),
        "n_inflow_skipped": n_skipped_inflow,
        "n_with_induced_catalyst": n_with_catalyst,
        "n_food_species": len(F0_names),
    }
    return crs, stats
