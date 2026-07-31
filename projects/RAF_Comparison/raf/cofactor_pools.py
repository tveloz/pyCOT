"""
cofactor_pools.py -- RAF-catalyst induction for BiGG-style flux-balance
networks via named cofactor pools, replacing the net-zero-in-one-reaction
heuristic in `biomodel_crs.py` for this network family.

WHY THIS EXISTS (see project memory / earlier conversation): the net-zero
heuristic requires the SAME species token to appear on both sides of ONE
reaction with near-zero net coefficient -- structurally impossible for a
redox/energy cofactor pair like nad_c/nadh_c, since the oxidized and
reduced forms are always two DIFFERENT species ids. That heuristic finds
ZERO catalysts on e_coli_core, not because catalysis is genuinely absent,
but because it is checking the wrong thing: RAF's catalyst notion
(Def 3.1(b), C: R -> 2^M) has NO stoichiometric requirement at all -- a
catalyst only needs to be *present*, and `crs.py`'s `Reaction.catalysts`
field is already free-standing, independent of `reactant_coeffs`/
`product_coeffs` (see raf_algo.py's `is_raf`, which only ever tests
`r.catalysts & W`).

This module instead flags a species as a catalyst of reaction r if r
performs a genuine INTERCONVERSION within a known cofactor pool: one pool
member consumed, a DIFFERENT pool member produced (e.g. nad_c -> nadh_c).
Both members are added to r.catalysts -- on TOP of, not instead of, their
real stoichiometric role, which is left untouched (nad_c stays a genuine
reactant, nadh_c a genuine product; the real chemistry is not hidden, only
annotated). This mirrors Sousa & Hordijk 2015's own "group cofactors into
pools" step (done there via curated UniProt cofactor/metal annotation);
here it is done via BiGG's standardized metabolite-id naming convention,
which recovers the organic-cofactor part of that annotation (not
metals/Fe-S clusters, which need real external per-enzyme curation) from
data already present locally, no new download needed.
"""
from __future__ import annotations

# Fixed pools: each is a closed set of alternate forms of the SAME cofactor.
# Restricted to the cytosolic ('_c') compartment, since that's where the
# large majority of redox/energy chemistry in these models happens; a
# reaction pairing a cytosolic and a periplasmic form of the same cofactor
# (rare) is simply not caught by this fixed-pool list.
COFACTOR_POOLS: list[frozenset[str]] = [
    frozenset({"nad_c", "nadh_c"}),
    frozenset({"nadp_c", "nadph_c"}),
    frozenset({"atp_c", "adp_c", "amp_c"}),
    frozenset({"q8_c", "q8h2_c"}),
    frozenset({"mql8_c", "mqn8_c"}),        # menaquinone/menaquinol, common in anaerobic E. coli models
    frozenset({"fad_c", "fadh2_c"}),
    frozenset({"fmn_c", "fmnh2_c"}),
    frozenset({"thf_c", "5mthf_c", "methf_c", "mlthf_c", "10fthf_c", "5fthf_c", "10fthf5glu3_c"}),
    frozenset({"gtp_c", "gdp_c", "gmp_c"}),
    frozenset({"utp_c", "udp_c", "ump_c"}),
    frozenset({"ctp_c", "cdp_c", "cmp_c"}),
    frozenset({"trdrd_c", "trdox_c"}),      # thioredoxin reduced/oxidized
    frozenset({"gthrd_c", "gthox_c"}),      # glutathione reduced/oxidized
]


def _coa_pool_members(species: set[str]) -> set[str]:
    """CoA is an open-ended pool: coa_c interconverts with ANY acyl-CoA
    thioester (accoa_c, succoa_c, ppcoa_c, malcoa_c, ...). Rather than
    enumerate every acyl-CoA by name, match the '...coa_c' naming suffix
    directly against whatever species the network actually has."""
    return {s for s in species if s.endswith("coa_c")}


def induce_catalysts(reactants: set[str], products: set[str], all_species: set[str]) -> set[str]:
    """Given the RAW reactant/product species-id sets of one reaction
    (coefficients irrelevant here -- only presence matters), return the
    set of species that should be flagged as RAF catalysts of this
    reaction because they participate in a genuine cofactor-pool
    interconversion (one member consumed, a different member produced)."""
    catalysts: set[str] = set()

    for pool in COFACTOR_POOLS:
        in_r = reactants & pool
        in_p = products & pool
        if in_r and in_p and (in_r != in_p or len(in_r) > 1):
            if any(a != b for a in in_r for b in in_p):
                catalysts |= in_r | in_p

    coa_r = _coa_pool_members(reactants)
    coa_p = _coa_pool_members(products)
    if coa_r and coa_p and any(a != b for a in coa_r for b in coa_p):
        catalysts |= coa_r | coa_p

    return catalysts & all_species
