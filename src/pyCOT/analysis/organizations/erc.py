"""
erc.py — ERC discovery, minimal bases, and per-ERC req/prod.

Public API
----------
compute_ercs(rn_data, counters=None) -> list[ERCData]
    Stage-1 main function.  For every non-trivial quotiented reaction r
    (supp_q[r] ≠ 0), compute erc_mask = closure_opt(supp_q[r]).
    Group reactions by their ERC.  For each group compute MinBas and req/prod.
    Also injects one extra ERC for E0 (the closure of the inflow), when
    E0_mask != 0 — see _e0_erc() below.
    Returns ERCData objects sorted by species_mask (smallest first).

    Metrics emitted (via counters, all prefixed "erc."):
      closures_computed   — number of closures run
      unique_ercs         — number of distinct ERCs found
      min_bases_total     — total MinBas entries across all ERCs
      req_prod_scans      — reactions scanned for req/prod per ERC (total)

Key assertions:
  • closure(min_base) == erc.species_mask for every min_base of every ERC
  • species_mask is closed (is_closed(supp_q, prod_q, mask))
  • If req_mask == 0: the ERC is a P-ERC (persistent)
  • MinBas forms an antichain (no element ⊆ another)

See oracles/erc_oracle.py for the brute-force reference.

E0 as a P-ERC
-------------
E0 (the closure of the union of all inflow species) is, by construction,
always closed and always semi-self-maintaining (req_mask == 0 relative to
the RAW network: the only reactions "active" within E0 are the inflows
themselves, which consume nothing). It is therefore always a genuine
persistent ERC in its own right whenever the network has an inflow — this
was previously left unrepresented, since every OTHER ERC is computed via
QUOTIENTED closure (supp_q/prod_q, which strip E0's species out entirely),
so E0's own species never appeared as their own ERC/EPM node.

E0 does not disqualify, or get disqualified by, any other ERC's
elementary/EPM status: every other ERC's species_mask lives in the
quotiented (E0-excluded) index space, so E0's bits are structurally never
set in any other ERC's mask — the two live in disjoint parts of the
bitset. That is the right behavior, not an approximation: E0's species are
unconditionally available, so no other ERC's req_mask can ever name one
(nothing is ever "missing" an E0 species), which is exactly why E0 never
participates in synergy or complementarity either (see synergy.py /
complementarity.py — both operate purely on quotiented req/prod masks).
E0 is "elementary in a generative sense" relative to every other ERC: the
only persistent ERC any other one could be said to contain is E0 itself,
and that containment doesn't count against being elementary, since E0 is
free by construction.
"""
from __future__ import annotations

from typing import Sequence

from .cot_types import RNData, ERCData
from .closure import closure_opt, build_inv_idx, is_closed


# ---------------------------------------------------------------------------
# MinBas computation
# ---------------------------------------------------------------------------

def _minimal_subset(masks: list[int]) -> list[int]:
    """
    Return the ⊆-antichain: ⊆-minimal elements of `masks`.

    Duplicates are removed first.  O(n²) — acceptable because the number
    of distinct supports per ERC is small in practice.
    """
    unique = list(dict.fromkeys(masks))  # deduplicate preserving order
    minimal = []
    for i, m_i in enumerate(unique):
        dominated = False
        for j, m_j in enumerate(unique):
            if i != j and (m_j & m_i) == m_j and m_j != m_i:
                # m_j ⊊ m_i → m_i is not minimal
                dominated = True
                break
        if not dominated:
            minimal.append(m_i)
    return minimal


# ---------------------------------------------------------------------------
# req / prod for a closed set
# ---------------------------------------------------------------------------

def _req_prod(supp_q: Sequence[int], prod_q: Sequence[int], mask: int):
    """
    Compute (req_mask, prod_mask) for the closed set `mask`.

    R_mask = {r : supp_q[r] ⊆ mask}
    prod_mask = ⋃ prod_q[r] for r ∈ R_mask
    req_mask  = (⋃ supp_q[r] for r ∈ R_mask) & ~prod_mask
    """
    agg_supp = 0
    agg_prod = 0
    scans = 0
    for s, p in zip(supp_q, prod_q):
        scans += 1
        if (s & ~mask) == 0:      # supp ⊆ mask
            agg_supp |= s
            agg_prod |= p
    req = agg_supp & ~agg_prod
    return req, agg_prod, scans


# ---------------------------------------------------------------------------
# E0 as a P-ERC (see module docstring)
# ---------------------------------------------------------------------------

def _e0_erc(rn_data: RNData, erc_id: int) -> ERCData:
    """
    Build the ERCData for E0 (the closure of the union of all inflow
    species). Always persistent (req_mask=0) and always closed, by
    construction of E0_mask itself (see cot_types.RNData / io_pyCOT).

    min_bases=[] deliberately, not [0]: E0 needs no reactant precondition
    at all, and an empty min_bases list makes every existing invariant
    check in _assert_erc_invariants (all of which are `for b in
    e.min_bases: ...` loops) hold vacuously, correctly — no special-casing
    needed. closure_opt (which computes ordinary ERCs' closures over
    QUOTIENTED supp_q/prod_q) cannot reconstruct E0_mask from any non-empty
    base anyway, since E0's species are exactly what quotienting strips out.
    """
    e0 = rn_data.E0_mask
    reaction_indices = [
        r for r in range(rn_data.n_reactions)
        if (rn_data.supp_raw[r] & ~e0) == 0
    ]
    return ERCData(
        erc_id=erc_id,
        species_mask=e0,
        reaction_indices=reaction_indices,
        min_bases=[],
        req_mask=0,
        prod_mask=e0,
    )


# ---------------------------------------------------------------------------
# Main Stage-1 function
# ---------------------------------------------------------------------------

def compute_ercs(
    rn_data: RNData,
    counters=None,          # optional cot_gen.metrics.Counters
    verify: bool = True,    # set False to skip closure-based assertions (faster on large nets)
) -> list[ERCData]:
    """
    Discover all ERCs and compute their minimal bases and req/prod.

    Steps
    -----
    1. Build inverted index from rn_data.species_to_reactions.
    2. For each reaction r with supp_q[r] ≠ 0:
       a. erc_mask = closure_opt(supp_q[r], starting from supp_q[r])
       b. Record (erc_mask, r_idx, supp_q[r]) in a dict.
    3. For each unique erc_mask group:
       a. MinBas = ⊆-minimal elements of {supp_q[r] : r in group}.
       b. req_mask, prod_mask from _req_prod(erc_mask).
       c. Build ERCData.
    4. Assert structural invariants.
    5. Return sorted list (by species_mask, for determinism).

    Parameters
    ----------
    rn_data  : RNData (Stage-0 output)
    counters : optional Counters for instrumentation

    Returns
    -------
    list[ERCData] sorted by species_mask (ascending).
    """
    supp_q = rn_data.supp_q
    prod_q = rn_data.prod_q
    inv_idx = list(rn_data.species_to_reactions)  # list of tuples → list of lists ok

    # --- Step 2: compute ERC for every non-trivial reaction -----------------
    # erc_groups: erc_mask → (list of reaction indices, list of their supp_q masks)
    erc_groups: dict[int, tuple[list[int], list[int]]] = {}

    closures_computed = 0
    for r_idx, s in enumerate(supp_q):
        if s == 0:
            continue   # trivial (inflow or fully-E0) — skip
        erc_mask = closure_opt(supp_q, prod_q, s, inv_idx, counters=counters)
        closures_computed += 1
        if erc_mask not in erc_groups:
            erc_groups[erc_mask] = ([], [])
        erc_groups[erc_mask][0].append(r_idx)
        erc_groups[erc_mask][1].append(s)

    if counters is not None:
        counters.inc("erc.closures_computed", closures_computed)
        counters.inc("erc.unique_ercs", len(erc_groups))

    # --- Step 3: build ERCData objects --------------------------------------
    ercs: list[ERCData] = []
    total_req_prod_scans = 0
    total_min_bases = 0

    for erc_mask, (r_indices, supports) in sorted(erc_groups.items()):
        min_bases = _minimal_subset(supports)
        req_mask, prod_mask, scans = _req_prod(supp_q, prod_q, erc_mask)
        total_req_prod_scans += scans
        total_min_bases += len(min_bases)

        ercs.append(ERCData(
            erc_id=-1,   # reassigned below, once E0 (if any) is in the list
            species_mask=erc_mask,
            reaction_indices=r_indices,
            min_bases=min_bases,
            req_mask=req_mask,
            prod_mask=prod_mask,
        ))

    # E0 (see module docstring): always closed, always persistent, never
    # disqualifies or is disqualified by any other ERC's elementary status
    # (disjoint quotiented bitset). Injected here, not among erc_groups
    # above, since it isn't discovered via closure_opt over supp_q/prod_q
    # at all -- it's read directly off rn_data.E0_mask.
    if rn_data.E0_mask != 0:
        ercs.append(_e0_erc(rn_data, erc_id=-1))

    ercs.sort(key=lambda e: e.species_mask)
    for erc_id, e in enumerate(ercs):
        e.erc_id = erc_id

    if counters is not None:
        counters.inc("erc.min_bases_total", total_min_bases)
        counters.inc("erc.req_prod_scans", total_req_prod_scans)

    # --- Step 4: structural assertions (skip with verify=False for large nets)
    if verify:
        _assert_erc_invariants(ercs, supp_q, prod_q, inv_idx)

    return ercs


# ---------------------------------------------------------------------------
# Invariant checks
# ---------------------------------------------------------------------------

def _assert_erc_invariants(
    ercs: list[ERCData],
    supp_q: Sequence[int],
    prod_q: Sequence[int],
    inv_idx: list,
) -> None:
    """
    Assert key ERC structural invariants.  Called after compute_ercs.

    Activation identity: r ∈ R_E  ⟺  supp_q[r] ⊆ E.species_mask
    Closure idempotent:  closure(min_base) == E.species_mask  ∀ min_base
    MinBas antichain:    no min_base is a strict subset of another
    Closure correctness: E.species_mask is closed
    """
    # Build a map erc_mask → ERCData for quick lookup
    mask_to_erc: dict[int, ERCData] = {e.species_mask: e for e in ercs}

    for e in ercs:
        mask = e.species_mask

        # Closed set check
        assert is_closed(supp_q, prod_q, mask), \
            f"ERC id={e.erc_id} species_mask={mask:#x} is not closed"

        # Each min_base must be contained in the ERC
        for b in e.min_bases:
            assert (b & ~mask) == 0, \
                f"ERC id={e.erc_id}: min_base {b:#x} not ⊆ species_mask {mask:#x}"

        # closure(min_base) == erc_mask
        for b in e.min_bases:
            got = closure_opt(supp_q, prod_q, b, inv_idx)
            assert got == mask, (
                f"ERC id={e.erc_id}: closure(min_base={b:#x})={got:#x} "
                f"≠ species_mask={mask:#x}"
            )

        # MinBas is an antichain (no element ⊆ another)
        for i, b_i in enumerate(e.min_bases):
            for j, b_j in enumerate(e.min_bases):
                if i != j:
                    assert not ((b_j & b_i) == b_j and b_j != b_i), (
                        f"ERC id={e.erc_id}: MinBas is not antichain — "
                        f"{b_j:#x} ⊊ {b_i:#x}"
                    )

        # Reaction membership consistency
        for r in e.reaction_indices:
            assert (supp_q[r] & ~mask) == 0, \
                f"ERC id={e.erc_id}: reaction {r} has supp_q not ⊆ species_mask"
