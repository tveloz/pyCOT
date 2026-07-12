"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  test_erc.py — Tests for Stage 0 (RNData) and Stage 1 (ERC discovery)     ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT IT CHECKS
--------------
  Stage 0 — RNData invariants on Biomodel 237 (real network from BioModels DB):
    • supp_q[r] = 0  ⟹  prod_q[r] = 0  (E0 quotient invariant)

  Stage 1 — ERC computation on gold networks (hand-crafted with known answers):
    • ERC count matches expected.
    • ERC masks match expected bitmasks.
    • Min-base counts match expected.
    • Persistence flag (req = ∅) matches expected.

  Cross-validation — efficient algorithm vs brute-force oracle:
    • ERC masks must be identical on all gold networks and on Biomodel 237.
    • Min-base counts per ERC must also match.

WHAT IT OUTPUTS
---------------
  PASS / FAIL status for each test case (shown by pytest with -v flag).
  With VERBOSE = True (see below) print statements inside tests are also shown.

HOW TO RUN
----------
  Option A — VS Code play button:  click ▶ on this file.
  Option B — terminal:  python -m pytest tests/test_erc.py -v -s
"""
from __future__ import annotations
import pytest

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# Filter: run only tests whose name contains this string.
# Leave empty ("") to run ALL tests in this file.
# Example: FILTER = "gold"   → only gold-network tests
#          FILTER = "biomd"  → only Biomodel-237 tests
#          FILTER = "oracle" → only oracle-vs-optimised tests
FILTER = ""

# Show print output from inside test functions (e.g. ERC details for BIOMD237)
VERBOSE = True
# ─────────────────────────────────────────────────────────────────────────────

from cot_gen.erc         import compute_ercs
from cot_gen.io_pyCOT    import build_rndata
from oracles.erc_oracle   import compute_ercs_oracle
from tests.gold_networks  import ALL_GOLD


# ── Helpers ───────────────────────────────────────────────────────────────────

def _build_rndata_from_gold(net):
    """Build a minimal RNData from a gold network descriptor."""
    from cot_gen.cot_types import RNData
    supp_q = tuple(s & ~net.E0_mask for s in net.supp)
    prod_q = tuple(p & ~net.E0_mask for p in net.prod)
    inv = tuple(
        tuple(r for r, s in enumerate(supp_q) if (s >> i) & 1)
        for i in range(net.n_species)
    )
    return RNData(
        n_species=net.n_species,
        species_names=tuple(net.species),
        species_index=tuple((name, i) for i, name in enumerate(net.species)),
        n_reactions=len(supp_q),
        reaction_names=tuple(f"r{i}" for i in range(len(supp_q))),
        supp_raw=tuple(net.supp),
        prod_raw=tuple(net.prod),
        E0_mask=net.E0_mask,
        supp_q=supp_q,
        prod_q=prod_q,
        species_to_reactions=inv,
    )


# ── Stage 1: ERC counts and structure on gold networks ────────────────────────

@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_erc_count_gold(net):
    """Number, masks, min-base counts, and persistence must match the gold table."""
    rnd  = _build_rndata_from_gold(net)
    ercs = compute_ercs(rnd)

    assert len(ercs) == len(net.expected_ercs), (
        f"[{net.name}] expected {len(net.expected_ercs)} ERCs, got {len(ercs)}\n"
        f"  got:      {[hex(e.species_mask) for e in ercs]}\n"
        f"  expected: {[hex(m) for m,_,_ in net.expected_ercs]}"
    )

    by_mask = {e.species_mask: e for e in ercs}
    for erc_mask, n_min_bases, persistent in net.expected_ercs:
        assert erc_mask in by_mask, (
            f"[{net.name}] ERC mask {erc_mask:#x} not found"
        )
        e = by_mask[erc_mask]
        assert len(e.min_bases) == n_min_bases, (
            f"[{net.name}] ERC {erc_mask:#x}: "
            f"expected {n_min_bases} min_bases, got {len(e.min_bases)}"
        )
        assert e.is_persistent() == persistent, (
            f"[{net.name}] ERC {erc_mask:#x}: is_persistent expected {persistent}"
        )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_erc_oracle_vs_opt_gold(net):
    """compute_ercs must produce same ERC masks and min-base counts as oracle."""
    from cot_gen.cot_types import RNData
    supp_q = tuple(s & ~net.E0_mask for s in net.supp)
    prod_q = tuple(p & ~net.E0_mask for p in net.prod)
    inv = tuple(
        tuple(r for r, s in enumerate(supp_q) if (s >> i) & 1)
        for i in range(net.n_species)
    )
    rnd = RNData(
        n_species=net.n_species,
        species_names=tuple(net.species),
        species_index=tuple((name, i) for i, name in enumerate(net.species)),
        n_reactions=len(supp_q),
        reaction_names=tuple(f"r{i}" for i in range(len(supp_q))),
        supp_raw=tuple(net.supp),
        prod_raw=tuple(net.prod),
        E0_mask=net.E0_mask,
        supp_q=supp_q, prod_q=prod_q,
        species_to_reactions=inv,
    )
    ercs_opt    = compute_ercs(rnd)
    ercs_oracle = compute_ercs_oracle(list(supp_q), list(prod_q))

    masks_opt = sorted(e.species_mask for e in ercs_opt)
    masks_orc = sorted(e["species_mask"] for e in ercs_oracle)
    assert masks_opt == masks_orc, (
        f"[{net.name}] ERC masks differ: opt={masks_opt} oracle={masks_orc}"
    )
    for e_opt in ercs_opt:
        e_orc = next(e for e in ercs_oracle if e["species_mask"] == e_opt.species_mask)
        assert len(e_opt.min_bases) == len(e_orc["min_bases"]), (
            f"[{net.name}] ERC {e_opt.species_mask:#x}: "
            f"min_bases opt={len(e_opt.min_bases)} oracle={len(e_orc['min_bases'])}"
        )


# ── Stage 0 + Stage 1: Biomodel 237 integration ───────────────────────────────

def test_stage0_invariants_biomd237(biomd237_rndata):
    """E0 quotient invariant: supp_q[r]=0 ⟹ prod_q[r]=0."""
    biomd237_rndata.assert_stage0_invariants()


def test_erc_count_biomd237(biomd237_rndata):
    """Biomodel 237 must have at least one ERC."""
    ercs = compute_ercs(biomd237_rndata)
    assert len(ercs) >= 1, "Expected at least one ERC"
    print(f"\nBIOMD237: {len(ercs)} ERCs")
    for e in ercs:
        print(f"  E{e.erc_id}: size={e.size()} rxns={len(e.reaction_indices)} "
              f"minbas={len(e.min_bases)} persistent={e.is_persistent()}")


def test_erc_oracle_vs_opt_biomd237(biomd237_rndata):
    """ERC masks for Biomodel 237 must match between oracle and optimised."""
    rn       = biomd237_rndata
    ercs_opt = compute_ercs(rn)
    ercs_orc = compute_ercs_oracle(list(rn.supp_q), list(rn.prod_q))

    masks_opt = sorted(e.species_mask for e in ercs_opt)
    masks_orc = sorted(e["species_mask"] for e in ercs_orc)
    assert masks_opt == masks_orc, "BIOMD237: ERC masks differ"


# ── Play-button entry point ───────────────────────────────────────────────────
if __name__ == "__main__":
    import pytest as _pytest
    _args = [__file__, "-v" if VERBOSE else "-q", "-s"]
    if FILTER:
        _args += ["-k", FILTER]
    raise SystemExit(_pytest.main([a for a in _args if a]))
