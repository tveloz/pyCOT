"""
test_closure.py — Tests for cot_gen/closure.py.

Checks:
  • oracle vs optimised: must agree on all gold networks and Biomodel 237
  • idempotence: closure(closure(X)) == closure(X)
  • monotonicity: X ⊆ closure(X)
  • is_closed / is_ssm agree with oracle variants
"""
from __future__ import annotations

import pytest

from pyCOT.analysis.organizations.closure import closure_opt, build_inv_idx, is_closed, is_ssm
from oracles.closure_oracle import closure_oracle, is_closed_oracle, is_ssm_oracle
from tests.gold_networks import ALL_GOLD


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _all_subsets(n: int):
    """Yield all bitmasks 0..(2^n - 1)."""
    for v in range(1 << n):
        yield v


# ---------------------------------------------------------------------------
# Gold network tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_closure_oracle_vs_opt_gold(net):
    """Optimized closure == oracle on every subset of every gold network."""
    inv = build_inv_idx(net.supp, net.n_species)
    for X in _all_subsets(net.n_species):
        expected = closure_oracle(net.supp, net.prod, X)
        got = closure_opt(net.supp, net.prod, X, inv)
        assert got == expected, (
            f"[{net.name}] X={X:#b}: oracle={expected:#b} opt={got:#b}"
        )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_closure_idempotent_gold(net):
    inv = build_inv_idx(net.supp, net.n_species)
    for X in _all_subsets(net.n_species):
        c = closure_opt(net.supp, net.prod, X, inv)
        c2 = closure_opt(net.supp, net.prod, c, inv)
        assert c == c2, f"[{net.name}] not idempotent at X={X:#b}"


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_closure_monotone_gold(net):
    inv = build_inv_idx(net.supp, net.n_species)
    for X in _all_subsets(net.n_species):
        c = closure_opt(net.supp, net.prod, X, inv)
        assert (c & X) == X, f"[{net.name}] X not ⊆ closure(X) at X={X:#b}"


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_is_closed_agrees_gold(net):
    for X in _all_subsets(net.n_species):
        assert is_closed(net.supp, net.prod, X) == is_closed_oracle(net.supp, net.prod, X)


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_is_ssm_agrees_gold(net):
    for X in _all_subsets(net.n_species):
        assert is_ssm(net.supp, net.prod, X) == is_ssm_oracle(net.supp, net.prod, X)


# ---------------------------------------------------------------------------
# Biomodel 237 tests
# ---------------------------------------------------------------------------

def test_closure_oracle_vs_opt_biomd237(biomd237_rndata):
    rn = biomd237_rndata
    inv = list(rn.species_to_reactions)
    # test on each reaction's supp_q as seed (not all 2^26 subsets — too large)
    for r in range(rn.n_reactions):
        X = rn.supp_q[r]
        expected = closure_oracle(rn.supp_q, rn.prod_q, X)
        got = closure_opt(rn.supp_q, rn.prod_q, X, inv)
        assert got == expected, f"reaction {r}: oracle≠opt"


def test_closure_idempotent_biomd237(biomd237_rndata):
    rn = biomd237_rndata
    inv = list(rn.species_to_reactions)
    for r in range(rn.n_reactions):
        X = rn.supp_q[r]
        c = closure_opt(rn.supp_q, rn.prod_q, X, inv)
        c2 = closure_opt(rn.supp_q, rn.prod_q, c, inv)
        assert c == c2, f"reaction {r}: closure not idempotent"
