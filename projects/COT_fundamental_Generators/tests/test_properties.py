"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  test_properties.py — Hypothesis-based property tests for cot_gen          ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT IT CHECKS
--------------
Using the Hypothesis property-testing library, this file generates hundreds
of random small reaction networks and verifies mathematical invariants that
must hold universally (not just on the gold networks):

  1. Closure oracle == optimised (Horn/Dowling-Gallier):
     Both algorithms must return the same closure on every seed set.

  2. Closure idempotence:
     cl(cl(X)) = cl(X)  for all X.

  3. Closure monotonicity:
     X ⊆ cl(X)  for all X.

  4. ERC masks are closed:
     For every ERC mask E, is_closed(supp_q, prod_q, E) must hold.

  5. MinBas antichain:
     No element of MinBas(E) is a strict subset of another element.

  6. E0 invariant:
     After correct E0 computation, supp_q[r] = 0  ⟹  prod_q[r] = 0.

REQUIRES
--------
  pip install hypothesis

  If hypothesis is not installed, ALL tests are skipped automatically
  (no failure — just "skipped" in the pytest output).

WHAT IT OUTPUTS
---------------
  PASS / FAIL / SKIP per property.  Hypothesis prints a minimal failing
  example if any property is violated.

HOW TO RUN
----------
  Option A — VS Code play button:  click ▶ on this file.
  Option B — terminal:  python -m pytest tests/test_properties.py -v -s

TUNING
------
  Increase MAX_EXAMPLES below to find rarer failures (slower).
  Decrease to speed up development runs.
"""
from __future__ import annotations
import pytest

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# Filter: run only tests whose name contains this string.
# Leave empty ("") to run ALL property tests.
FILTER = ""

# Show output (Hypothesis statistics summary printed at end of run with -s flag)
VERBOSE = True

# Number of random examples Hypothesis generates per property.
# More examples → higher confidence, but slower.
MAX_EXAMPLES = 300
# ─────────────────────────────────────────────────────────────────────────────

try:
    from hypothesis import given, settings, HealthCheck
    from hypothesis import strategies as st
    HAS_HYPOTHESIS = True
except ImportError:
    HAS_HYPOTHESIS = False

from pyCOT.analysis.organizations.closure    import closure_opt, build_inv_idx, is_closed
from oracles.closure_oracle import closure_oracle
from pyCOT.analysis.organizations.erc        import compute_ercs
from pyCOT.analysis.organizations.cot_types  import RNData

pytestmark = pytest.mark.skipif(
    not HAS_HYPOTHESIS,
    reason="hypothesis not installed — run: pip install hypothesis",
)


# ── Strategy: random reaction network ─────────────────────────────────────────

@st.composite
def random_net(draw, max_s=6, max_r=8):
    """Generate (supp, prod, n_species, seed) for a small random network."""
    n    = draw(st.integers(min_value=1, max_value=max_s))
    nr   = draw(st.integers(min_value=1, max_value=max_r))
    full = (1 << n) - 1
    supp = draw(st.lists(st.integers(0, full), min_size=nr, max_size=nr))
    prod = draw(st.lists(st.integers(0, full), min_size=nr, max_size=nr))
    seed = draw(st.integers(0, full))
    return supp, prod, n, seed


def _make_rndata(supp, prod, n, E0_mask=0):
    """Package raw lists into a minimal RNData for compute_ercs."""
    supp_q = [s & ~E0_mask for s in supp]
    prod_q = [p & ~E0_mask for p in prod]
    inv = tuple(
        tuple(r for r, s in enumerate(supp_q) if (s >> i) & 1)
        for i in range(n)
    )
    return RNData(
        n_species=n,
        species_names=tuple(f"s{i}" for i in range(n)),
        species_index=tuple((f"s{i}", i) for i in range(n)),
        n_reactions=len(supp),
        reaction_names=tuple(f"r{i}" for i in range(len(supp))),
        supp_raw=tuple(supp), prod_raw=tuple(prod),
        E0_mask=E0_mask, supp_q=tuple(supp_q), prod_q=tuple(prod_q),
        species_to_reactions=inv,
    )


# ── Property tests ────────────────────────────────────────────────────────────

if HAS_HYPOTHESIS:

    @given(random_net())
    @settings(max_examples=MAX_EXAMPLES, suppress_health_check=[HealthCheck.too_slow])
    def test_property_closure_oracle_vs_opt(args):
        """The optimised closure must equal the oracle on every random input."""
        supp, prod, n, seed = args
        inv      = build_inv_idx(supp, n)
        expected = closure_oracle(supp, prod, seed)
        got      = closure_opt(supp, prod, seed, inv)
        assert got == expected, (
            f"oracle={expected:#b} opt={got:#b} seed={seed:#b} n={n}"
        )

    @given(random_net())
    @settings(max_examples=MAX_EXAMPLES, suppress_health_check=[HealthCheck.too_slow])
    def test_property_closure_idempotent(args):
        """cl(cl(X)) = cl(X) for all seeds."""
        supp, prod, n, seed = args
        inv = build_inv_idx(supp, n)
        c   = closure_opt(supp, prod, seed, inv)
        c2  = closure_opt(supp, prod, c, inv)
        assert c == c2

    @given(random_net())
    @settings(max_examples=MAX_EXAMPLES, suppress_health_check=[HealthCheck.too_slow])
    def test_property_closure_monotone(args):
        """X ⊆ cl(X): the seed is always contained in its own closure."""
        supp, prod, n, seed = args
        inv = build_inv_idx(supp, n)
        c   = closure_opt(supp, prod, seed, inv)
        assert (c & seed) == seed

    @given(random_net())
    @settings(max_examples=MAX_EXAMPLES, suppress_health_check=[HealthCheck.too_slow])
    def test_property_erc_masks_are_closed(args):
        """Every ERC mask must be a fixed point of the closure operator."""
        supp, prod, n, _ = args
        inflow_prod = 0
        for s, p in zip(supp, prod):
            if s == 0:
                inflow_prod |= p
        E0  = closure_oracle(supp, prod, inflow_prod)
        rnd = _make_rndata(supp, prod, n, E0_mask=E0)
        for e in compute_ercs(rnd):
            assert is_closed(rnd.supp_q, rnd.prod_q, e.species_mask), (
                f"ERC mask {e.species_mask:#b} is not closed"
            )

    @given(random_net())
    @settings(max_examples=MAX_EXAMPLES, suppress_health_check=[HealthCheck.too_slow])
    def test_property_minbas_antichain(args):
        """MinBas(E) must be an antichain: no element strictly contains another."""
        supp, prod, n, _ = args
        inflow_prod = 0
        for s, p in zip(supp, prod):
            if s == 0:
                inflow_prod |= p
        E0  = closure_oracle(supp, prod, inflow_prod)
        rnd = _make_rndata(supp, prod, n, E0_mask=E0)
        for e in compute_ercs(rnd):
            bases = e.min_bases
            for i, b_i in enumerate(bases):
                for j, b_j in enumerate(bases):
                    if i != j:
                        assert not ((b_j & b_i) == b_j and b_j != b_i), (
                            f"MinBas not antichain in ERC {e.erc_id}: "
                            f"{b_j:#b} ⊊ {b_i:#b}"
                        )

    @given(random_net())
    @settings(max_examples=MAX_EXAMPLES, suppress_health_check=[HealthCheck.too_slow])
    def test_property_stage0_E0_invariant(args):
        """After E0 computation: supp_q[r]=0  ⟹  prod_q[r]=0."""
        supp, prod, n, _ = args
        inflow_prod = 0
        for s, p in zip(supp, prod):
            if s == 0:
                inflow_prod |= p
        E0 = closure_oracle(supp, prod, inflow_prod)
        for r, (s, p) in enumerate(zip(supp, prod)):
            sq = s & ~E0
            pq = p & ~E0
            if sq == 0:
                assert pq == 0, (
                    f"Reaction {r}: supp_q=0 but prod_q={pq:#b} (E0={E0:#b})"
                )


# ── Play-button entry point ───────────────────────────────────────────────────
if __name__ == "__main__":
    import pytest as _pytest
    _args = [__file__, "-v" if VERBOSE else "-q", "-s"]
    if FILTER:
        _args += ["-k", FILTER]
    raise SystemExit(_pytest.main([a for a in _args if a]))
