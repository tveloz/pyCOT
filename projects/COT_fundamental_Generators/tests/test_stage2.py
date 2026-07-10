"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  test_stage2.py — Tests for Stage H (hierarchy) and Stage S (synergy)     ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT IT CHECKS
--------------
  Stage H — ERC containment hierarchy on all gold networks:
    • The hierarchy is a valid partial order (antisymmetric + transitive).
    • mask_i ⊊ mask_j  ⟺  j ∈ ancestors[i]  (containment matches bits).
    • Hasse edges are the transitive reduction (no skipped intermediate ERC).

  Stage S — Synergy computation on all gold networks:
    • compute_basic_synergies matches the brute-force synergy oracle exactly.
    • Every maximal synergy is also a basic synergy.
    • Every fundamental synergy is also a maximal synergy.

  Integration — Biomodel 237 (real network, no fixed expected values):
    • hierarchy builds without error; all parent edges are true ancestors.
    • synergy computation completes without error.

  Hypothesis property tests (if hypothesis is installed):
    • oracle == optimised on 200 random small networks.
    • hierarchy antisymmetry holds on 200 random networks.

WHAT IT OUTPUTS
---------------
  PASS / FAIL per test.  For BIOMD237 synergy, counts are printed.

HOW TO RUN
----------
  Option A — VS Code play button:  click ▶ on this file.
  Option B — terminal:  python -m pytest tests/test_stage2.py -v -s
"""

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# Filter: run only tests whose name contains this string.
# Leave empty ("") to run ALL tests.
# Examples:  "hierarchy"  |  "synergy"  |  "oracle"  |  "biomd"  |  "gold2"
FILTER = ""

# Show print output from inside tests (e.g. synergy counts for BIOMD237)
VERBOSE = True
# ─────────────────────────────────────────────────────────────────────────────

from __future__ import annotations
import pytest
from itertools import combinations

from cot_gen.hierarchy  import build_hierarchy
from cot_gen.synergy    import compute_synergies, compute_basic_synergies
from oracles.synergy_oracle import synergy_set as oracle_synergy_set
from tests.gold_networks    import ALL_GOLD

try:
    from hypothesis import given, settings, HealthCheck
    from hypothesis import strategies as st
    HAS_HYPOTHESIS = True
except ImportError:
    HAS_HYPOTHESIS = False


# ── Helpers ───────────────────────────────────────────────────────────────────

def _build_ercs_from_gold(net):
    from cot_gen.erc       import compute_ercs
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
        supp_raw=tuple(net.supp), prod_raw=tuple(net.prod),
        E0_mask=net.E0_mask, supp_q=supp_q, prod_q=prod_q,
        species_to_reactions=inv,
    )
    return compute_ercs(rnd)


# ── Stage H: hierarchy tests ──────────────────────────────────────────────────

@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_hierarchy_partial_order_gold(net):
    """Hierarchy must be a valid partial order: antisymmetric, transitive, mask-consistent."""
    ercs  = _build_ercs_from_gold(net)
    hier  = build_hierarchy(ercs)
    masks = [e.species_mask for e in ercs]
    n     = len(ercs)

    for i, j in combinations(range(n), 2):
        in_anc_i = j in hier.ancestors[i]
        in_anc_j = i in hier.ancestors[j]
        assert not (in_anc_i and in_anc_j), (
            f"[{net.name}] antisymmetry violated: {i}<{j} and {j}<{i}"
        )

    for i in range(n):
        for k in hier.ancestors[i]:
            for j in hier.ancestors[k]:
                assert j in hier.ancestors[i], (
                    f"[{net.name}] transitivity: {i}<{k}<{j} but {j} not ancestor of {i}"
                )

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            mi, mj = masks[i], masks[j]
            is_strict = ((mi & mj) == mi) and mi != mj
            in_anc    = j in hier.ancestors[i]
            assert is_strict == in_anc, (
                f"[{net.name}] mask containment vs ancestors mismatch: i={i} j={j}"
            )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_hasse_is_transitive_reduction(net):
    """Hasse parent edges must be the transitive reduction (no skipped intermediate)."""
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)

    for i in range(hier.n):
        for p in hier.parents[i]:
            for k in hier.ancestors[i]:
                if k != p and p in hier.ancestors[k]:
                    pytest.fail(
                        f"[{net.name}] Hasse edge {i}→{p} is not direct: "
                        f"intermediate {k} exists"
                    )


# ── Stage S: synergy tests ────────────────────────────────────────────────────

@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_basic_synergy_oracle_vs_opt_gold(net):
    """compute_basic_synergies must match the oracle exactly on every gold network."""
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)

    opt_set = {(s.i, s.j, s.k) for s in compute_basic_synergies(ercs, hier)}
    orc_set = oracle_synergy_set(ercs)

    assert opt_set == orc_set, (
        f"[{net.name}] basic synergy mismatch:\n"
        f"  opt only:    {opt_set - orc_set}\n"
        f"  oracle only: {orc_set - opt_set}"
    )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_maximal_is_subset_of_basic_gold(net):
    """Every maximal synergy must appear in the basic synergy set."""
    ercs   = _build_ercs_from_gold(net)
    hier   = build_hierarchy(ercs)
    result = compute_synergies(ercs, hier, level="maximal")
    basic_set = {(s.i, s.j, s.k) for s in result.basic}
    for s in result.maximal:
        assert (s.i, s.j, s.k) in basic_set, (
            f"[{net.name}] maximal ({s.i},{s.j})→{s.k} not in basic"
        )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_fundamental_is_subset_of_maximal_gold(net):
    """Every fundamental synergy must appear in the maximal synergy set."""
    ercs   = _build_ercs_from_gold(net)
    hier   = build_hierarchy(ercs)
    result = compute_synergies(ercs, hier, level="fundamental")
    maximal_set = {(s.i, s.j, s.k) for s in result.maximal}
    for s in result.fundamental:
        assert (s.i, s.j, s.k) in maximal_set, (
            f"[{net.name}] fundamental ({s.i},{s.j})→{s.k} not in maximal"
        )


# ── Integration: Biomodel 237 ─────────────────────────────────────────────────

def test_hierarchy_biomd237(biomd237_rndata):
    """Hierarchy builds without error; parent edges are true ancestors."""
    from cot_gen.erc import compute_ercs
    ercs = compute_ercs(biomd237_rndata)
    hier = build_hierarchy(ercs)
    assert hier.n == len(ercs)
    for i in range(hier.n):
        for p in hier.parents[i]:
            assert p in hier.ancestors[i]


def test_synergy_biomd237(biomd237_rndata):
    """Synergy computation completes on Biomodel 237 without crashing."""
    from cot_gen.erc import compute_ercs
    ercs   = compute_ercs(biomd237_rndata)
    hier   = build_hierarchy(ercs)
    result = compute_synergies(ercs, hier, level="maximal")
    print(f"\nBIOMD237: {len(result.basic)} basic, {len(result.maximal)} maximal")
    assert len(result.basic) >= 0   # just must not crash


# ── Hypothesis property tests ─────────────────────────────────────────────────

pytestmark_hyp = pytest.mark.skipif(
    not HAS_HYPOTHESIS, reason="hypothesis not installed"
)

@st.composite
def _random_net_ercs(draw, max_s=5, max_r=6):
    from cot_gen.erc       import compute_ercs
    from cot_gen.cot_types import RNData
    from oracles.closure_oracle import closure_oracle

    n    = draw(st.integers(min_value=1, max_value=max_s))
    nr   = draw(st.integers(min_value=1, max_value=max_r))
    full = (1 << n) - 1
    supp = draw(st.lists(st.integers(0, full), min_size=nr, max_size=nr))
    prod = draw(st.lists(st.integers(0, full), min_size=nr, max_size=nr))

    inflow_prod = 0
    for s, p in zip(supp, prod):
        if s == 0:
            inflow_prod |= p
    E0     = closure_oracle(supp, prod, inflow_prod)
    supp_q = [s & ~E0 for s in supp]
    prod_q = [p & ~E0 for p in prod]
    inv    = tuple(
        tuple(r for r, s in enumerate(supp_q) if (s >> i) & 1)
        for i in range(n)
    )
    rnd = RNData(
        n_species=n, species_names=tuple(f"s{i}" for i in range(n)),
        species_index=tuple((f"s{i}", i) for i in range(n)),
        n_reactions=nr, reaction_names=tuple(f"r{i}" for i in range(nr)),
        supp_raw=tuple(supp), prod_raw=tuple(prod),
        E0_mask=E0, supp_q=tuple(supp_q), prod_q=tuple(prod_q),
        species_to_reactions=inv,
    )
    return compute_ercs(rnd)


if HAS_HYPOTHESIS:
    @given(_random_net_ercs())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_property_synergy_oracle_vs_opt(ercs):
        hier = build_hierarchy(ercs)
        opt  = {(s.i, s.j, s.k) for s in compute_basic_synergies(ercs, hier)}
        orc  = oracle_synergy_set(ercs)
        assert opt == orc

    @given(_random_net_ercs())
    @settings(max_examples=200, suppress_health_check=[HealthCheck.too_slow])
    def test_property_hierarchy_antisymmetry(ercs):
        hier = build_hierarchy(ercs)
        for i, j in combinations(range(hier.n), 2):
            assert not (j in hier.ancestors[i] and i in hier.ancestors[j])


# ── Play-button entry point ───────────────────────────────────────────────────
if __name__ == "__main__":
    import pytest as _pytest
    _args = [__file__, "-v" if VERBOSE else "-q", "-s"]
    if FILTER:
        _args += ["-k", FILTER]
    raise SystemExit(_pytest.main([a for a in _args if a]))
