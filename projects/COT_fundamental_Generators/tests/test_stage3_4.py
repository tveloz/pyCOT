"""
╔══════════════════════════════════════════════════════════════════════════════╗
║  test_stage3_4.py — Tests for Stage C (complementarity) and Stage G       ║
║                     (generators / primitive ERCs) and the MetaNetwork      ║
╚══════════════════════════════════════════════════════════════════════════════╝

WHAT IT CHECKS
--------------
  PART 1 — Complementarity  (cot_gen/complementarity.py):
    1a. Oracle vs optimised — exact set equality on all gold networks.
    1b. No complementarity when all ERC pairs are comparable.
    1c. GOLD2 specific known results: pair (E_s0s1, E_s2s3) has Type 1 and Type 3.
    1d. Custom networks with hand-crafted Type 1 and Type 3 results.
    1e. Mathematical invariants:
        • All reported pairs are incomparable.
        • Type 2 is always empty (subsumed by Type 1 — see paper for proof).

  PART 2 — Generators  (cot_gen/generators.py):
    2a. Gold network results (GOLD1: single primitive, GOLD2: two primitives +
        synergy reach covers all, GOLD5: all primitive).
    2b. Mathematical properties:
        • No primitive ERC is a fundamental synergy target.
        • All primitives are in the basis reach.
        • Coverage ∈ [0, 1].
        • Unreachable and basis_reach are disjoint and together cover all ERCs.
    2c. reachable_from helper: empty seed, no synergies, chain of synergies.

  PART 3 — MetaNetwork  (cot_gen/metanetwork.py):
    3a. Builds on all gold networks without error.
    3b. Node list has one entry per ERC with required fields.
    3c. Edge list contains only valid rel_type values.
    3d. Stats dict is consistent with underlying component counts.
    3e. Integration on Biomodel 237.

WHAT IT OUTPUTS
---------------
  PASS / FAIL per test, plus printed stats for Biomodel 237 integration tests.

HOW TO RUN
----------
  Option A — VS Code play button:  click ▶ on this file.
  Option B — terminal:  python -m pytest tests/test_stage3_4.py -v -s
"""

# ── CONFIGURATION ─────────────────────────────────────────────────────────────
# Filter: run only tests whose name contains this string.
# Leave empty ("") to run ALL tests in this file.
# Examples:  "comp"  |  "generators"  |  "metanetwork"  |  "gold2"  |  "biomd"
FILTER = ""

# Show print output from inside tests (e.g. BIOMD237 MetaNetwork stats)
VERBOSE = True
# ─────────────────────────────────────────────────────────────────────────────

from __future__ import annotations
import pytest
from itertools import combinations   # noqa: F401 (used by hypothesis tests)

from cot_gen.erc              import compute_ercs
from cot_gen.hierarchy        import build_hierarchy
from cot_gen.synergy          import compute_synergies
from cot_gen.complementarity  import compute_complementarities, CompTuple   # noqa: F401
from cot_gen.generators       import compute_generators, reachable_from, primitive_ercs
from cot_gen.metanetwork      import build_metanetwork
from cot_gen.cot_types        import RNData
from oracles.complementarity_oracle import comp_oracle, comp_set as oracle_comp_set
from tests.gold_networks import ALL_GOLD, GOLD1_LOOP, GOLD2_HIERARCHY, GOLD5_NONPERSISTENT


# ── Helpers ───────────────────────────────────────────────────────────────────

def _build_rndata(n_species, species, supp, prod, E0_mask=0):
    supp_q = tuple(s & ~E0_mask for s in supp)
    prod_q = tuple(p & ~E0_mask for p in prod)
    inv = tuple(
        tuple(r for r, s in enumerate(supp_q) if (s >> i) & 1)
        for i in range(n_species)
    )
    return RNData(
        n_species=n_species,
        species_names=tuple(species),
        species_index=tuple((name, i) for i, name in enumerate(species)),
        n_reactions=len(supp_q),
        reaction_names=tuple(f"r{i}" for i in range(len(supp_q))),
        supp_raw=tuple(supp), prod_raw=tuple(prod),
        E0_mask=E0_mask, supp_q=supp_q, prod_q=prod_q,
        species_to_reactions=inv,
    )


def _build_ercs_from_gold(net):
    rnd = _build_rndata(net.n_species, net.species, net.supp, net.prod, net.E0_mask)
    return compute_ercs(rnd)


def _all_stages(net):
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    comp = compute_complementarities(ercs, hier)
    gen  = compute_generators(ercs, hier, syn)
    return ercs, hier, syn, comp, gen


# ── Custom mini-networks ──────────────────────────────────────────────────────

def _two_loops_ercs():
    """
    Two independent loops → two incomparable persistent ERCs → Type 3.
    Species: s0(0), s1(1), s2(2), s3(3)
    r0: s0→s1, r1: s1→s0  →  ERC_0={s0,s1}, req=∅, prod={s0,s1}
    r2: s2→s3, r3: s3→s2  →  ERC_1={s2,s3}, req=∅, prod={s2,s3}
    """
    rnd = _build_rndata(
        4, ["s0","s1","s2","s3"],
        supp=[0b0001, 0b0010, 0b0100, 0b1000],
        prod=[0b0010, 0b0001, 0b1000, 0b0100],
    )
    return compute_ercs(rnd)


def _type1_ercs():
    """
    Two ERCs with shared requirements → Type 1.
    r0: {s0,s1}→{s2}  →  ERC_0={s0,s1,s2}, req={s0,s1}
    r1: {s1,s3}→{s2}  →  ERC_1={s1,s2,s3}, req={s1,s3}
    Joint covers both: joint_req=∅ < req_0+req_1=4 → Type 1.
    """
    rnd = _build_rndata(
        4, ["s0","s1","s2","s3"],
        supp=[0b0011, 0b1010],   # {s0,s1}=3, {s1,s3}=10
        prod=[0b0100, 0b0100],   # both produce s2=4
    )
    return compute_ercs(rnd)


# ══════════════════════════════════════════════════════════════════════════════
# PART 1: COMPLEMENTARITY
# ══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_comp_oracle_vs_opt_gold(net):
    """compute_complementarities must match comp_oracle exactly on every gold network."""
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)

    orc = oracle_comp_set(ercs)
    opt = {(c.i, c.j, c.comp_type, c.direction) for c in compute_complementarities(ercs, hier).all}

    assert opt == orc, (
        f"[{net.name}] mismatch:\n  opt only: {opt-orc}\n  oracle only: {orc-opt}"
    )


def test_comp_empty_on_single_erc():
    """Only one ERC → no pairs → no complementarity."""
    ercs = _build_ercs_from_gold(GOLD1_LOOP)
    comp = compute_complementarities(ercs, build_hierarchy(ercs))
    assert len(comp) == 0


def test_comp_empty_on_comparable_pair():
    """GOLD5: two comparable ERCs → no complementarity."""
    ercs = _build_ercs_from_gold(GOLD5_NONPERSISTENT)
    comp = compute_complementarities(ercs, build_hierarchy(ercs))
    assert len(comp) == 0


def test_comp_gold2_type1():
    """
    GOLD2 pair (E_s0s1, E_s2s3): joint_req=0 < req_i+req_j → Type 1.
    """
    ercs = _build_ercs_from_gold(GOLD2_HIERARCHY)
    comp = compute_complementarities(ercs, build_hierarchy(ercs))

    masks   = {e.species_mask: idx for idx, e in enumerate(ercs)}
    i_s0s1  = masks[0b0011]
    i_s2s3  = masks[0b1100]
    i_p, j_p = min(i_s0s1, i_s2s3), max(i_s0s1, i_s2s3)

    assert (i_p, j_p) in {(c.i, c.j) for c in comp.type1}, (
        f"Expected Type 1 for ({i_p},{j_p}), got {comp.type1}"
    )


def test_comp_gold2_type3():
    """
    GOLD2 pair (E_s0s1, E_s2s3): E_s2s3 extends E_s0s1's production → Type 3.
    """
    ercs = _build_ercs_from_gold(GOLD2_HIERARCHY)
    comp = compute_complementarities(ercs, build_hierarchy(ercs))

    masks   = {e.species_mask: idx for idx, e in enumerate(ercs)}
    i_s0s1  = masks[0b0011]
    i_s2s3  = masks[0b1100]
    i_p, j_p = min(i_s0s1, i_s2s3), max(i_s0s1, i_s2s3)

    assert any(c.i == i_p and c.j == j_p for c in comp.type3), (
        f"Expected Type 3 for ({i_p},{j_p}), got {comp.type3}"
    )


def test_comp_two_loops_type3():
    """Two independent persistent loops → Type 3 (both directions) only."""
    ercs = _two_loops_ercs()
    hier = build_hierarchy(ercs)
    comp = compute_complementarities(ercs, hier)

    assert len(ercs) == 2
    assert not hier.is_comparable(0, 1)
    assert len(comp.type3) == 1, f"Expected 1 Type-3, got {comp.type3}"
    assert len(comp.type1) == 0
    assert comp.type3[0].direction == 0   # both directions hold simultaneously


def test_comp_type1_example():
    """Custom Type-1 network: shared requirement covered by joint closure."""
    ercs = _type1_ercs()
    hier = build_hierarchy(ercs)
    comp = compute_complementarities(ercs, hier)

    assert len(ercs) == 2
    assert not hier.is_comparable(0, 1)
    assert len(comp.type1) >= 1


def test_comp_type1_oracle_vs_opt_custom():
    ercs = _type1_ercs()
    hier = build_hierarchy(ercs)
    orc  = oracle_comp_set(ercs)
    opt  = {(c.i, c.j, c.comp_type, c.direction)
            for c in compute_complementarities(ercs, hier).all}
    assert opt == orc


def test_comp_two_loops_oracle_vs_opt():
    ercs = _two_loops_ercs()
    hier = build_hierarchy(ercs)
    orc  = oracle_comp_set(ercs)
    opt  = {(c.i, c.j, c.comp_type, c.direction)
            for c in compute_complementarities(ercs, hier).all}
    assert opt == orc


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_comp_only_incomparable_pairs(net):
    """All reported complementarity pairs must be incomparable."""
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)
    comp = compute_complementarities(ercs, hier)

    for c in comp.all:
        assert not hier.is_comparable(c.i, c.j), (
            f"[{net.name}] comparable pair ({c.i},{c.j}) has complementarity"
        )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_comp_type2_always_empty_gold(net):
    """
    Type 2 is mathematically unreachable:
    req_i ∩ req_j ≠ ∅  always implies  |joint_req| < |req_i|+|req_j|  → Type 1 fires first.
    """
    ercs = _build_ercs_from_gold(net)
    comp = compute_complementarities(ercs, build_hierarchy(ercs))
    assert len(comp.type2) == 0, f"[{net.name}] unexpected Type-2: {comp.type2}"


def test_comp_type2_always_empty_custom():
    ercs = _two_loops_ercs()
    assert len(compute_complementarities(ercs, build_hierarchy(ercs)).type2) == 0


# ══════════════════════════════════════════════════════════════════════════════
# PART 2: GENERATORS
# ══════════════════════════════════════════════════════════════════════════════

def test_generators_empty_net():
    """GOLD4 (no ERCs): empty GeneratorResult with coverage = 1.0."""
    from tests.gold_networks import GOLD4_INFLOW
    ercs, hier, syn, comp, gen = _all_stages(GOLD4_INFLOW)
    assert gen.primitive_indices == []
    assert gen.basis_reach == frozenset()
    assert gen.coverage == 1.0
    assert gen.is_complete()


def test_generators_single_erc():
    """GOLD1 (one ERC, no synergies): that ERC is the sole primitive."""
    ercs = _build_ercs_from_gold(GOLD1_LOOP)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    gen  = compute_generators(ercs, hier, syn)

    assert gen.primitive_indices == [0]
    assert gen.basis_reach == frozenset({0})
    assert gen.coverage == 1.0
    assert gen.is_complete()


def test_generators_gold2():
    """
    GOLD2: fundamental synergy (E1,E2)→E3.
    Primitives = {E1, E2}, basis reach = all 3, coverage = 1.0.
    """
    ercs = _build_ercs_from_gold(GOLD2_HIERARCHY)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    gen  = compute_generators(ercs, hier, syn)

    masks     = {e.species_mask: idx for idx, e in enumerate(ercs)}
    idx_s0s1  = masks[0b0011]
    idx_s2s3  = masks[0b1100]
    idx_full  = masks[0b1111]

    fund_targets = {s.k for s in syn.fundamental}
    assert idx_full in fund_targets, "E3 should be a fundamental synergy target"
    assert sorted(gen.primitive_indices) == sorted([idx_s0s1, idx_s2s3])
    assert idx_full in gen.basis_reach
    assert gen.coverage == 1.0
    assert gen.is_complete()


def test_generators_gold5_all_primitive():
    """GOLD5: two comparable ERCs, no synergies → both are primitive."""
    ercs = _build_ercs_from_gold(GOLD5_NONPERSISTENT)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    gen  = compute_generators(ercs, hier, syn)

    assert len(syn.fundamental) == 0
    assert sorted(gen.primitive_indices) == [0, 1]
    assert gen.coverage == 1.0


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_generators_primitives_not_syn_targets(net):
    """No primitive ERC may be the target of a fundamental synergy."""
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    gen  = compute_generators(ercs, hier, syn)

    fund_targets = {s.k for s in syn.fundamental}
    for p in gen.primitive_indices:
        assert p not in fund_targets, (
            f"[{net.name}] primitive ERC {p} is a fundamental synergy target"
        )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_generators_primitives_in_basis_reach(net):
    """Every primitive ERC must be in the basis reach."""
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    gen  = compute_generators(ercs, hier, syn)

    for p in gen.primitive_indices:
        assert p in gen.basis_reach, f"[{net.name}] primitive {p} not in basis_reach"


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_generators_coverage_in_01(net):
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    gen  = compute_generators(ercs, hier, syn)
    assert 0.0 <= gen.coverage <= 1.0


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_generators_unreachable_disjoint_from_reach(net):
    """Unreachable ∪ basis_reach = all ERCs, and they are disjoint."""
    ercs = _build_ercs_from_gold(net)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    gen  = compute_generators(ercs, hier, syn)

    all_idx = frozenset(range(len(ercs)))
    assert gen.basis_reach | gen.unreachable == all_idx
    assert gen.basis_reach & gen.unreachable == frozenset()


# ── reachable_from helper ─────────────────────────────────────────────────────

def test_reachable_from_empty_seed():
    from cot_gen.synergy import SynergyTuple
    triples = [SynergyTuple(i=0, j=1, k=2, level="fundamental")]
    assert reachable_from([], triples) == frozenset()


def test_reachable_from_no_synergies():
    assert reachable_from([0, 1], []) == frozenset({0, 1})


def test_reachable_from_chain():
    """Chain (0,1)→2, (2,1)→3: starting from {0,1} reaches {0,1,2,3}."""
    from cot_gen.synergy import SynergyTuple
    triples = [
        SynergyTuple(i=0, j=1, k=2, level="fundamental"),
        SynergyTuple(i=2, j=1, k=3, level="fundamental"),
    ]
    assert reachable_from([0, 1], triples) == frozenset({0, 1, 2, 3})


# ══════════════════════════════════════════════════════════════════════════════
# PART 3: METANETWORK
# ══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_metanetwork_builds_on_gold(net):
    ercs, hier, syn, comp, gen = _all_stages(net)
    assert build_metanetwork(ercs, hier, syn, comp, gen) is not None


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_metanetwork_node_list_length(net):
    """Node list must have exactly one entry per ERC."""
    ercs, hier, syn, comp, gen = _all_stages(net)
    mn    = build_metanetwork(ercs, hier, syn, comp, gen)
    nodes = mn.to_node_list()
    assert len(nodes) == len(ercs), f"[{net.name}] {len(nodes)} nodes ≠ {len(ercs)} ERCs"


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_metanetwork_node_fields(net):
    """Each node dict must contain the required keys."""
    ercs, hier, syn, comp, gen = _all_stages(net)
    mn       = build_metanetwork(ercs, hier, syn, comp, gen)
    required = {
        "idx", "erc_id", "size", "is_persistent",
        "is_primitive", "in_basis_reach",
        "n_reactions", "n_min_bases", "req_popcount", "prod_popcount",
    }
    for node in mn.to_node_list():
        assert required.issubset(node.keys()), (
            f"[{net.name}] missing fields: {required - node.keys()}"
        )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_metanetwork_edge_list_types(net):
    """All edge rel_type values must be in the known set."""
    valid = {"hasse", "syn_basic", "syn_maximal", "syn_fundamental",
             "comp_type1", "comp_type2", "comp_type3"}
    ercs, hier, syn, comp, gen = _all_stages(net)
    mn = build_metanetwork(ercs, hier, syn, comp, gen)
    for edge in mn.to_edge_list():
        assert edge["rel_type"] in valid, (
            f"[{net.name}] unknown rel_type: {edge['rel_type']}"
        )


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_metanetwork_stats_consistent(net):
    """Stats dict values must match the underlying component counts."""
    ercs, hier, syn, comp, gen = _all_stages(net)
    mn = build_metanetwork(ercs, hier, syn, comp, gen)
    s  = mn.stats()

    assert s["n_ercs"]            == len(ercs)
    assert s["n_basic_syn"]       == len(syn.basic)
    assert s["n_maximal_syn"]     == len(syn.maximal)
    assert s["n_fundamental_syn"] == len(syn.fundamental)
    assert s["n_comp_type1"]      == len(comp.type1)
    assert s["n_comp_type2"]      == len(comp.type2)
    assert s["n_comp_type3"]      == len(comp.type3)
    assert s["n_comp_total"]      == len(comp)
    assert s["n_primitives"]      == len(gen.primitive_indices)
    assert 0.0 <= s["coverage"]   <= 1.0
    assert s["n_persistent_ercs"] + s["n_non_persistent"] == s["n_ercs"]
    assert s["n_comparable"]      + s["n_incomparable"]   == s["n_pairs"]


@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_metanetwork_hasse_edges_count(net):
    """Hasse edge count in stats must match the sum of parent-edge counts."""
    ercs, hier, syn, comp, gen = _all_stages(net)
    mn = build_metanetwork(ercs, hier, syn, comp, gen)
    s  = mn.stats()
    direct = sum(len(hier.parents[i]) for i in range(hier.n))
    assert s["n_hasse_edges"] == direct


def test_metanetwork_biomd237(biomd237_rndata):
    """Build a complete MetaNetwork for Biomodel 237 without crashing."""
    ercs = compute_ercs(biomd237_rndata)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    comp = compute_complementarities(ercs, hier)
    gen  = compute_generators(ercs, hier, syn)
    mn   = build_metanetwork(ercs, hier, syn, comp, gen)

    s = mn.stats()
    print(f"\nBIOMD237: {s['n_ercs']} ERCs | {s['n_fundamental_syn']} fundamental syn | "
          f"{s['n_comp_total']} comp | {s['n_primitives']} primitives | "
          f"coverage={s['coverage']:.1%}")
    assert s["n_ercs"] > 0
    assert s["coverage"] >= 0.0
    assert len(mn.to_node_list()) == s["n_ercs"]


def test_metanetwork_print_summary_biomd237(biomd237_rndata):
    """print_summary should not raise on Biomodel 237."""
    ercs = compute_ercs(biomd237_rndata)
    hier = build_hierarchy(ercs)
    syn  = compute_synergies(ercs, hier, level="fundamental")
    comp = compute_complementarities(ercs, hier)
    gen  = compute_generators(ercs, hier, syn)
    build_metanetwork(ercs, hier, syn, comp, gen).print_summary()


# ── Play-button entry point ───────────────────────────────────────────────────
if __name__ == "__main__":
    import pytest as _pytest
    _args = [__file__, "-v" if VERBOSE else "-q", "-s"]
    if FILTER:
        _args += ["-k", FILTER]
    raise SystemExit(_pytest.main([a for a in _args if a]))
