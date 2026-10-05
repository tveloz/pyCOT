"""
test_so_search.py — Tests for pyCOT.analysis.organizations.so_search
(elementary SOs / SO hierarchy) and fundamental_graph.py (Mode-1/Mode-2
traversal), plus the conjectures/ strategy comparison.

WHAT IT CHECKS
--------------
  The worked-example network from the companion algorithm paper
  ("A Fully Worked Network", see worked_example.py): 8 reactions over 13
  species, with a sterile ERC (E_dagger, requires a never-produced species)
  and a higher-order SO (M3) that is only reachable via a hierarchy vertical
  lift, not via synergy or minimal-producer complementarity. This is a
  precise regression test for a bug where compute_so_hierarchy's Mode-2 only
  ever extended a found SSM via fundamental-synergy edges, so it silently
  returned zero higher-order SOs on any network whose novelty was
  complementarity- or lift-driven (which the companion paper's own
  statistics say is the common case in real biological networks).

  It is ALSO the flagship illustration of what the paper's two §6.4
  conjectures disagree about (see test_conjecture_independent_seed_misses_M3
  below): run_contained_bfs (Conjecture 1, lift enabled) finds M3;
  run_independent_seed (Conjecture 2, lift disabled) provably cannot,
  because M3 requires absorbing an ancestor ERC that never appears as a
  fundamental synergy/complementarity partner in its own right.

  compute_elementary_sos / compute_so_hierarchy are cross-validated against
  the brute-force oracle (oracles/so_oracle.py) on the gold networks and on
  real BioModels networks small enough for brute force (|ERCs| <= ~24).
  compute_so_hierarchy deliberately excludes pure "latent joins" (unions of
  already-persistent modules with zero synergy/complementarity between
  them) from its search — see so_search.py's latent_join() docstring — so
  the oracle comparison allows for exactly those, reconciled by repeatedly
  applying latent_join() and requiring the result to close the gap
  completely.

HOW TO RUN
----------
  python -m pytest tests/test_so_search.py -v
"""
from __future__ import annotations

import os
import sys
from itertools import combinations

import pytest

from pyCOT.analysis.organizations.cot_types import RNData
from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity import compute_complementarities
from pyCOT.analysis.organizations.so_search import (
    compute_elementary_sos, compute_so_hierarchy, latent_join,
)
from oracles.so_oracle import elementary_so_oracle, so_hierarchy_oracle

from tests.gold_networks import ALL_GOLD
from worked_example import worked_rndata, wmask, names_of
from conjectures.strategies import run_contained_bfs, run_independent_seed

_here = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
_BIOMD91 = os.path.join(_repo, "data", "biochemical_databases", "BioMD_other", "BIOMD0000000091.txt")
_BIOMD999 = os.path.join(_repo, "data", "biochemical_databases", "BioMD_signaling", "BIOMD0000000999.txt")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run_pipeline(rn):
    ercs = compute_ercs(rn)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn_result=syn)
    elem = compute_elementary_sos(rn, ercs, hier, syn_result=syn, comp_result=comp, verbose=False)
    so_hier = compute_so_hierarchy(rn, ercs, hier, syn, comp, elem, max_order=30, verbose=False)
    return ercs, hier, syn, comp, elem, so_hier


def _assert_matches_oracle_modulo_latent_joins(rn, ercs, got_elementary, got_all):
    """
    got_elementary/got_all must be sound (⊆ oracle) and every oracle SO
    absent from got_all must be reconstructible by repeatedly joining pairs
    already in got_all via latent_join() — i.e. the only things
    compute_so_hierarchy omits are pure latent joins, never a genuine
    generated-core member.
    """
    oracle_elementary = set(elementary_so_oracle(rn, ercs))
    oracle_all = set(so_hierarchy_oracle(rn, ercs)["all_sos"])

    assert got_elementary == oracle_elementary, (
        f"elementary-SO mismatch: ours-only={got_elementary - oracle_elementary}  "
        f"oracle-only={oracle_elementary - got_elementary}"
    )
    assert got_all <= oracle_all, f"unsound: found SOs not in oracle: {got_all - oracle_all}"

    missing = oracle_all - got_all
    pool = set(got_all)
    reconstructed: set[int] = set()
    changed = True
    while changed:
        changed = False
        for a, b in combinations(sorted(pool), 2):
            if (a & b) == 0:
                continue
            j = latent_join(rn, a, b)
            if j is not None and j not in pool and j in missing:
                pool.add(j)
                reconstructed.add(j)
                changed = True

    still_missing = missing - reconstructed
    assert not still_missing, (
        f"compute_so_hierarchy omitted genuine generated-core members "
        f"(not explainable as latent joins): {sorted(still_missing)}"
    )


def _build_rndata_from_gold(net):
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
        supp_raw=tuple(net.supp), prod_raw=tuple(net.prod),
        E0_mask=net.E0_mask, supp_q=supp_q, prod_q=prod_q,
        species_to_reactions=inv,
    )


# ---------------------------------------------------------------------------
# Worked-example network (companion algorithm paper, "A Fully Worked
# Network"): E_A, E_B, E_C, E*, E', E_D, E_E, E_dagger.
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def worked():
    rn = worked_rndata()
    ercs, hier, syn, comp, elem, so_hier = _run_pipeline(rn)
    return rn, ercs, hier, syn, comp, elem, so_hier


def test_worked_example_elementary_sos_exact(worked):
    rn, ercs, hier, syn, comp, elem, so_hier = worked
    got = {frozenset(names_of(rn, m)) for m in elem.all_elementary_masks}
    expected = {
        frozenset(["a", "b", "f", "p"]),  # Y1 = E_A join E_B
        frozenset(["a", "c", "f", "p"]),  # Y2 = E_A join E_C
        frozenset(["d", "e", "g", "q"]),  # Y3 = E_D join E_E
    }
    assert got == expected


def test_worked_example_so_hierarchy_matches_oracle_exactly(worked):
    """This network has no pure latent joins, so the generated core equals
    all of SO exactly (no reconciliation needed)."""
    rn, ercs, hier, syn, comp, elem, so_hier = worked
    got_all = set(so_hier.all_so_masks)
    oracle_all = set(so_hierarchy_oracle(rn, ercs)["all_sos"])
    assert got_all == oracle_all


def test_worked_example_vertical_lift_M3_found(worked):
    """M3 = Y1 join E' (vertical lift absorbing E_B).  Before the Mode-2 fix,
    compute_so_hierarchy found zero higher-order SOs on this network because
    M3/M1/M2 all require a non-synergy extension move."""
    rn, ercs, hier, syn, comp, elem, so_hier = worked
    m3 = wmask("a", "b", "f", "p", "z")
    assert m3 in set(so_hier.all_so_masks)


def test_worked_example_discovery_M1_M2_found(worked):
    """M1 = Y1 join E*, M2 = Y2 join E* (complementarity-consumer growth)."""
    rn, ercs, hier, syn, comp, elem, so_hier = worked
    m1 = wmask("a", "b", "f", "h", "p")
    m2 = wmask("a", "c", "f", "h", "p")
    got_all = set(so_hier.all_so_masks)
    assert m1 in got_all
    assert m2 in got_all


def test_worked_example_sterile_erc_excluded(worked):
    """E_dagger = clos{k,w} requires w externally, which nothing produces —
    its species must never appear in any SO."""
    rn, ercs, hier, syn, comp, elem, so_hier = worked
    from worked_example import WORKED_NAMES
    w_bit = WORKED_NAMES.index("w")
    k_bit = WORKED_NAMES.index("k")
    for m in so_hier.all_so_masks:
        assert not ((m >> w_bit) & 1), "sterile species w leaked into a persistent module"
        assert not ((m >> k_bit) & 1), "sterile ERC's species k leaked into a persistent module"


def test_worked_example_orders(worked):
    rn, ercs, hier, syn, comp, elem, so_hier = worked

    def order_of(*names):
        m = wmask(*names)
        if m in set(elem.all_elementary_masks):
            return 0
        for k, masks in so_hier.so_by_order.items():
            if m in masks:
                return k
        return None

    assert order_of("a", "b", "c", "f", "p") == 1
    assert order_of("a", "b", "f", "h", "p") == 1
    assert order_of("a", "c", "f", "h", "p") == 1
    assert order_of("a", "b", "f", "p", "z") == 1
    assert order_of("a", "b", "c", "f", "h", "p") == 2
    assert order_of("a", "b", "c", "f", "p", "z") == 2
    assert order_of("a", "b", "f", "h", "p", "z") == 2
    assert order_of("a", "b", "c", "f", "h", "p", "z") == 3


# ---------------------------------------------------------------------------
# Conjecture 1 vs Conjecture 2 — the concrete divergence
# ---------------------------------------------------------------------------

def test_conjecture_contained_bfs_finds_M3(worked):
    """Conjecture 1 (vertical lift enabled) must find M3 — sanity check that
    the strategies.py wrapper reproduces the underlying so_search result."""
    rn, ercs, hier, syn, comp, elem, so_hier = worked
    run = run_contained_bfs(rn, ercs, hier, syn, comp, max_order=30)
    m3 = wmask("a", "b", "f", "p", "z")
    assert m3 in set(run.all_so_masks)
    assert run.lift_candidates > 0, "contained_bfs should have exercised vertical lift"


def test_conjecture_independent_seed_misses_M3(worked):
    """Conjecture 2 (no vertical lift, independent per-ERC seeding) CANNOT
    find M3: M3 = Y1 ∨ E' requires absorbing E_B's ancestor, and that
    ancestor never appears as a fundamental synergy or complementarity
    partner in its own right, so no purely horizontal chain reaches it.
    This is the concrete, reproducible gap the paper's own §6.4 discussion
    anticipates for the "independent seeding" conjecture."""
    rn, ercs, hier, syn, comp, elem, so_hier = worked
    run = run_independent_seed(rn, ercs, hier, syn, comp, max_order=30)
    m3 = wmask("a", "b", "f", "p", "z")
    assert m3 not in set(run.all_so_masks)
    assert run.lift_candidates == 0, "independent_seed must never use vertical lift"


def test_conjecture_independent_seed_subset_of_contained_bfs(worked):
    """Every SO conjecture 2 finds must also be found by conjecture 1 — lift
    can only ever ADD reachable SOs relative to the horizontal-only search,
    never remove any (both share identical Mode-1 DFS + horizontal Mode-2
    moves; conjecture 1 is a strict superset of conjecture 2's move set)."""
    rn, ercs, hier, syn, comp, elem, so_hier = worked
    lifted = run_contained_bfs(rn, ercs, hier, syn, comp, max_order=30)
    unlifted = run_independent_seed(rn, ercs, hier, syn, comp, max_order=30)
    assert set(unlifted.all_so_masks) <= set(lifted.all_so_masks)
    assert set(unlifted.all_so_masks) < set(lifted.all_so_masks), (
        "expected a strict gap on the worked example (M3 et al.)"
    )


# ---------------------------------------------------------------------------
# Oracle cross-validation — gold networks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_so_oracle_vs_opt_gold(net):
    rn = _build_rndata_from_gold(net)
    ercs = compute_ercs(rn)
    if not ercs:
        pytest.skip(f"{net.name} has no ERCs")
    _, _, _, _, elem, so_hier = _run_pipeline(rn)
    # compute_elementary_sos additionally reports E0 itself as an extra
    # elementary SO whenever E0_mask != 0 (erc.py's module docstring) -- a
    # deliberate departure from Def 30's original "E_∅ = ∅" convention,
    # which the oracle still encodes literally (E0-alone fails its own
    # _is_reactive check, by construction, since inflow reactions have
    # supp_q=0). Excluded here rather than taught to the oracle, same as
    # test_erc.py's test_erc_oracle_vs_opt_gold.
    got_elementary = {m for m in elem.all_elementary_masks if m != rn.E0_mask}
    got_all = {m for m in so_hier.all_so_masks if m != rn.E0_mask}
    _assert_matches_oracle_modulo_latent_joins(rn, ercs, got_elementary, got_all)


# ---------------------------------------------------------------------------
# Oracle cross-validation — a real biological network
# ---------------------------------------------------------------------------

def test_so_oracle_vs_opt_biomd91():
    """
    BIOMD0000000091 (12 ERCs, oracle-feasible): the network on which the
    Mode-2 fix was originally caught (compute_so_hierarchy was silently
    returning fewer higher-order SOs than exist, and briefly — mid-fix —
    returning invalid ones, before the species-level synergy-ignition fix
    in erc_syn_close).
    """
    if not os.path.exists(_BIOMD91):
        pytest.skip(f"BIOMD0000000091 not found at {_BIOMD91}")

    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations.io_pyCOT import build_rndata

    rn_pycot = read_txt(_BIOMD91)
    rn = build_rndata(rn_pycot, network_id="BIOMD0000000091")
    ercs, hier, syn, comp, elem, so_hier = _run_pipeline(rn)

    assert len(ercs) <= 24, "keep this fixture small enough for the brute-force oracle"
    _assert_matches_oracle_modulo_latent_joins(
        rn, ercs, set(elem.all_elementary_masks), set(so_hier.all_so_masks)
    )


def test_so_oracle_vs_opt_biomd999_orphaned_higher_order_carryover():
    """
    BIOMD0000000999 (TGF-beta/SMAD signalling, 8 ERCs, oracle-feasible): the
    network on which the current_layer Mode-2 BFS-seeding bug was caught.

    compute_elementary_sos's own Mode-1 DFS can land directly on an
    order->=1 SSM in a single hop whenever a fundamental complementarity (or
    synergy) pulls in an ERC that -- by itself, as a single ERC -- already
    properly contains a smaller known SO (no intermediate elementary state
    is ever visited along that DFS path). Before the fix, compute_so_hierarchy
    seeded its first BFS round ONLY from elementary_result.all_elementary_masks
    (order 0), silently orphaning any such higher-order carryover from ALL
    further Mode-2 extension -- including vertical lift to its own hierarchy
    parent. On this network: a single fundamental complementarity pulls in an
    11-species P-ERC that is itself already order 2 on arrival (it properly
    contains a smaller order-1, order-0 chain); its hierarchy parent is a
    15-species P-ERC (order 3) reachable ONLY by vertical lift from it. The
    bug made both conjectures (contained_bfs AND independent_seed) silently
    miss the order-3 SO -- this is NOT the expected Conjecture-1-vs-2 gap
    (the tell was both strategies disagreeing with the oracle while still
    agreeing with EACH OTHER).
    """
    if not os.path.exists(_BIOMD999):
        pytest.skip(f"BIOMD0000000999 not found at {_BIOMD999}")

    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations.io_pyCOT import build_rndata

    rn_pycot = read_txt(_BIOMD999)
    rn = build_rndata(rn_pycot, network_id="BIOMD0000000999")
    ercs, hier, syn, comp, elem, so_hier = _run_pipeline(rn)

    assert len(ercs) <= 24, "keep this fixture small enough for the brute-force oracle"
    e0 = rn.E0_mask
    got_elementary = {m for m in elem.all_elementary_masks if m != e0}
    got_all = {m for m in so_hier.all_so_masks if m != e0}
    _assert_matches_oracle_modulo_latent_joins(rn, ercs, got_elementary, got_all)

    # The specific orphaned SO itself, named explicitly so a future
    # regression fails loudly and close to its root cause rather than only
    # via the generic oracle-diff assertion above.
    target = 0x6ed6b5
    assert target in got_all, (
        "the order-3 P-ERC reachable only by vertical lift from an "
        "orphaned order-2 Phase-1 carry-over was not found"
    )
