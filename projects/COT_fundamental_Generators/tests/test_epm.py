"""
test_epm.py — Tests for cot_gen/epm.py (EPMs and ESPMs) and
cot_gen/fundamental_graph.py (Mode-1/Mode-2 traversal).

WHAT IT CHECKS
--------------
  The worked-example network from the companion algorithm paper
  ("A Fully Worked Network"): 8 reactions over 13 species, with a sterile
  ERC (E_dagger, requires a never-produced species) and an ESPM (M3) that
  is only reachable via a hierarchy vertical lift, not via synergy or
  minimal-producer complementarity.  This is a precise regression test for
  a bug where compute_espm's Mode-2 only ever extended a found SSM via
  fundamental-synergy edges, so it silently returned zero ESPMs on any
  network whose novelty was complementarity- or lift-driven (which the
  companion paper's own statistics say is the common case in real
  biological networks).

  compute_epms / compute_espm are cross-validated against the brute-force
  oracle (oracles/epm_oracle.py) on the gold networks and on real BioModels
  networks small enough for brute force (|ERCs| <= ~24).  compute_espm
  deliberately excludes pure "latent joins" (unions of already-persistent
  modules with zero synergy/complementarity between them) from its search —
  see epm.py's latent_join() docstring — so the oracle comparison allows
  for exactly those, reconciled by repeatedly applying latent_join() and
  requiring the result to close the gap completely.

HOW TO RUN
----------
  python -m pytest tests/test_epm.py -v
"""
from __future__ import annotations

import os
import sys
from itertools import combinations

import pytest

from cot_gen.cot_types import RNData
from cot_gen.erc import compute_ercs
from cot_gen.hierarchy import build_hierarchy
from cot_gen.synergy import compute_synergies_basis_first
from cot_gen.complementarity import compute_complementarities
from cot_gen.epm import compute_epms, compute_espm, latent_join
from oracles.epm_oracle import epm_oracle, espm_oracle

from tests.gold_networks import ALL_GOLD

_here = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
_BIOMD91 = os.path.join(_repo, "data", "biomodels", "BioMD_other", "BIOMD0000000091.txt")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run_pipeline(rn):
    ercs = compute_ercs(rn)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn_result=syn)
    epm = compute_epms(rn, ercs, hier, syn_result=syn, comp_result=comp, verbose=False)
    espm = compute_espm(rn, ercs, hier, syn, comp, epm, max_order=30, verbose=False)
    return ercs, hier, syn, comp, epm, espm


def _assert_matches_oracle_modulo_latent_joins(rn, ercs, got_epm, got_all):
    """
    got_epm/got_all must be sound (⊆ oracle) and every oracle SO absent from
    got_all must be reconstructible by repeatedly joining pairs already in
    got_all via latent_join() — i.e. the only things compute_espm omits are
    pure latent joins, never a genuine generated-core member.
    """
    oracle_epm = set(epm_oracle(rn, ercs))
    oracle_all = set(espm_oracle(rn, ercs)["all_sos"])

    assert got_epm == oracle_epm, (
        f"EPM mismatch: ours-only={got_epm - oracle_epm}  oracle-only={oracle_epm - got_epm}"
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
        f"compute_espm omitted genuine generated-core members "
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

_WORKED_NAMES = ["a", "b", "c", "d", "e", "f", "g", "h", "k", "p", "q", "w", "z"]
_WORKED_IDX = {n: i for i, n in enumerate(_WORKED_NAMES)}


def _wmask(*species):
    m = 0
    for s in species:
        m |= 1 << _WORKED_IDX[s]
    return m


def _worked_rndata() -> RNData:
    # r_A: f+a->2a+p   r_B: p+b->2b+f   r_C: p+c->2c+f
    # r_D: g+d->2d+q   r_E: q+e->2e+g   r_F: p+h->2h
    # r_H: p+z->z+b    r_G: k+w->2k+p
    reactions = [
        ("r_A", _wmask("f", "a"), _wmask("a", "p")),
        ("r_B", _wmask("p", "b"), _wmask("b", "f")),
        ("r_C", _wmask("p", "c"), _wmask("c", "f")),
        ("r_D", _wmask("g", "d"), _wmask("d", "q")),
        ("r_E", _wmask("q", "e"), _wmask("e", "g")),
        ("r_F", _wmask("p", "h"), _wmask("h")),
        ("r_H", _wmask("p", "z"), _wmask("z", "b")),
        ("r_G", _wmask("k", "w"), _wmask("k", "p")),
    ]
    supp = [r[1] for r in reactions]
    prod = [r[2] for r in reactions]
    names = [r[0] for r in reactions]
    n_sp = len(_WORKED_NAMES)
    inv = tuple(
        tuple(r for r, s in enumerate(supp) if (s >> i) & 1)
        for i in range(n_sp)
    )
    return RNData(
        n_species=n_sp,
        species_names=tuple(_WORKED_NAMES),
        species_index=tuple((n, i) for i, n in enumerate(_WORKED_NAMES)),
        n_reactions=len(reactions),
        reaction_names=tuple(names),
        supp_raw=tuple(supp), prod_raw=tuple(prod),
        E0_mask=0, supp_q=tuple(supp), prod_q=tuple(prod),
        species_to_reactions=inv,
    )


def _names_of(rn, mask):
    return sorted(rn.species_name(i) for i in range(rn.n_species) if (mask >> i) & 1)


@pytest.fixture(scope="module")
def worked():
    rn = _worked_rndata()
    ercs, hier, syn, comp, epm, espm = _run_pipeline(rn)
    return rn, ercs, hier, syn, comp, epm, espm


def test_worked_example_epms_exact(worked):
    rn, ercs, hier, syn, comp, epm, espm = worked
    got = {frozenset(_names_of(rn, m)) for m in epm.all_epm_masks}
    expected = {
        frozenset(["a", "b", "f", "p"]),  # Y1 = E_A join E_B
        frozenset(["a", "c", "f", "p"]),  # Y2 = E_A join E_C
        frozenset(["d", "e", "g", "q"]),  # Y3 = E_D join E_E
    }
    assert got == expected


def test_worked_example_espms_match_oracle_exactly(worked):
    """This network has no pure latent joins, so the generated core equals
    all of Sos exactly (no reconciliation needed)."""
    rn, ercs, hier, syn, comp, epm, espm = worked
    got_all = set(espm.all_so_masks)
    oracle_all = set(espm_oracle(rn, ercs)["all_sos"])
    assert got_all == oracle_all


def test_worked_example_vertical_lift_M3_found(worked):
    """M3 = Y1 join E' (vertical lift absorbing E_B).  Before the Mode-2 fix,
    compute_espm found zero ESPMs on this network because M3/M1/M2 all
    require a non-synergy extension move."""
    rn, ercs, hier, syn, comp, epm, espm = worked
    m3 = _wmask("a", "b", "f", "p", "z")
    assert m3 in set(espm.all_so_masks)


def test_worked_example_discovery_M1_M2_found(worked):
    """M1 = Y1 join E*, M2 = Y2 join E* (complementarity-consumer growth)."""
    rn, ercs, hier, syn, comp, epm, espm = worked
    m1 = _wmask("a", "b", "f", "h", "p")
    m2 = _wmask("a", "c", "f", "h", "p")
    got_all = set(espm.all_so_masks)
    assert m1 in got_all
    assert m2 in got_all


def test_worked_example_sterile_erc_excluded(worked):
    """E_dagger = clos{k,w} requires w externally, which nothing produces —
    its species must never appear in any EPM or ESPM."""
    rn, ercs, hier, syn, comp, epm, espm = worked
    w_bit = _WORKED_IDX["w"]
    k_bit = _WORKED_IDX["k"]
    for m in espm.all_so_masks:
        assert not ((m >> w_bit) & 1), "sterile species w leaked into a persistent module"
        assert not ((m >> k_bit) & 1), "sterile ERC's species k leaked into a persistent module"


def test_worked_example_orders(worked):
    rn, ercs, hier, syn, comp, epm, espm = worked

    def order_of(*names):
        m = _wmask(*names)
        if m in set(epm.all_epm_masks):
            return 0
        for k, masks in espm.espm_by_order.items():
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
# Oracle cross-validation — gold networks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("net", ALL_GOLD, ids=[n.name for n in ALL_GOLD])
def test_epm_espm_oracle_vs_opt_gold(net):
    rn = _build_rndata_from_gold(net)
    ercs = compute_ercs(rn)
    if not ercs:
        pytest.skip(f"{net.name} has no ERCs")
    _, _, _, _, epm, espm = _run_pipeline(rn)
    # compute_epms additionally reports E0 itself as an extra EPM whenever
    # E0_mask != 0 (erc.py's module docstring) -- a deliberate departure
    # from Def 29's original "E_∅ = ∅" convention, which epm_oracle/
    # espm_oracle still encode literally (E0-alone fails their own
    # _is_reactive check, by construction, since inflow reactions have
    # supp_q=0). Excluded here rather than taught to the oracle, same as
    # test_erc.py's test_erc_oracle_vs_opt_gold.
    got_epm = {m for m in epm.all_epm_masks if m != rn.E0_mask}
    got_all = {m for m in espm.all_so_masks if m != rn.E0_mask}
    _assert_matches_oracle_modulo_latent_joins(rn, ercs, got_epm, got_all)


# ---------------------------------------------------------------------------
# Oracle cross-validation — a real biological network
# ---------------------------------------------------------------------------

def test_epm_espm_oracle_vs_opt_biomd91():
    """
    BIOMD0000000091 (12 ERCs, oracle-feasible): the network on which the
    Mode-2 fix was originally caught (compute_espm was silently returning
    fewer ESPMs than exist, and briefly — mid-fix — returning invalid ones,
    before the species-level synergy-ignition fix in erc_syn_close).
    """
    if not os.path.exists(_BIOMD91):
        pytest.skip(f"BIOMD0000000091 not found at {_BIOMD91}")

    from pyCOT.io.functions import read_txt
    from cot_gen.io_pyCOT import build_rndata

    rn_pycot = read_txt(_BIOMD91)
    rn = build_rndata(rn_pycot, network_id="BIOMD0000000091")
    ercs, hier, syn, comp, epm, espm = _run_pipeline(rn)

    assert len(ercs) <= 24, "keep this fixture small enough for the brute-force oracle"
    _assert_matches_oracle_modulo_latent_joins(
        rn, ercs, set(epm.all_epm_masks), set(espm.all_so_masks)
    )
