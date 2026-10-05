"""
benchmark_suite.py — network corpus for the Conjecture 1 vs Conjecture 2
comparison (conjectures/compare_strategies.py).

Two tiers:

  small_networks()  — the 5 gold networks (tests/gold_networks.py) plus the
    companion-paper worked example (worked_example.py). Small enough (<=13
    species) for the brute-force oracle (oracles/so_oracle.py), so these
    give an exact correctness anchor, not just a scaling data point.

  bio_networks(max_reactions=100, ...) — a filtered scan of the real
    BioModels + BiGG corpus already on disk under data/biochemical_databases/
    (the same dataset the companion paper's own statistics are drawn from),
    capped to networks with <= max_reactions reactions to give a real
    size-scaling spread up to the paper's own "regular desktop computer"
    target. A cheap line-count pre-filter (count "=>" occurrences) avoids
    paying a full pyCOT parse for every oversized network in the corpus
    before even checking its size.
"""
from __future__ import annotations

import os
import glob
from dataclasses import dataclass

from pyCOT.analysis.organizations.cot_types import RNData

from tests.gold_networks import ALL_GOLD
from worked_example import worked_rndata

_here = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
_DATA_ROOT = os.path.join(_repo, "data", "biochemical_databases")

# Same category folders Generative_Structure_Orgs/Cleaned/config.py scans —
# referenced, not imported, since that project's SCAN_FOLDERS also encodes
# plotting/grouping concerns this module doesn't need.
_CATEGORY_FOLDERS = [
    "BioMD_metabolic", "BioMD_cell_cycle", "BioMD_circadian", "BioMD_signaling",
    "BioMD_gene_regulation", "BioMD_apoptosis", "BioMD_immune", "BioMD_other",
    "BiGG", "Other",
]


@dataclass
class BenchmarkNetwork:
    name: str
    rn_data: RNData
    n_reactions: int
    n_species: int
    category: str   # "gold" | "worked_example" | folder name for bio networks


def _build_rndata_from_gold(net) -> RNData:
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


def small_networks() -> list[BenchmarkNetwork]:
    """5 gold networks + the companion-paper worked example — oracle-checkable."""
    out = []
    for net in ALL_GOLD:
        rn = _build_rndata_from_gold(net)
        out.append(BenchmarkNetwork(
            name=f"gold_{net.name}", rn_data=rn,
            n_reactions=rn.n_reactions, n_species=rn.n_species, category="gold",
        ))
    wrn = worked_rndata()
    out.append(BenchmarkNetwork(
        name="worked_example", rn_data=wrn,
        n_reactions=wrn.n_reactions, n_species=wrn.n_species, category="worked_example",
    ))
    return out


def _cheap_reaction_count(path: str) -> int:
    """Count reaction lines without a full pyCOT parse: one non-empty,
    non-comment line per reaction in this text format, each containing '=>'."""
    n = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            if "=>" in line:
                n += 1
    return n


def discover_bio_networks(max_reactions: int = 100, min_reactions: int = 4) -> list[tuple[str, str, int]]:
    """
    Scan data/biochemical_databases/{category}/*.txt, deduplicated by
    basename (many networks are cross-listed, e.g. BiGG networks also
    appear under a generic 'Other'/'biomodels_all_txt' export), filtered to
    [min_reactions, max_reactions] by the cheap line-count pre-filter.

    Returns [(name, path, n_reactions_estimate), ...] sorted by size ascending.
    """
    seen: dict[str, tuple[str, int]] = {}
    for cat in _CATEGORY_FOLDERS:
        folder = os.path.join(_DATA_ROOT, cat)
        if not os.path.isdir(folder):
            continue
        for path in sorted(glob.glob(os.path.join(folder, "*.txt"))):
            name = os.path.splitext(os.path.basename(path))[0]
            if name.endswith("_manyOrgs"):
                continue   # variant files, not independent networks
            if name.startswith("bigg_"):
                name = name[len("bigg_"):]
            if name in seen:
                continue
            n_rxn = _cheap_reaction_count(path)
            if min_reactions <= n_rxn <= max_reactions:
                seen[name] = (path, n_rxn)
    return sorted(((n, p, r) for n, (p, r) in seen.items()), key=lambda t: t[2])


def bio_networks(
    max_reactions: int = 100,
    min_reactions: int = 4,
    sample_cap: int | None = 50,
) -> list[BenchmarkNetwork]:
    """
    Load up to `sample_cap` real BioModels/BiGG networks with
    min_reactions <= n_reactions <= max_reactions, evenly spread across the
    size range (not just the smallest ones) so the comparison covers the
    full requested scale, not only its bottom end.
    """
    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations.io_pyCOT import build_rndata

    candidates = discover_bio_networks(max_reactions=max_reactions, min_reactions=min_reactions)
    if sample_cap is not None and len(candidates) > sample_cap:
        step = len(candidates) / sample_cap
        candidates = [candidates[int(i * step)] for i in range(sample_cap)]

    out = []
    for name, path, _n_est in candidates:
        try:
            rn_pycot = read_txt(path)
            rn = build_rndata(rn_pycot, network_id=name)
        except Exception:
            continue
        category = os.path.basename(os.path.dirname(path))
        out.append(BenchmarkNetwork(
            name=name, rn_data=rn,
            n_reactions=rn.n_reactions, n_species=rn.n_species, category=category,
        ))
    out.sort(key=lambda b: b.n_reactions)
    return out


def full_suite(max_reactions: int = 100, sample_cap: int | None = 50) -> list[BenchmarkNetwork]:
    return small_networks() + bio_networks(max_reactions=max_reactions, sample_cap=sample_cap)
