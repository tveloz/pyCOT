"""
benchmarks/stage1_biomd237.py — Stage 1 metrics on Biomodel 237.

Run from the repo root:
    python benchmarks/stage1_biomd237.py

Outputs:
  • Stage-0 network statistics
  • Stage-1 ERC table (one row per ERC)
  • StageContext summary (wall time, peak memory, counters)
"""
from __future__ import annotations

import os
import sys

# Layout: pyCOT/projects/COT_Fundamental_Generators_Complex/benchmarks/
# _proj = .../COT_Fundamental_Generators_Complex/  (cot_gen + oracles importable from here)
# _repo = .../pyCOT/                       (pyCOT importable from _repo/src)
_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from pyCOT.analysis.organizations.io_pyCOT import load_rndata
from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.metrics import StageContext, Counters, print_summary

_BIOMD237 = os.path.join(
    _repo, "data", "biomodels", "BioMD_other", "BIOMD0000000237.txt"
)


def main():
    # ---- Stage 0 -----------------------------------------------------------
    print("=" * 70)
    print("STAGE 0 — Load & E0 quotient")
    print("=" * 70)

    with StageContext("stage0", "BIOMD0000000237", append_log=True) as ctx0:
        rn = load_rndata(_BIOMD237, network_id="BIOMD0000000237")

    ctx0.row.n_species = rn.n_species
    ctx0.row.n_reactions = rn.n_reactions

    n_E0 = bin(rn.E0_mask).count('1')
    n_nontrivial = sum(1 for s in rn.supp_q if s != 0)

    print(f"  species:            {rn.n_species}")
    print(f"  reactions:          {rn.n_reactions}")
    print(f"  E0 species:         {n_E0}  ({rn.bitset_to_names(rn.E0_mask)})")
    print(f"  non-trivial rxns:   {n_nontrivial}  (supp_q != 0)")
    print(f"  wall:               {ctx0.row.wall_s*1000:.2f} ms")
    print(f"  peak:               {ctx0.row.peak_mb:.3f} MB")

    # ---- Stage 1 -----------------------------------------------------------
    print()
    print("=" * 70)
    print("STAGE 1 — ERC discovery")
    print("=" * 70)

    ctr = Counters()
    with StageContext(
        "stage1", "BIOMD0000000237",
        n_species=rn.n_species, n_reactions=rn.n_reactions,
        counters=ctr, append_log=True,
    ) as ctx1:
        ercs = compute_ercs(rn, counters=ctr)

    ctx1.row.n_ercs = len(ercs)

    print(f"  ERCs found:         {len(ercs)}")
    print(f"  persistent ERCs:    {sum(1 for e in ercs if e.is_persistent())}")
    print(f"  wall:               {ctx1.row.wall_s*1000:.2f} ms")
    print(f"  peak:               {ctx1.row.peak_mb:.3f} MB")
    print(f"  closures computed:  {ctr.get('erc.closures_computed', 0)}")
    print(f"  unique ERCs:        {ctr.get('erc.unique_ercs', 0)}")
    print(f"  total min_bases:    {ctr.get('erc.min_bases_total', 0)}")
    print()

    # ---- ERC table ---------------------------------------------------------
    print(f"{'id':>3}  {'|E|':>5}  {'|R_E|':>6}  {'|MinBas|':>8}  "
          f"{'persistent':>10}  {'req species'}")
    print("-" * 65)
    for e in sorted(ercs, key=lambda x: x.size()):
        req_names = rn.bitset_to_names(e.req_mask) if e.req_mask else []
        print(
            f"{e.erc_id:>3}  {e.size():>5}  {len(e.reaction_indices):>6}  "
            f"{len(e.min_bases):>8}  {str(e.is_persistent()):>10}  "
            f"{', '.join(req_names[:3])}{'...' if len(req_names)>3 else ''}"
        )

    # ---- Metrics summary ---------------------------------------------------
    print()
    print("=" * 70)
    print("METRICS LOG")
    print("=" * 70)
    print_summary()


if __name__ == "__main__":
    main()
