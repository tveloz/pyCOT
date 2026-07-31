"""
generate_heteropolymer_network.py -- build a binary-string heteropolymer
network (raf/heteropolymer.py), save it as a pyCOT .txt file, and run it
through the SAME RAF + single-set-decomposition pipeline used elsewhere in
this project (raf/biomodel_crs.py's net-zero catalyst induction,
raf/raf_algo.py's compute_maxRAF, raf/maxraf_decomp.py's
analyze_maxraf_decomposition) -- no new analysis code, just a new network
source.

Edit the CONFIGURATION block and run:
  python projects/RAF_Comparison/scripts/generate_heteropolymer_network.py
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_proj = os.path.normpath(os.path.join(_here, ".."))
_decomp_proj = os.path.normpath(os.path.join(_here, "..", "..", "Decomposition_Theorem"))
_cot_gen_proj = os.path.normpath(os.path.join(_here, "..", "..", "COT_fundamental_Generators"))
_repo = os.path.normpath(os.path.join(_here, "..", "..", ".."))
for _p in (_proj, _decomp_proj, _cot_gen_proj, os.path.join(_repo, "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8", line_buffering=True)

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION -- edit this block, then run                                 ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# LOCKED-IN CONFIG for the paper's Sec. 4.1 network. Full reversibility
# (cleavage_prob=1.0) made every semi-organization trivially an organization
# (any reaction paired with its exact reverse gets a free Sv=0 witness for
# that pair, in any closed set containing it -- regardless of what the
# circuit's own origin is), so the semi-organization/organization gap the
# theory is built to distinguish never showed up. cleavage_prob=0.5 removes
# that universal escape hatch while keeping the necessity of SOME
# length-decreasing reaction (a purely uncatalyzed ligation-only network
# provably cannot support a fragile circuit at all -- but catalysis itself,
# assigned without regard to species length in 'uniform' mode, is a SECOND,
# independent source of circularity: a longer species can gate production
# of a shorter one it depends on, exactly like a RAF catalytic 2-cycle).
MAX_LENGTH = 3
CATALYSIS_MODE = "uniform"      # 'uniform' | 'template' | 'preferential'
P_CATALYST = 0.2                # only used by 'uniform' / 'preferential'
INCLUDE_CLEAVAGE = True
CLEAVAGE_PROB = 0.5             # independent per-ligation coin flip for its reverse
REACTION_KEEP_PROB = 1.0
SEED = 3

OUTPUT_DIR = os.path.join(_proj, "data", "heteropolymer")

# ╚══════════════════════════════════════════════════════════════════════════════╝

from raf.heteropolymer import generate_heteropolymer_network  # noqa: E402

os.makedirs(OUTPUT_DIR, exist_ok=True)
net = generate_heteropolymer_network(
    max_length=MAX_LENGTH,
    catalysis_mode=CATALYSIS_MODE,
    p_catalyst=P_CATALYST,
    include_cleavage=INCLUDE_CLEAVAGE,
    cleavage_prob=CLEAVAGE_PROB,
    reaction_keep_prob=REACTION_KEEP_PROB,
    seed=SEED,
)

out_name = f"heteropolymer_L{MAX_LENGTH}_{CATALYSIS_MODE}_cleave{CLEAVAGE_PROB:g}.txt"
out_path = os.path.join(OUTPUT_DIR, out_name)
with open(out_path, "w", encoding="utf-8") as f:
    f.write(net.txt)

SEP = "=" * 100
print(f"{SEP}\nHeteropolymer network  (max_length={MAX_LENGTH}, catalysis_mode={CATALYSIS_MODE})\n{SEP}")
print(f"  wrote {out_path}")
for k, v in net.stats.items():
    print(f"  {k}: {v}")

# ── Load it back exactly as any other network in this project ─────────────
from pyCOT.io.functions import read_txt  # noqa: E402
from cot_gen.io_pyCOT import build_rndata  # noqa: E402
from raf.biomodel_crs import crs_from_biomodel  # noqa: E402
from raf.maxraf_decomp import analyze_maxraf_decomposition  # noqa: E402

rn_pycot = read_txt(out_path, exact_names=True)
rn_data = build_rndata(rn_pycot, network_id=out_name)
print(f"\n  parsed back: species={rn_data.n_species}  reactions={rn_data.n_reactions}")

crs, crs_stats = crs_from_biomodel(rn_pycot, rn_data)
print(f"  CRS induction (net_zero): {crs_stats}")

a = analyze_maxraf_decomposition(crs, net.food)
r = a["result"]
sizes = sorted((c.size() for c in r.circuits), reverse=True)
print(f"\n  maxRAF = {len(a['maxraf'])} reactions   gen(F0,maxRAF) = {len(a['X_raf'])} species")
print(f"  decomposition of X_raf: |E|={bin(r.E_mask).count('1')}  |F|={bin(r.F_mask).count('1')}  "
      f"circuits={len(r.circuits)}  sizes={sizes}  is_organization={r.is_organization}")
if r.circuits:
    print(f"  circuit species (up to 15 each):")
    for c in r.circuits:
        names = sorted(a["shim"].bitset_to_names(c.species_mask))
        shown = names if len(names) <= 15 else names[:15] + ["..."]
        print(f"    {shown}  self_maintaining={c.is_self_maintaining}")
