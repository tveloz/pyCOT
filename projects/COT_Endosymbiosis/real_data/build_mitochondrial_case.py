"""
build_mitochondrial_case.py -- real-data "mitochondrial-type" endosymbiosis
merger (THEORY.md worked plan, item 3).

Host: e_coli_core with its three explicitly oxidative-phosphorylation
reactions removed (ATPS4r fwd+rev = ATP synthase; CYTBD = cytochrome oxidase
bd, the O2-consuming terminal step; NADH16 = NADH dehydrogenase / Complex-I
equivalent, which feeds electrons from NADH into the quinone pool). What
remains is a real, self-consistent glycolysis + TCA + fermentation network
that can extract only substrate-level-phosphorylation ATP -- a structural
stand-in for a fermentative/microaerophilic host lineage, built by deletion
from a real curated model rather than invented from scratch.

Symbiont: the full, unmodified e_coli_core (same reconstruction), standing
in structurally for a free-living aerobic respirer -- not a literal
alphaproteobacterium (BiGG has none in the catalogue used by this project),
stated explicitly as a modelling simplification per THEORY.md's honesty
notes. Reusing the same organism's own reconstruction for both host and
symbiont additionally means there is no cross-organism species-identity
ambiguity to resolve -- every "atp_c" in both files means the same real
metabolite, so the relabelling step is doing exactly and only the job
THEORY.md assigns it (keeping two physically distinct compartments/pools
distinct), not silently reconciling different organisms' namespacing
conventions too.

Interface (curated, stated explicitly): the fermentative host exports
pyruvate to the symbiont (its main fermentation substrate, now available for
full oxidation instead), and the symbiont exports ATP back -- the textbook
mitochondrial endosymbiosis trade, built the same way as the toy model but
now with the real ~50-reaction e_coli_core glycolysis/TCA network on both
sides instead of two lumped reactions.
"""
from __future__ import annotations

import os
import sys

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_here, '..', 'scripts'))
from merge_networks import merge_networks
sys.path.insert(0, os.path.join(_repo_root, 'projects',
                                 'COT_Fundamental_Generators_Exploration', 'scripts'))
from inflow_regime_analysis import build_scenario_txt

ECOLI_CORE = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG', 'bigg_e_coli_core.txt')
CORE_FOOD = ['co2_e', 'glc__D_e', 'h2o_e', 'h_e', 'nh4_e', 'o2_e', 'pi_e']

# Reaction names to delete to build the fermentative-only host (oxidative
# phosphorylation reactions, identified by their own BiGG annotation -- see
# module docstring).
OXPHOS_REACTIONS = {'R34_fwd', 'R35_rev', 'R50', 'R134'}


def build_fermentative_host(out_path: str) -> str:
    """Fermentative host: e_coli_core, oxidative-phosphorylation reactions
    deleted, native inflows stripped and replaced with the standard 7-species
    aerobic core food set (same convention as the companion paper), so this
    file is directly comparable to that paper's own e_coli_core numbers."""
    with open(ECOLI_CORE, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    kept = []
    removed = []
    for line in lines:
        body = line.split(';', 1)[0]
        name = body.split(':', 1)[0].strip() if ':' in body else None
        if name in OXPHOS_REACTIONS:
            removed.append(line.strip())
            continue
        kept.append(line)
    assert len(removed) == len(OXPHOS_REACTIONS), \
        f"expected to remove {len(OXPHOS_REACTIONS)} reactions, removed {len(removed)}: {removed}"
    txt = build_scenario_txt("".join(kept), CORE_FOOD)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(txt)
    print(f"[build_fermentative_host] removed {len(removed)} oxphos reactions:")
    for r in removed:
        print(f"    {r}")
    print(f"[build_fermentative_host] wrote {out_path}")
    return out_path


def build_full_symbiont(out_path: str) -> str:
    """The unmodified e_coli_core (still with the same standard 7-species
    aerobic core food set applied, for the same reason)."""
    with open(ECOLI_CORE, 'r', encoding='utf-8') as f:
        raw = f.read()
    txt = build_scenario_txt(raw, CORE_FOOD)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(txt)
    return out_path


INTERFACE_REACTIONS = [
    # Host pyruvate -> symbiont (its main fermentation substrate, now fully
    # oxidisable by the symbiont's intact oxidative phosphorylation).
    "XPORT_PYR: 1 pyr_c => 1 pyr_c__endo;",
    # Symbiont ATP -> host (the payoff: far more ATP per unit carbon than
    # substrate-level phosphorylation alone).
    "XPORT_ATP: 1 atp_c__endo => 1 atp_c;",
]


def main():
    host_path = build_fermentative_host(os.path.join(_here, 'mito_host_fermentative.txt'))
    symb_path = build_full_symbiont(os.path.join(_here, 'mito_symbiont_aerobic.txt'))
    merged_path = merge_networks(
        host_path=host_path,
        symbiont_path=symb_path,
        interface_reactions=INTERFACE_REACTIONS,
        out_path=os.path.join(_here, 'mito_merged.txt'),
    )
    print(f"[build_mitochondrial_case] wrote {merged_path}")


if __name__ == "__main__":
    main()
