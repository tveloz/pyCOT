"""
build_chloroplast_case.py -- real-data "chloroplast-type" endosymbiosis
merger (THEORY.md worked plan, item 3), genuinely cross-organism this time.

Host: Saccharomyces cerevisiae (iND750) -- a real heterotrophic eukaryote,
no photosynthetic capacity at all.

Symbiont: Synechocystis sp. PCC 6803 (iJN678) -- a real cyanobacterium,
oxygenic photosynthesis, already used (as a free-living organism) in the
companion paper's cross-organism EPM comparison. Directly mirrors the real
evolutionary origin of the plastid: a photosynthetic cyanobacterium engulfed
by a eukaryotic host.

Both organisms use standard BiGG cytoplasm/extracellular compartment tags
("_c"/"_e"), so unlike the mitochondrial-type case (same reconstruction
twice), this is the first real test of merge_networks.py's relabelling on
two genuinely independently-curated reconstructions with overlapping but
NOT identical reaction-ID and species-ID conventions.

Interface (curated, stated explicitly): the real, textbook mechanism is the
chloroplast triose-phosphate translocator (TPT) -- fixed carbon leaves the
plastid as a triose phosphate (glyceraldehyde-3-phosphate, g3p, or
dihydroxyacetone phosphate, dhap), not as free glucose; we use g3p, the
canonical TPT substrate. The reverse direction -- host CO2 back into the
symbiont's Calvin cycle, recycling respired carbon -- is the other
textbook-standard half of the exchange.

    XPORT_G3P: symbiont g3p -> host g3p   (fixed carbon out, real TPT route)
    XPORT_CO2: host co2 -> symbiont co2   (respired carbon back in)
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

YEAST = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG', 'bigg_iND750.txt')
SYNECHOCYSTIS = os.path.join(_repo_root, 'data', 'biochemical_databases', 'BiGG', 'bigg_iJN678.txt')

# Matched core used throughout this project and the companion paper. Note
# Synechocystis's own native medium normally includes photon_e (light) --
# deliberately NOT included in this matched-core regime, so that in the
# free-living "symbiont alone" baseline it is heterotrophic/dark and cannot
# fix carbon at all; the point of the merger is to test whether the
# INTERFACE alone (not photon supply) drives complexification. A second
# variant with photon_e included is run separately for comparison.
CORE_FOOD = ['co2_e', 'glc__D_e', 'h2o_e', 'h_e', 'nh4_e', 'o2_e', 'pi_e']
CORE_FOOD_WITH_LIGHT = CORE_FOOD + ['photon_e']

INTERFACE_REACTIONS = [
    "XPORT_G3P: 1 g3p_c__endo => 1 g3p_c;",
    "XPORT_CO2: 1 co2_c => 1 co2_c__endo;",
]


def build_variant(label: str, symbiont_food: list[str]):
    out_dir = _here
    with open(YEAST, 'r', encoding='utf-8') as f:
        host_txt = build_scenario_txt(f.read(), CORE_FOOD)
    with open(SYNECHOCYSTIS, 'r', encoding='utf-8') as f:
        symb_txt = build_scenario_txt(f.read(), symbiont_food)

    host_path = os.path.join(out_dir, f'chloro_host_yeast.txt')
    symb_path = os.path.join(out_dir, f'chloro_symbiont_synecho_{label}.txt')
    with open(host_path, 'w', encoding='utf-8') as f:
        f.write(host_txt)
    with open(symb_path, 'w', encoding='utf-8') as f:
        f.write(symb_txt)

    merged_path = merge_networks(
        host_path=host_path,
        symbiont_path=symb_path,
        interface_reactions=INTERFACE_REACTIONS,
        out_path=os.path.join(out_dir, f'chloro_merged_{label}.txt'),
    )
    print(f"[{label}] host={host_path}\n  symbiont={symb_path}\n  merged={merged_path}")
    return host_path, symb_path, merged_path


def main():
    # Variant A: dark/heterotrophic symbiont baseline (no light) -- isolates
    # the interface's own effect from the (much larger, expected) effect of
    # simply switching the symbiont's own regime to photoautotrophic.
    build_variant('dark', CORE_FOOD)
    # Variant B: symbiont supplied with light too -- the biologically
    # complete picture (a real chloroplast/cyanobiont DOES get light).
    build_variant('light', CORE_FOOD_WITH_LIGHT)


if __name__ == "__main__":
    main()
