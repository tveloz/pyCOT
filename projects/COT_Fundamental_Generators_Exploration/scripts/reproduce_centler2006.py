"""
reproduce_centler2006.py

Computes the FULL organization hierarchy (closed + self-maintaining species
sets, in the classical Dittrich/Speroni di Fenizio sense -- NOT the
ERC/SO0/SOi engine in pyCOT.analysis.organizations) for all 5 scenarios of:

  Centler, Speroni di Fenizio, Matsumaru, Dittrich (2006)
  "Chemical Organizations in the Central Sugar Metabolism of Escherichia Coli"

and checks the result against the paper's own Fig 1.1 / Table 1.1 (exactly
4 organizations per scenario, specific species-set sizes).

Uses src/pyCOT/analysis/Persistent_Modules_Generator.py's
brute_force_organizations -- the exhaustive, provably-complete LP-based
enumerator (as opposed to compute_all_organizations, which combines
elementary organizations pairwise only and is not guaranteed complete).

──────────────────────────────────────────────────────────────────────────
WHY THIS ISN'T JUST brute_force_organizations() ALONE
──────────────────────────────────────────────────────────────────────────
A first pass using brute_force_organizations() alone found only 1, 1, 2, 2,
and 4 organizations for the five scenarios (expected: 4 in every scenario).
Root cause, precisely diagnosed (not assumed) by direct testing:

  Glcex, Lacex, Glyex are the only three species in this network with NO
  decay reaction and NO reaction triggered by their presence alone -- the
  paper's own text confirms this exactly: "The remaining species that do
  not decay are: all 21 promoter species, RNAP, Tscription, Glcex, Lacex,
  and Glyex." Any one of them can therefore be freely added to an already-
  valid organization without disturbing closure (nothing new gets
  triggered) or self-maintenance (their own net production is trivially
  zero, since no reaction touches them). This is exactly the paper's own
  explanation for why its Org.2/Org.3 exist: "they are isolated nodes in
  the reaction network ... fulfilling closure and self-maintenance."

  brute_force_organizations() enumerates closures reachable as "closure of
  a union of ERC closures". In this network EVERY ERC's closure already
  spans 64+ of the 92 species (confirmed by direct inspection), so no ERC
  isolates "just Glyex" cleanly -- the free-species-addition states are
  real organizations but are not reachable via ERC-combination alone.

  Fix: after the brute-force pass, for every organization found, test
  adding every subset of {Glcex, Lacex, Glyex} not already contained in
  it, and verify EACH candidate for real via the same LP self-maintenance
  check used everywhere else (never assumed just because the reasoning
  above sounds right).

This same principle (species with no decay reaction and no self-triggered
reaction are freely addable) generalizes beyond this one paper -- watch for
it whenever brute_force_organizations() is applied to a new network with a
similar structure (e.g. an un-metabolized alternate carbon/energy source
sitting inert in a genome-scale model).

──────────────────────────────────────────────────────────────────────────
RUNTIME
──────────────────────────────────────────────────────────────────────────
Not cheap -- the brute-force ERC-combination search is exponential in the
number of ERCs. Observed wall times (this machine): all_sugars ~3s,
glucose ~1min, lactose/glycerol ~6min, starvation ~52min (23 ERCs is
already enough to make this slow). This does NOT scale to genome-scale
networks like iAF692 (hundreds of ERCs) -- see network_structure_profile.py
and the project notes on what a genome-scale-appropriate approach would
need to look like.

Usage (from repo root):
    python projects/COT_Fundamental_Generators_Exploration/scripts/build_centler2006_network.py   # once
    python projects/COT_Fundamental_Generators_Exploration/scripts/reproduce_centler2006.py
"""
import sys, os, time, json, itertools

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
sys.path.insert(0, os.path.join(_repo_root, 'src'))

from pyCOT.io.functions import read_txt
from pyCOT.analysis.Persistent_Modules_Generator import brute_force_organizations, check_self_maintenance
from pyCOT.analysis.SORN_Generators import is_semi_self_maintaining
from pyCOT.analysis.ERC_Hierarchy import closure

NET_DIR = os.path.join(_repo_root, 'data', 'Examples_tests', 'Centler2006_EcoliSugar')
OUT_DIR = os.path.join(_repo_root, 'projects', 'COT_Fundamental_Generators_Exploration', 'outputs', 'centler2006_reproduction')
os.makedirs(OUT_DIR, exist_ok=True)

SCENARIOS = ['starvation', 'glucose', 'lactose', 'glycerol', 'all_sugars']
FREE_CANDIDATES = ['Glcex', 'Lacex', 'Glyex']

# Expected sizes per organization, per scenario, derived from Centler et al.
# Table 1.1 (Genes+Enzymes=63, Metabolites=12, LacSpecies=6, GlySpecies=8,
# Metabolites*=Metabolites-{Glc}=11), sorted ascending -- used as the
# executable pass/fail check below.
EXPECTED_SIZES = {
    'starvation': [63, 64, 64, 65],
    'glucose':    [76, 77, 77, 78],
    'lactose':    [64, 65, 82, 83],
    'glycerol':   [64, 65, 83, 84],
    'all_sugars': [78, 84, 86, 92],
}


def extend_with_free_species(orgs, RN):
    """For each found organization, test adding every subset of the
    not-yet-present free-candidate species, verifying via real LP checks
    (see module docstring for why this is needed and why it's valid)."""
    extended = set(orgs)
    for base in list(orgs):
        missing = [s for s in FREE_CANDIDATES if s not in base]
        for r in range(1, len(missing) + 1):
            for combo in itertools.combinations(missing, r):
                candidate_names = set(base) | set(combo)
                sp_objects = [sp for sp in RN.species() if sp.name in candidate_names]
                cl = closure(RN, sp_objects)
                cl_frozen = frozenset(sp.name for sp in cl)
                if cl_frozen in extended:
                    continue
                # must still equal exactly base+combo -- otherwise closure
                # pulled in something else and this isn't a "free" addition
                if cl_frozen != frozenset(candidate_names):
                    continue
                cl_sp = [sp for sp in RN.species() if sp.name in cl_frozen]
                if not is_semi_self_maintaining(RN, cl_sp):
                    continue
                if check_self_maintenance(cl_sp, RN)[0]:
                    extended.add(cl_frozen)
    return extended


def main():
    all_results = {}
    for name in SCENARIOS:
        path = os.path.join(NET_DIR, f'centler_{name}.txt')
        print(f'\n{"#"*70}\n# SCENARIO: {name}\n{"#"*70}')
        RN = read_txt(path, exact_names=True)

        t0 = time.perf_counter()
        res = brute_force_organizations(RN, max_combo_size=None, early_stop=True, verbose=False)
        base_orgs = res['organizations']
        print(f'  brute-force (early_stop) found: {len(base_orgs)} organizations')

        extended = extend_with_free_species(base_orgs, RN)
        wall = time.perf_counter() - t0
        print(f'  after free-species extension:   {len(extended)} organizations  ({wall:.1f}s)')

        orgs_sorted = sorted(extended, key=lambda s: (len(s), sorted(s)))
        all_results[name] = {
            'wall_s': wall,
            'n_organizations': len(orgs_sorted),
            'organizations': [sorted(o) for o in orgs_sorted],
        }
        for i, org in enumerate(orgs_sorted):
            print(f'    Org {i+1}: {len(org)} species')

    with open(os.path.join(OUT_DIR, 'results.json'), 'w') as f:
        json.dump(all_results, f, indent=2)

    print(f'\n{"="*70}')
    print('SUMMARY vs. Centler et al. (2006) Fig 1.1 / Table 1.1')
    print(f'{"="*70}')
    all_pass = True
    report_lines = ['# Centler et al. (2006) reproduction — results\n']
    for name in SCENARIOS:
        got_sizes = sorted(len(o) for o in all_results[name]['organizations'])
        expected = EXPECTED_SIZES[name]
        ok = got_sizes == expected
        all_pass &= ok
        line = (f'  {name:12s}: sizes {got_sizes} vs expected {expected}  '
                f'-> {"MATCH" if ok else "MISMATCH"}')
        print(line)
        report_lines.append(line)
    print(f'\n{"ALL SCENARIOS MATCH" if all_pass else "SOME SCENARIOS DID NOT MATCH"}')
    report_lines.append(f'\n{"ALL SCENARIOS MATCH" if all_pass else "SOME SCENARIOS DID NOT MATCH"}')

    with open(os.path.join(OUT_DIR, 'summary.md'), 'w') as f:
        f.write('\n'.join(report_lines) + '\n')
    print(f'\nSaved -> {os.path.join(OUT_DIR, "results.json")}')
    print(f'Saved -> {os.path.join(OUT_DIR, "summary.md")}')


if __name__ == '__main__':
    main()
