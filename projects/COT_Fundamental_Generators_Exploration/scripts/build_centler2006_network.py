"""
build_centler2006_network.py

Reconstructs the E. coli sugar-metabolism reaction network from:
  Centler, Speroni di Fenizio, Matsumaru, Dittrich (2006)
  "Chemical Organizations in the Central Sugar Metabolism of Escherichia Coli"
  Appendix B (reaction list) + Table 1.1 (species-set labels).

Transcribed directly and mechanically from the paper's Appendix B patterns
(not hand-typed reaction-by-reaction) specifically to minimize transcription
error risk, with the paper's own stated counts (92 species, 168 reactions)
AND its explicit Table 1.1 species partition used as executable checks below
-- both must pass before any output file is written.

Output: one .txt file per scenario (pyCOT reaction-network format), in
  data/Examples_tests/Centler2006_EcoliSugar/centler_{scenario}.txt

Usage (from repo root):
    python projects/COT_fundamental_Generators/scripts/build_centler2006_network.py

See reproduce_centler2006.py for computing the organization hierarchy on
these networks and comparing against the paper's Fig 1.1 / Table 1.1.
"""
import os

_here = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.normpath(os.path.join(_here, '..', '..', '..'))
OUT_DIR = os.path.join(_repo_root, 'data', 'Examples_tests', 'Centler2006_EcoliSugar')
os.makedirs(OUT_DIR, exist_ok=True)

reactions: list[tuple[str, str, str]] = []  # (lhs_terms, rhs_terms, comment)


def add(lhs, rhs, comment=''):
    reactions.append((lhs, rhs, comment))


# ─────────────────────────────────────────────────────────────────────────
# Section 1 (Appendix B.1): 18 genes with identical synthesis/decay pattern
#   RNAP + PromX -> Tscription + PromX + XmRNA
#   XmRNA -> XmRNA + X
#   XmRNA -> (decay)
#   X -> (decay)
# ─────────────────────────────────────────────────────────────────────────
SIMPLE_GENES = [
    'Crp', 'Cya', 'EIIA', 'EIIBC', 'EI', 'Fbp', 'Fda', 'Gap', 'GlcT', 'Glk',
    'GlpR', 'Gpm', 'HPr', 'LacI', 'Pfk', 'Pgi', 'Pyk', 'Tpi',
]
for g in SIMPLE_GENES:
    add(f'RNAP + Prom{g}', f'Tscription + Prom{g} + {g}mRNA', f'transcription of {g}')
    add(f'{g}mRNA', f'{g}mRNA + {g}', f'translation of {g}')
    add(f'{g}mRNA', '', f'decay {g}mRNA')
    add(g, '', f'decay {g}')

# ─────────────────────────────────────────────────────────────────────────
# Section 2 (Appendix B.2): inducible operons lacZY, glpFK, glpD
# ─────────────────────────────────────────────────────────────────────────
add('RNAP + PromLacZY + Allo + Crp + cAMP',
    'Tscription + PromLacZY + LacZYmRNA + Allo + Crp + cAMP',
    'transcription of lacZY (inducer Allo, activator Crp-cAMP)')
add('LacZYmRNA', 'LacZYmRNA1 + LacZ', 'translation: LacZ')
add('LacZYmRNA1', 'LacZYmRNA + LacY', 'translation: LacY (bicistronic)')
add('LacZYmRNA', '', 'decay LacZYmRNA')
add('LacZYmRNA1', '', 'decay LacZYmRNA1')
add('LacZ', '', 'decay LacZ')
add('LacY', '', 'decay LacY')

add('RNAP + PromGlpFK + G3P + Crp + cAMP',
    'Tscription + PromGlpFK + GlpFKmRNA + G3P + Crp + cAMP',
    'transcription of glpFK (inducer G3P, activator Crp-cAMP)')
add('GlpFKmRNA', 'GlpFKmRNA1 + GlpF', 'translation: GlpF')
add('GlpFKmRNA1', 'GlpFKmRNA + GlpK', 'translation: GlpK (bicistronic)')
add('GlpFKmRNA', '', 'decay GlpFKmRNA')
add('GlpFKmRNA1', '', 'decay GlpFKmRNA1')
add('GlpF', '', 'decay GlpF')
add('GlpK', '', 'decay GlpK')

add('RNAP + PromGlpD + G3P + Crp + cAMP',
    'Tscription + PromGlpD + GlpDmRNA + G3P + Crp + cAMP',
    'transcription of glpD (inducer G3P, activator Crp-cAMP)')
add('GlpDmRNA', 'GlpDmRNA + GlpD', 'translation: GlpD (monocistronic)')
add('GlpDmRNA', '', 'decay GlpDmRNA')
add('GlpD', '', 'decay GlpD')

# ─────────────────────────────────────────────────────────────────────────
# Section 3 (Appendix B.3): RNAP unbinding
# ─────────────────────────────────────────────────────────────────────────
add('Tscription', 'RNAP', 'RNAP unbinding')

# ─────────────────────────────────────────────────────────────────────────
# Section 4 (Appendix B.4): signal transduction, transport, metabolism
# ─────────────────────────────────────────────────────────────────────────
METABOLIC = [
    ('ATP + Cya', 'cAMP + Cya'),
    ('PEP + EI + HPr', 'Pyr + EI + HPrP'),
    ('Pyr + EI + HPrP', 'PEP + EI + HPr'),
    ('EIIA + HPrP', 'EIIAP + HPr'),
    ('EIIAP + HPr', 'EIIA + HPrP'),
    ('Glcex + EIIAP + EIIBC', 'Glc6P + EIIA + EIIBC'),
    ('Glc + EIIAP + EIIBC', 'Glc6P + EIIA + EIIBC'),
    ('Glcex + GlcT', 'Glc + GlcT'),
    ('Lacex + LacY', 'Lac + LacY'),
    ('Lac + LacZ', 'Allo + LacZ'),
    ('Lac + LacZ', 'Glc + Glc6P + LacZ'),
    ('Allo + LacZ', 'Glc + Glc6P + LacZ'),
    ('Glc + Glk', 'Glc6P + Glk'),
    ('Glc6P + Pgi', 'Fru6P + Pgi'),
    ('Fru6P + Pgi', 'Glc6P + Pgi'),
    ('Fru6P + Fbp', 'FBP + Fbp'),
    ('FBP + Fbp', 'Fru6P + Fbp'),
    ('Fru6P + Pfk', 'FBP + Pfk'),
    ('FBP + Fda', 'T3P + DHAP + Fda'),
    ('T3P + DHAP + Fda', 'FBP + Fda'),
    ('Glyex + GlpF', 'Gly + GlpF'),
    ('Gly + GlpF', 'Glyex + GlpF'),
    ('Gly + GlpK', 'G3P + GlpK'),
    ('G3P + GlpD', 'DHAP + GlpD'),
    ('DHAP + Tpi', 'T3P + Tpi'),
    ('T3P + Tpi', 'DHAP + Tpi'),
    ('T3P + Gap', '3PG + Gap'),
    ('3PG + Gap', 'T3P + Gap'),
    ('3PG + Gpm', 'PEP + Gpm'),
    ('PEP + Gpm', '3PG + Gpm'),
    ('PEP + FBP + Pyk', 'Pyr + FBP + Pyk'),
    ('Pyr', 'Metabolism'),
]
for lhs, rhs in METABOLIC:
    add(lhs, rhs, 'metabolic/signal transduction')

# ─────────────────────────────────────────────────────────────────────────
# Section 5 (Appendix B.5): decay reactions
# ─────────────────────────────────────────────────────────────────────────
DECAY_SPECIES = [
    'ATP', 'ADP', 'AMP', 'cAMP', 'EIIAP', 'HPrP', 'Glc', 'Gly', 'Lac', 'Allo',
    'Glc6P', 'G3P', 'Fru6P', 'FBP', 'DHAP', 'T3P', '3PG', 'PEP', 'Pyr', 'Metabolism',
]
for s in DECAY_SPECIES:
    add(s, '', f'decay {s}')

# ─────────────────────────────────────────────────────────────────────────
# Section 6 (Appendix B.6): unconditional input reactions
# NOTE: notably absent from this list, and from the decay list above --
# Glcex, Lacex, Glyex (the paper's own text: "The remaining species that do
# not decay are: all 21 promoter species, RNAP, Tscription, Glcex, Lacex,
# and Glyex.") -- these three have NO reaction touching them unless their
# specific uptake enzyme is ALSO present, which is exactly what makes them
# behave as "free" extras in the organization analysis (see
# reproduce_centler2006.py).
# ─────────────────────────────────────────────────────────────────────────
INPUT_SPECIES = [
    'ATP', 'ADP', 'AMP', 'RNAP',
    'PromCrp', 'PromCya', 'PromEIIA', 'PromEIIBC', 'PromEI', 'PromFbp', 'PromFda',
    'PromGap', 'PromGlcT', 'PromGlk', 'PromGlpD', 'PromGlpR', 'PromGpm', 'PromHPr',
    'PromLacI', 'PromPfk', 'PromPgi', 'PromPyk', 'PromTpi', 'PromGlpFK', 'PromLacZY',
]
for s in INPUT_SPECIES:
    add('', s, f'unconditional input {s}')

BASE_REACTIONS = list(reactions)

# ═════════════════════════════════════════════════════════════════════════
# Verification against the paper's own stated counts
# ═════════════════════════════════════════════════════════════════════════
def collect_species(rxns):
    sp = set()
    for lhs, rhs, _ in rxns:
        for side in (lhs, rhs):
            for term in side.split('+'):
                term = term.strip()
                if term:
                    sp.add(term)
    return sp


def main():
    all_species = collect_species(BASE_REACTIONS)
    print(f'Base network (no sugar inputs): {len(BASE_REACTIONS)} reactions, {len(all_species)} species')
    assert len(BASE_REACTIONS) == 168, f'expected 168 reactions, got {len(BASE_REACTIONS)}'
    assert len(all_species) == 92, f'expected 92 species, got {len(all_species)}'
    print('MATCHES paper: 92 species, 168 reactions.')

    # Cross-check against Table 1.1's explicit species-set partition
    GENES_ENZYMES = {
        'PromCrp', 'PromCya', 'PromEIIA', 'PromEIIBC', 'PromEI', 'PromFbp', 'PromFda',
        'PromGap', 'PromGlcT', 'PromGlk', 'PromGlpD', 'PromGlpFK', 'PromGlpR', 'PromGpm',
        'PromHPr', 'PromLacI', 'PromLacZY', 'PromPfk', 'PromPgi', 'PromPyk', 'PromTpi',
        'RNAP', 'Tscription', 'CrpmRNA', 'CyamRNA', 'EIIAmRNA', 'EIIBCmRNA', 'EImRNA',
        'FbpmRNA', 'FdamRNA', 'GapmRNA', 'GlcTmRNA', 'GlkmRNA', 'GlpRmRNA', 'GpmmRNA',
        'HPrmRNA', 'LacImRNA', 'PfkmRNA', 'PgimRNA', 'PykmRNA', 'TpimRNA', 'Crp', 'Cya',
        'EIIA', 'EIIBC', 'EI', 'Fbp', 'Fda', 'Gap', 'GlcT', 'Glk', 'GlpR', 'Gpm', 'HPr',
        'LacI', 'Pfk', 'Pgi', 'Pyk', 'Tpi', 'AMP', 'ATP', 'ADP', 'cAMP',
    }
    METABOLITES = {'Glc', 'Glc6P', 'Fru6P', 'FBP', 'DHAP', 'T3P', '3PG', 'PEP', 'Pyr',
                    'Metabolism', 'EIIAP', 'HPrP'}
    LAC_SPECIES = {'Lac', 'Allo', 'LacZYmRNA', 'LacZYmRNA1', 'LacZ', 'LacY'}
    GLY_SPECIES = {'Gly', 'G3P', 'GlpDmRNA', 'GlpFKmRNA', 'GlpFKmRNA1', 'GlpD', 'GlpF', 'GlpK'}
    EX_SPECIES = {'Glcex', 'Lacex', 'Glyex'}

    table_species = GENES_ENZYMES | METABOLITES | LAC_SPECIES | GLY_SPECIES | EX_SPECIES
    print(f'\nTable 1.1 partition total: {len(table_species)} species')
    missing = table_species - all_species
    extra = all_species - table_species
    if missing:
        print('Missing from our network (in Table 1.1 but not found):', missing)
    if extra:
        print('Extra in our network (found but not in Table 1.1):', extra)
    assert table_species == all_species, 'MISMATCH vs Table 1.1 species partition!'
    print('MATCHES Table 1.1 species partition exactly.')

    # ═════════════════════════════════════════════════════════════════════
    # Write pyCOT-format .txt files, one per scenario
    # ═════════════════════════════════════════════════════════════════════
    def write_network(rxns, sugar_inputs, out_path):
        lines = []
        idx = 0
        for lhs, rhs, comment in rxns:
            lhs_fmt = ' + '.join(f'1 {t.strip()}' for t in lhs.split('+') if t.strip())
            rhs_fmt = ' + '.join(f'1 {t.strip()}' for t in rhs.split('+') if t.strip())
            lines.append(f'R{idx}: {lhs_fmt} => {rhs_fmt}; {comment}')
            idx += 1
        for s in sugar_inputs:
            lines.append(f'R{idx}: => 1 {s}; unconditional input {s} (scenario)')
            idx += 1
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
        return idx

    SCENARIOS = {
        'starvation': [],
        'glucose':    ['Glcex'],
        'lactose':    ['Lacex'],
        'glycerol':   ['Glyex'],
        'all_sugars': ['Glcex', 'Lacex', 'Glyex'],
    }

    for name, sugars in SCENARIOS.items():
        out_path = os.path.join(OUT_DIR, f'centler_{name}.txt')
        n = write_network(BASE_REACTIONS, sugars, out_path)
        print(f'  wrote {out_path}  ({n} reactions)')

    print('\nDone.')


if __name__ == '__main__':
    main()
