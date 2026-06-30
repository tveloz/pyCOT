#!/usr/bin/env python3
"""
classify_biomd_models.py
========================
Fetch metadata from the BioModels REST API for every BIOMD*.txt file in
data/biomodels/biomodels_all_txt/, classify each model into one of:

    metabolic | cell_cycle | signaling | gene_regulation |
    apoptosis | immune | other

based on GO term annotations, then write a classification CSV and
create per-category symlinked subdirectories under data/biomodels/.

Classification priority (first match wins):
    metabolic > cell_cycle > signaling > gene_regulation >
    apoptosis > immune > other

The BiGG models (bigg_*.txt) and central_ecoli.txt are classified
separately without API calls.

Output
------
  projects/Generative_Structure_Orgs/outputs/biomd_classification.csv
  data/biomodels/BioMD_metabolic/         ← symlinks to .txt + .pkl
  data/biomodels/BioMD_cell_cycle/
  data/biomodels/BioMD_signaling/
  data/biomodels/BioMD_gene_regulation/
  data/biomodels/BioMD_apoptosis/
  data/biomodels/BioMD_immune/
  data/biomodels/BioMD_other/
  data/biomodels/BiGG/
  data/biomodels/Other/
"""

import os
import sys
import json
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT  = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
_DATA_DIR    = os.path.join(_PYCOT_ROOT, 'data', 'biomodels', 'biomodels_all_txt')
_BIOMD_DIR   = os.path.join(_PYCOT_ROOT, 'data', 'biomodels')
_OUT_DIR     = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', 'outputs'))
_CSV_OUT     = os.path.join(_OUT_DIR, 'biomd_classification.csv')

os.makedirs(_OUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# GO-term → category mapping
# ---------------------------------------------------------------------------
# Each entry: (go_accession_prefix_or_exact, category)
# Matched left-to-right; first match wins.
GO_CATEGORY_RULES = [
    # Metabolic
    ('GO:0008152', 'metabolic'),   # metabolic process
    ('GO:0044237', 'metabolic'),   # cellular metabolic process
    ('GO:0044281', 'metabolic'),   # small molecule metabolic process
    ('GO:0006096', 'metabolic'),   # glycolytic process
    ('GO:0006098', 'metabolic'),   # pentose-phosphate shunt
    ('GO:0006099', 'metabolic'),   # tricarboxylic acid cycle
    ('GO:0006629', 'metabolic'),   # lipid metabolic process
    ('GO:0006520', 'metabolic'),   # amino acid metabolic process
    ('GO:0009117', 'metabolic'),   # nucleotide metabolic process
    ('GO:0006144', 'metabolic'),   # purine nucleobase metabolic process
    ('GO:0006163', 'metabolic'),   # purine nucleotide metabolic process
    ('GO:0006139', 'metabolic'),   # nucleobase metabolic process
    ('GO:0006749', 'metabolic'),   # glutathione metabolic process
    ('GO:0019320', 'metabolic'),   # hexose catabolic process
    ('GO:0006090', 'metabolic'),   # pyruvate metabolic process
    ('GO:0006091', 'metabolic'),   # generation of precursor metabolites and energy
    ('GO:0015980', 'metabolic'),   # energy derivation by oxidation
    ('GO:0045333', 'metabolic'),   # cellular respiration
    ('GO:0006119', 'metabolic'),   # oxidative phosphorylation
    ('GO:0022904', 'metabolic'),   # respiratory electron transport chain
    ('GO:0006094', 'metabolic'),   # gluconeogenesis
    ('GO:0006071', 'metabolic'),   # glycerol metabolic process
    ('GO:0006066', 'metabolic'),   # alcohol metabolic process
    ('GO:0019362', 'metabolic'),   # pyridine metabolic process
    ('GO:0006733', 'metabolic'),   # oxidoreduction coenzyme metabolic process
    ('GO:0006732', 'metabolic'),   # coenzyme metabolic process
    ('GO:0006725', 'metabolic'),   # cellular aromatic compound metabolic process
    ('GO:0009108', 'metabolic'),   # coenzyme biosynthetic process
    ('GO:0044255', 'metabolic'),   # cellular lipid metabolic process
    ('GO:0005975', 'metabolic'),   # carbohydrate metabolic process
    ('GO:0006040', 'metabolic'),   # amino sugar metabolic process
    ('GO:0006706', 'metabolic'),   # steroid catabolic process
    ('GO:0008202', 'metabolic'),   # steroid metabolic process
    # Cell cycle
    ('GO:0007049', 'cell_cycle'),  # cell cycle
    ('GO:0000278', 'cell_cycle'),  # mitotic cell cycle
    ('GO:0051726', 'cell_cycle'),  # regulation of cell cycle
    ('GO:0022402', 'cell_cycle'),  # cell cycle process
    ('GO:0045787', 'cell_cycle'),  # positive regulation of cell cycle
    ('GO:0045786', 'cell_cycle'),  # negative regulation of cell cycle
    ('GO:0000087', 'cell_cycle'),  # mitotic M phase
    ('GO:0000080', 'cell_cycle'),  # mitotic G1 phase
    # Signaling
    ('GO:0007165', 'signaling'),   # signal transduction
    ('GO:0007166', 'signaling'),   # cell surface receptor signaling pathway
    ('GO:0007267', 'signaling'),   # cell-cell signaling
    ('GO:0023052', 'signaling'),   # signaling
    ('GO:0019226', 'signaling'),   # transmission of nerve impulse
    ('GO:0007274', 'signaling'),   # neuromuscular synaptic transmission
    ('GO:0048583', 'signaling'),   # regulation of response to stimulus
    ('GO:0007154', 'signaling'),   # cell communication
    ('GO:0000165', 'signaling'),   # MAPK cascade
    ('GO:0038061', 'signaling'),   # NIK/NF-kappaB signaling
    ('GO:0016055', 'signaling'),   # Wnt signaling pathway
    ('GO:0007178', 'signaling'),   # transmembrane receptor protein serine/threonine kinase
    ('GO:0007173', 'signaling'),   # epidermal growth factor receptor signaling
    ('GO:0007259', 'signaling'),   # JAK-STAT cascade
    ('GO:0043410', 'signaling'),   # positive regulation of MAPK cascade
    ('GO:0040013', 'signaling'),   # negative regulation of locomotion
    ('GO:0007186', 'signaling'),   # G protein-coupled receptor signaling
    ('GO:0007187', 'signaling'),   # G protein-coupled receptor signaling pathway, coupled to cAMP
    ('GO:0007188', 'signaling'),   # adenylate cyclase-modulating G protein-coupled receptor
    ('GO:0030168', 'signaling'),   # platelet activation
    ('GO:0033674', 'signaling'),   # positive regulation of kinase activity
    ('GO:0045859', 'signaling'),   # regulation of protein kinase activity
    # Gene regulation
    ('GO:0006355', 'gene_regulation'),  # regulation of transcription, DNA-templated
    ('GO:0010468', 'gene_regulation'),  # regulation of gene expression
    ('GO:0006351', 'gene_regulation'),  # transcription, DNA-templated
    ('GO:0006366', 'gene_regulation'),  # transcription from RNA pol II promoter
    ('GO:0045893', 'gene_regulation'),  # positive regulation of transcription
    ('GO:0045892', 'gene_regulation'),  # negative regulation of transcription
    ('GO:0006353', 'gene_regulation'),  # DNA-templated transcription, termination
    ('GO:0043161', 'gene_regulation'),  # proteasome-mediated ubiquitin-dependent protein catabolic
    ('GO:0006605', 'gene_regulation'),  # protein targeting
    ('GO:0000184', 'gene_regulation'),  # nuclear-transcribed mRNA catabolic process
    # Metabolic (additional specific terms)
    ('GO:0019253', 'metabolic'),   # reductive pentose-phosphate cycle (Calvin)
    ('GO:0005986', 'metabolic'),   # sucrose biosynthetic process
    ('GO:0009401', 'metabolic'),   # PTS (phosphotransferase system)
    ('GO:0046655', 'metabolic'),   # folic acid metabolic process
    ('GO:1901575', 'metabolic'),   # organic substance catabolic process
    ('GO:0046364', 'metabolic'),   # monosaccharide biosynthetic process
    ('GO:0000162', 'metabolic'),   # tryptophan biosynthetic process
    ('GO:0006167', 'metabolic'),   # AMP biosynthetic process
    ('GO:0006110', 'metabolic'),   # regulation of glycolysis
    ('GO:0006109', 'metabolic'),   # regulation of carbohydrate metabolic process
    ('GO:0016692', 'metabolic'),   # NADH peroxidase activity
    ('GO:0018205', 'metabolic'),   # peptidyl-lysine modification
    ('GO:0005518', 'metabolic'),   # collagen binding (placeholder, keeps Ferreira here)
    ('GO:0019362', 'metabolic'),   # pyridine metabolic process
    ('GO:0006006', 'metabolic'),   # glucose metabolic process
    ('GO:0006007', 'metabolic'),   # glucose catabolic process
    ('GO:0004022', 'metabolic'),   # alcohol dehydrogenase activity
    ('GO:0006085', 'metabolic'),   # acetyl-CoA biosynthetic process
    ('GO:0009116', 'metabolic'),   # nucleoside metabolic process
    # Circadian / oscillatory
    ('GO:0042752', 'circadian'),   # regulation of circadian rhythm
    ('GO:0032922', 'circadian'),   # circadian regulation of gene expression
    ('GO:0007622', 'circadian'),   # rhythmic behavior
    ('GO:0042753', 'circadian'),   # positive regulation of circadian rhythm
    ('GO:0042754', 'circadian'),   # negative regulation of circadian rhythm
    ('GO:0007623', 'circadian'),   # circadian rhythm
    # Signaling (calcium + electrophysiology additions)
    ('GO:0019722', 'signaling'),   # calcium-mediated signaling
    ('GO:0048016', 'signaling'),   # inositol phosphate-mediated signaling
    ('GO:0050848', 'signaling'),   # regulation of calcium-mediated signaling
    ('GO:0005249', 'signaling'),   # voltage-gated potassium channel activity
    ('GO:0019227', 'signaling'),   # neuronal action potential propagation
    ('GO:0050796', 'signaling'),   # regulation of insulin secretion
    ('GO:0010018', 'signaling'),   # far-red light signaling pathway
    ('GO:0050804', 'signaling'),   # modulation of synaptic transmission
    # Gene regulation (more)
    ('GO:0040029', 'gene_regulation'),  # regulation of gene expression, epigenetic
    ('GO:0051726', 'cell_cycle'),  # regulation of cell cycle (duplicate — already above)
    # Apoptosis
    ('GO:0006915', 'apoptosis'),   # apoptotic process
    ('GO:0012501', 'apoptosis'),   # programmed cell death
    ('GO:0008219', 'apoptosis'),   # cell death
    ('GO:0043065', 'apoptosis'),   # positive regulation of apoptotic process
    ('GO:0043066', 'apoptosis'),   # negative regulation of apoptotic process
    ('GO:0097190', 'apoptosis'),   # apoptotic signaling pathway
    # Immune
    ('GO:0006955', 'immune'),      # immune response
    ('GO:0002376', 'immune'),      # immune system process
    ('GO:0050776', 'immune'),      # regulation of immune response
    ('GO:0006954', 'immune'),      # inflammatory response
    ('GO:0045087', 'immune'),      # innate immune response
    ('GO:0002250', 'immune'),      # adaptive immune response
]

# Name-keyword fallback (applied when GO terms give no match)
NAME_KEYWORD_RULES = [
    # Metabolic
    ('metabol', 'metabolic'),
    ('glycoly', 'metabolic'),
    ('tca cycle', 'metabolic'),
    ('krebs', 'metabolic'),
    ('oxidative phosphorylation', 'metabolic'),
    ('fatty acid', 'metabolic'),
    ('amino acid', 'metabolic'),
    ('nucleotide', 'metabolic'),
    ('calvin', 'metabolic'),
    ('sucrose', 'metabolic'),
    ('purine', 'metabolic'),
    ('folate', 'metabolic'),
    ('glycogen', 'metabolic'),
    ('pentose', 'metabolic'),
    # Cell cycle
    ('cell cycle', 'cell_cycle'),
    ('mitosis', 'cell_cycle'),
    ('cdk', 'cell_cycle'),
    ('cyclin', 'cell_cycle'),
    ('cellcycle', 'cell_cycle'),
    # Circadian
    ('circadian', 'circadian'),
    ('circclock', 'circadian'),
    ('oscillat', 'circadian'),
    # Signaling
    ('signal', 'signaling'),
    ('mapk', 'signaling'),
    ('nfkb', 'signaling'),
    ('wnt', 'signaling'),
    ('egfr', 'signaling'),
    ('receptor', 'signaling'),
    ('kinase', 'signaling'),
    ('camp', 'signaling'),
    ('calcium', 'signaling'),
    ('ca oscillat', 'signaling'),
    ('hodgkin', 'signaling'),
    # Gene regulation
    ('transcri', 'gene_regulation'),
    ('gene expression', 'gene_regulation'),
    ('mrna', 'gene_regulation'),
    ('repressilat', 'gene_regulation'),
    ('operon', 'gene_regulation'),
    # Apoptosis
    ('apoptos', 'apoptosis'),
    ('caspase', 'apoptosis'),
    # Immune
    ('immune', 'immune'),
    ('nk cell', 'immune'),
    ('cytokine', 'immune'),
    ('inflamm', 'immune'),
]

CATEGORY_ORDER = ['metabolic', 'cell_cycle', 'circadian', 'signaling',
                  'gene_regulation', 'apoptosis', 'immune', 'other']

# ---------------------------------------------------------------------------
# GO-term accession → category (fast lookup dict)
# ---------------------------------------------------------------------------
_GO_LOOKUP = {acc: cat for acc, cat in GO_CATEGORY_RULES}


def classify_by_go(go_terms):
    """Return category for a list of (accession, name) GO-term pairs.
    Returns 'other' if no rule matches."""
    for cat in CATEGORY_ORDER[:-1]:   # skip 'other'
        for acc, _ in go_terms:
            if _GO_LOOKUP.get(acc) == cat:
                return cat
    return None  # no GO match — caller will try name fallback


def classify_by_name(name):
    name_lc = name.lower()
    for kw, cat in NAME_KEYWORD_RULES:
        if kw in name_lc:
            return cat
    return 'other'


# ---------------------------------------------------------------------------
# Fetch metadata from BioModels REST API
# ---------------------------------------------------------------------------
def fetch_meta(model_id):
    """Fetch name + GO terms for one BIOMD model ID.
    Returns (model_id, name, go_terms, error_str)."""
    # Strip variant suffix for API call (e.g. BIOMD0000000237_manyOrgs → BIOMD0000000237)
    api_id = model_id.split('_')[0] if '_' in model_id else model_id
    url = f'https://www.ebi.ac.uk/biomodels/{api_id}?format=json'
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=20) as r:
                data = json.loads(r.read())
            name = data.get('name', '')
            go_terms = [
                (ann.get('accession', ''), ann.get('name', ''))
                for ann in data.get('modelLevelAnnotations', [])
                if ann.get('resource') == 'Gene Ontology'
            ]
            return model_id, name, go_terms, None
        except Exception as e:
            if attempt < 2:
                time.sleep(1 + attempt)
            else:
                return model_id, '', [], str(e)


# ---------------------------------------------------------------------------
# Main classification routine
# ---------------------------------------------------------------------------
def main():
    # Collect all BIOMD IDs from the data directory
    files = os.listdir(_DATA_DIR)
    biomd_ids = sorted(set(
        f.replace('.txt', '').replace('.pkl', '')
        for f in files
        if f.startswith('BIOMD')
    ))
    print(f"Found {len(biomd_ids)} BIOMD model IDs to classify.")

    # Check for existing CSV to avoid re-fetching
    existing = {}
    if os.path.exists(_CSV_OUT):
        df_ex = pd.read_csv(_CSV_OUT)
        existing = dict(zip(df_ex['model_id'], df_ex['category']))
        print(f"  Loaded {len(existing)} previously classified models from cache.")

    to_fetch = [bid for bid in biomd_ids if bid not in existing]
    print(f"  Fetching metadata for {len(to_fetch)} models...")

    results = []
    n_done = 0
    with ThreadPoolExecutor(max_workers=12) as ex:
        futs = {ex.submit(fetch_meta, bid): bid for bid in to_fetch}
        for fut in as_completed(futs):
            mid, name, go_terms, err = fut.result()
            if err:
                cat = classify_by_name(mid)   # fall back on ID
                go_str = ''
            else:
                cat = classify_by_go(go_terms) or classify_by_name(name)
                go_str = '; '.join(f"{a}:{n}" for a, n in go_terms)
            results.append({'model_id': mid, 'name': name,
                            'category': cat, 'go_terms': go_str,
                            'api_error': err or ''})
            n_done += 1
            if n_done % 50 == 0:
                print(f"    {n_done}/{len(to_fetch)} fetched...")

    # Merge with existing
    df_new = pd.DataFrame(results)
    if existing:
        df_ex = pd.read_csv(_CSV_OUT)
        df_all = pd.concat([df_ex, df_new], ignore_index=True).drop_duplicates('model_id')
    else:
        df_all = df_new

    df_all.to_csv(_CSV_OUT, index=False)
    print(f"\nClassification saved -> {_CSV_OUT}")

    # Summary
    print("\nCategory distribution:")
    print(df_all['category'].value_counts().to_string())
    return df_all


# ---------------------------------------------------------------------------
# Folder creation and file organization
# ---------------------------------------------------------------------------
CATEGORY_FOLDER = {
    'metabolic':       'BioMD_metabolic',
    'cell_cycle':      'BioMD_cell_cycle',
    'circadian':       'BioMD_circadian',
    'signaling':       'BioMD_signaling',
    'gene_regulation': 'BioMD_gene_regulation',
    'apoptosis':       'BioMD_apoptosis',
    'immune':          'BioMD_immune',
    'other':           'BioMD_other',
}


def organize_folders(df_all):
    """Create subfolders and copy files (symlinks on POSIX, copies on Windows)."""
    import shutil

    # Create all category dirs + BiGG + Other
    all_dirs = list(CATEGORY_FOLDER.values()) + ['BiGG', 'Other']
    for d in all_dirs:
        os.makedirs(os.path.join(_BIOMD_DIR, d), exist_ok=True)

    placed = 0

    # BIOMD models
    cat_map = dict(zip(df_all['model_id'], df_all['category']))
    for mid, cat in cat_map.items():
        folder = CATEGORY_FOLDER.get(cat, 'BioMD_other')
        for ext in ('.txt', '.pkl'):
            src = os.path.join(_DATA_DIR, mid + ext)
            dst = os.path.join(_BIOMD_DIR, folder, mid + ext)
            if os.path.exists(src) and not os.path.exists(dst):
                shutil.copy2(src, dst)
                placed += 1

    # BiGG models
    bigg_files = [f for f in os.listdir(_DATA_DIR) if f.startswith('bigg_')]
    for f in bigg_files:
        src = os.path.join(_DATA_DIR, f)
        dst = os.path.join(_BIOMD_DIR, 'BiGG', f)
        if not os.path.exists(dst):
            shutil.copy2(src, dst)
            placed += 1

    # central_ecoli.txt
    for f in ['central_ecoli.txt']:
        src = os.path.join(_DATA_DIR, f)
        if os.path.exists(src):
            dst = os.path.join(_BIOMD_DIR, 'Other', f)
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
                placed += 1

    print(f"\nOrganized {placed} files into subfolders under {_BIOMD_DIR}")
    # Count per folder
    for d in all_dirs:
        n = len([f for f in os.listdir(os.path.join(_BIOMD_DIR, d))
                 if f.endswith('.txt')])
        print(f"  {d:30s}: {n} networks")


# ---------------------------------------------------------------------------
if __name__ == '__main__':
    df = main()
    organize_folders(df)
