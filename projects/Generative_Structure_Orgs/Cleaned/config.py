"""
config.py — shared configuration for all Cleaned/ pipeline scripts.

Edit this file to change scan folders, computation thresholds, or plot styles.
All other scripts import from here so changes propagate automatically.
"""

import os

# ── Paths ──────────────────────────────────────────────────────────────────────
_HERE      = os.path.dirname(os.path.abspath(__file__))
PYCOT_ROOT = os.path.normpath(os.path.join(_HERE, '..', '..', '..'))
_BIOMD     = os.path.join(PYCOT_ROOT, 'data', 'biochemical_databases')

OUT_DIR = os.path.join(_HERE, 'outputs')
VIZ_DIR = os.path.join(_HERE, 'visualizations')

# ── Scan folders ───────────────────────────────────────────────────────────────
SCAN_FOLDERS = {
    os.path.join(_BIOMD, 'BioMD_metabolic'):       'BioMD_metabolic',
    os.path.join(_BIOMD, 'BioMD_cell_cycle'):       'BioMD_cell_cycle',
    os.path.join(_BIOMD, 'BioMD_circadian'):        'BioMD_circadian',
    os.path.join(_BIOMD, 'BioMD_signaling'):        'BioMD_signaling',
    os.path.join(_BIOMD, 'BioMD_gene_regulation'):  'BioMD_gene_regulation',
    os.path.join(_BIOMD, 'BioMD_apoptosis'):        'BioMD_apoptosis',
    os.path.join(_BIOMD, 'BioMD_immune'):           'BioMD_immune',
    os.path.join(_BIOMD, 'BioMD_other'):            'BioMD_other',
    os.path.join(_BIOMD, 'BiGG'):                   'BiGG',
    os.path.join(_BIOMD, 'Other'):                  'Other',
}

# ── Computation thresholds ─────────────────────────────────────────────────────
MAX_REACTIONS    = 1000   # skip networks with more reactions (pre-ERC filter)
MAX_ERCS         = 900    # skip networks with more ERCs (post-ERC filter)
MAX_ERCS_TERNARY = 900    # ternary synergies only for networks <= this many ERCs
MIN_ERCS         = 4      # skip networks with fewer ERCs
MIN_TREE_NODES   = 3      # skip ERC hierarchy trees with fewer nodes

# ── Group assignment ───────────────────────────────────────────────────────────
def assign_group(dataset):
    if str(dataset).startswith('BioMD_'):
        return 'BioModels'
    if dataset == 'BiGG':
        return 'BiGG'
    return 'Other'

# ── Group plot styles ──────────────────────────────────────────────────────────
GROUP_STYLE = {
    'BioModels': {'color': '#E74C3C', 'marker': 'o', 'label': 'BioModels', 'lw': 0.3},
    'BiGG':      {'color': '#2C3E50', 'marker': '*', 'label': 'BiGG',      'lw': 1.5},
    'Other':     {'color': '#7F8C8D', 'marker': 's', 'label': 'Other',     'lw': 0.3},
}

# ── Font and figure sizes (hierarchy + growth plots) ──────────────────────────
FONT_SCALE = 2.0    # multiply all base font sizes by this factor
FBASE_LAB  = 12     # axis label base fontsize (pt, before scaling)
FBASE_TIT  = 11     # title base fontsize
FBASE_LEG  = 9      # legend base fontsize
FBASE_TCK  = 10     # tick label base fontsize
FIG_W      = 14     # figure width (inches)
FIG_H      = 10     # figure height (inches)
