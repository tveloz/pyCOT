#!/usr/bin/env python3
"""
script_erc_synergy_stats.py
============================
Batch-process a folder of reaction network .txt files, compute ERC synergy
statistics for each, save the results to a CSV, and plot scatter plots showing
how synergy counts scale with network size.

Statistics collected per network:
  n_species, n_reactions, n_ercs,
  n_basic, n_maximal, n_fundamental,
  ratio_basic    = n_basic    / C(n_ercs,2)   (fraction of pairs with basic synergy)
  ratio_maximal  = n_maximal  / C(n_ercs,2)
  ratio_fundamental = n_fundamental / C(n_ercs,2)
  time_ercs, time_synergies

Synergy detection uses the corrected implementations from script_erc_synergy_viz.py
(the library's get_maximal_synergies has an inverted containment check; see
the docstring in script_erc_synergy_viz.py for details).

Default folders scanned:
  1. data/biomodels/biomodels_interesting/
  2. networks/testing/performance_benchmark/
  3. networks/testing/   (other .txt files)
Add or change SCAN_FOLDERS below.

Networks that fail (parse error, timeout, too many ERCs) are skipped and logged.
"""

import os
import sys
import time
import traceback
from itertools import combinations
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import networkx as nx

# -- Path setup ----------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))
sys.path.insert(0, _SCRIPT_DIR)

from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from utils_ercs import load_ercs

# -- Configuration -------------------------------------------------------------
_BIOMD = os.path.join(_PYCOT_ROOT, 'data', 'biomodels')
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

OUT_DIR  = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', 'outputs', 'synergy_stats'))
CSV_FILE = os.path.join(OUT_DIR, 'synergy_stats.csv')

MAX_REACTIONS  = 1000 # skip networks with more reactions than this (pre-ERC filter)
MAX_ERCS       = 900  # skip networks with more than this many ERCs (post-ERC filter)
MIN_ERCS       = 4    # skip networks with fewer than this many ERCs (post-ERC filter)
MAX_ERCS_3SYN  = 900   # only compute ternary synergies for networks this size or smaller
MAX_TIME_ERCS  = 600  # max seconds for ERC computation per network

os.makedirs(OUT_DIR, exist_ok=True)


# ===========================================================================
# Corrected synergy functions (same logic as script_erc_synergy_viz.py)
# ===========================================================================

def _compute_basic(erc1, erc2, hierarchy, RN):
    """Basic synergies: (erc1,erc2)→T for each jointly-coverable-but-not-individually generator."""
    if (erc1 in hierarchy.get_contain(erc2) or
            erc2 in hierarchy.get_contain(erc1)):
        return []
    cl1   = erc1.get_closure_names(RN)
    cl2   = erc2.get_closure_names(RN)
    joint = cl1 | cl2
    contained_by_1 = {e.label for e in hierarchy.get_contain(erc1)}
    contained_by_2 = {e.label for e in hierarchy.get_contain(erc2)}

    result = []
    for target in hierarchy.ercs:
        if target is erc1 or target is erc2:
            continue
        if target.label in contained_by_1 or target.label in contained_by_2:
            continue
        for gen in target.min_generators:
            gen_sp = set(species_list_to_names(gen))
            if not gen_sp.issubset(joint):
                continue
            if gen_sp.issubset(cl1) or gen_sp.issubset(cl2):
                continue
            result.append((erc1, erc2, target))
            break
    return result


def _filter_maximal(basics, RN):
    """
    Keep only maximal synergies (Definition 19): for each pair (E1,E2), keep
    targets T such that no other target T' in the basic synergies for that pair
    has cl(T') ⊋ cl(T).
    """
    if not basics:
        return []
    by_pair = defaultdict(list)
    for e1, e2, target in basics:
        key = tuple(sorted([e1.label, e2.label]))
        by_pair[key].append((e1, e2, target))

    result = []
    for syns in by_pair.values():
        closures = {s[2].label: s[2].get_closure_names(RN) for s in syns}
        for e1, e2, T in syns:
            cl_T = closures[T.label]
            dominated = any(
                cl_T < closures[T2.label]     # T strictly smaller than T2
                for _, _, T2 in syns
                if T2.label != T.label
            )
            if not dominated:
                result.append((e1, e2, T))
    return result


def _filter_fundamental(maximals, hierarchy, RN):
    """
    Keep only fundamental synergies (Definition 20): maximal synergy (E1,E2)→T
    is fundamental iff no strictly smaller pair (E1'⊆E1, E2'⊆E2, at least one
    strict) has a maximal synergy to T.
    """
    if not maximals:
        return []

    label_to_erc = {erc.label: erc for erc in hierarchy.ercs}

    def _desc_plus_self(erc):
        if erc.label not in hierarchy.graph:
            return {erc.label}
        return {erc.label} | set(nx.descendants(hierarchy.graph, erc.label))

    desc_cache = {}
    for e1, e2, _ in maximals:
        for e in (e1, e2):
            if e.label not in desc_cache:
                desc_cache[e.label] = _desc_plus_self(e)

    result = []
    for e1, e2, target in maximals:
        is_fund = True
        d1, d2  = desc_cache[e1.label], desc_cache[e2.label]

        for l1 in d1:
            if not is_fund:
                break
            for l2 in d2:
                if l1 == l2:
                    continue
                orig = {tuple(sorted([e1.label, e2.label]))}
                if tuple(sorted([l1, l2])) in orig:
                    continue
                sub1 = label_to_erc.get(l1)
                sub2 = label_to_erc.get(l2)
                if sub1 is None or sub2 is None:
                    continue
                sub_bas = _compute_basic(sub1, sub2, hierarchy, RN)
                sub_max = _filter_maximal(sub_bas, RN)
                if any(s[2].label == target.label for s in sub_max):
                    is_fund = False
                    break

        if is_fund:
            result.append((e1, e2, target))
    return result


# ===========================================================================
# Ternary (degree-3) synergy functions
# A triple (E1,E2,E3) is a basic 3-synergy if some reaction r is triggered by
# the triple but NOT by any of the three pairwise joins.
# ===========================================================================

def _compute_basic_3(erc1, erc2, erc3, hierarchy, RN):
    """Basic ternary synergies: triples that jointly trigger a reaction not
    triggered by any singleton or pair sub-collection."""
    # Containment between any two members makes the triple reducible
    for a, b in [(erc1, erc2), (erc1, erc3), (erc2, erc3)]:
        if (a in hierarchy.get_contain(b) or b in hierarchy.get_contain(a)):
            return []

    cl1 = erc1.get_closure_names(RN)
    cl2 = erc2.get_closure_names(RN)
    cl3 = erc3.get_closure_names(RN)
    triple  = cl1 | cl2 | cl3
    pair12  = cl1 | cl2
    pair13  = cl1 | cl3
    pair23  = cl2 | cl3

    contained_by_any = ({e.label for e in hierarchy.get_contain(erc1)} |
                        {e.label for e in hierarchy.get_contain(erc2)} |
                        {e.label for e in hierarchy.get_contain(erc3)})

    result = []
    for target in hierarchy.ercs:
        if target is erc1 or target is erc2 or target is erc3:
            continue
        if target.label in contained_by_any:
            continue
        for gen in target.min_generators:
            gen_sp = set(species_list_to_names(gen))
            if not gen_sp.issubset(triple):
                continue
            # Must NOT be coverable by any pair (that would be a 2-synergy)
            if (gen_sp.issubset(pair12) or
                gen_sp.issubset(pair13) or
                gen_sp.issubset(pair23)):
                continue
            result.append((erc1, erc2, erc3, target))
            break
    return result


def _filter_maximal_3(basics_3, RN):
    """Keep only maximal ternary synergies for each triple."""
    if not basics_3:
        return []
    by_triple = defaultdict(list)
    for e1, e2, e3, target in basics_3:
        key = tuple(sorted([e1.label, e2.label, e3.label]))
        by_triple[key].append((e1, e2, e3, target))

    result = []
    for syns in by_triple.values():
        closures = {s[3].label: s[3].get_closure_names(RN) for s in syns}
        for e1, e2, e3, T in syns:
            cl_T = closures[T.label]
            dominated = any(
                cl_T < closures[T2.label]
                for _, _, _, T2 in syns if T2.label != T.label
            )
            if not dominated:
                result.append((e1, e2, e3, T))
    return result


def _filter_fundamental_3(maximals_3, hierarchy, RN):
    """Keep only fundamental ternary synergies.
    A maximal 3-synergy (E1,E2,E3)→T is fundamental iff no strictly smaller
    triple (replacing any ERC by a contained one) OR any pair of descendants
    yields the same maximal synergy to T.
    """
    if not maximals_3:
        return []

    label_to_erc = {erc.label: erc for erc in hierarchy.ercs}

    def _desc_plus_self(erc):
        if erc.label not in hierarchy.graph:
            return {erc.label}
        return {erc.label} | set(nx.descendants(hierarchy.graph, erc.label))

    desc_cache = {}
    for e1, e2, e3, _ in maximals_3:
        for e in (e1, e2, e3):
            if e.label not in desc_cache:
                desc_cache[e.label] = _desc_plus_self(e)

    result = []
    for e1, e2, e3, target in maximals_3:
        is_fund = True
        d1 = desc_cache[e1.label]
        d2 = desc_cache[e2.label]
        d3 = desc_cache[e3.label]
        orig_key = tuple(sorted([e1.label, e2.label, e3.label]))

        # Check: does any strictly smaller triple achieve the same maximal target?
        for l1 in d1:
            if not is_fund: break
            for l2 in d2:
                if not is_fund: break
                for l3 in d3:
                    sub_key = tuple(sorted([l1, l2, l3]))
                    if sub_key == orig_key:
                        continue
                    sub1 = label_to_erc.get(l1)
                    sub2 = label_to_erc.get(l2)
                    sub3 = label_to_erc.get(l3)
                    if None in (sub1, sub2, sub3):
                        continue
                    sub_bas = _compute_basic_3(sub1, sub2, sub3, hierarchy, RN)
                    sub_max = _filter_maximal_3(sub_bas, RN)
                    if any(s[3].label == target.label for s in sub_max):
                        is_fund = False
                        break

        # Check: does any pair of descendants achieve the same target as a 2-synergy?
        if is_fund:
            for da, db in [(d1, d2), (d1, d3), (d2, d3)]:
                if not is_fund: break
                for la in da:
                    if not is_fund: break
                    for lb in db:
                        if la == lb:
                            continue
                        sa = label_to_erc.get(la)
                        sb = label_to_erc.get(lb)
                        if None in (sa, sb):
                            continue
                        bas2 = _compute_basic(sa, sb, hierarchy, RN)
                        max2 = _filter_maximal(bas2, RN)
                        if any(s[2].label == target.label for s in max2):
                            is_fund = False
                            break

        if is_fund:
            result.append((e1, e2, e3, target))
    return result


def compute_synergy_stats(ercs, hierarchy, RN):
    """
    Returns a dict with binary (degree-2) and ternary (degree-3) synergy counts.

    Binary:
      n_basic_pairs / n_maximal_pairs / n_fundamental_pairs  — unique (E1,E2) pairs
      n_basic_targets / n_maximal_targets / n_fundamental_targets — unique target ERCs

    Ternary (only if len(ercs) <= MAX_ERCS_3SYN, else -1):
      n_basic_3syn_triples / n_maximal_3syn_triples / n_fundamental_3syn_triples
      n_basic_3syn_targets / n_maximal_3syn_targets / n_fundamental_3syn_targets
      syn3_computed — bool flag
    """
    b_pairs, m_pairs, f_pairs       = set(), set(), set()
    b_targets, m_targets, f_targets = set(), set(), set()

    for e1, e2 in combinations(ercs, 2):
        b = _compute_basic(e1, e2, hierarchy, RN)
        if not b:
            continue
        pair_key = tuple(sorted([e1.label, e2.label]))
        b_pairs.add(pair_key)
        b_targets.update(t.label for _, _, t in b)

        m = _filter_maximal(b, RN)
        if m:
            m_pairs.add(pair_key)
            m_targets.update(t.label for _, _, t in m)

        f = _filter_fundamental(m, hierarchy, RN)
        if f:
            f_pairs.add(pair_key)
            f_targets.update(t.label for _, _, t in f)

    stats = {
        'n_basic_pairs':         len(b_pairs),
        'n_maximal_pairs':       len(m_pairs),
        'n_fundamental_pairs':   len(f_pairs),
        'n_basic_targets':       len(b_targets),
        'n_maximal_targets':     len(m_targets),
        'n_fundamental_targets': len(f_targets),
    }

    # Ternary synergies — only for small enough networks
    if len(ercs) <= MAX_ERCS_3SYN:
        b3_triples, m3_triples, f3_triples = set(), set(), set()
        b3_targets, m3_targets, f3_targets = set(), set(), set()

        for e1, e2, e3 in combinations(ercs, 3):
            b3 = _compute_basic_3(e1, e2, e3, hierarchy, RN)
            if not b3:
                continue
            triple_key = tuple(sorted([e1.label, e2.label, e3.label]))
            b3_triples.add(triple_key)
            b3_targets.update(t.label for _, _, _, t in b3)

            m3 = _filter_maximal_3(b3, RN)
            if m3:
                m3_triples.add(triple_key)
                m3_targets.update(t.label for _, _, _, t in m3)

            f3 = _filter_fundamental_3(m3, hierarchy, RN)
            if f3:
                f3_triples.add(triple_key)
                f3_targets.update(t.label for _, _, _, t in f3)

        stats.update({
            'n_basic_3syn_triples':       len(b3_triples),
            'n_maximal_3syn_triples':     len(m3_triples),
            'n_fundamental_3syn_triples': len(f3_triples),
            'n_basic_3syn_targets':       len(b3_targets),
            'n_maximal_3syn_targets':     len(m3_targets),
            'n_fundamental_3syn_targets': len(f3_targets),
            'syn3_computed':              True,
        })
    else:
        stats.update({
            'n_basic_3syn_triples': -1, 'n_maximal_3syn_triples': -1,
            'n_fundamental_3syn_triples': -1, 'n_basic_3syn_targets': -1,
            'n_maximal_3syn_targets': -1, 'n_fundamental_3syn_targets': -1,
            'syn3_computed': False,
        })

    return stats


# ===========================================================================
# File discovery
# ===========================================================================

def collect_files(folders):
    """folders is a dict {folder_path: dataset_label}.
    Returns list of (abs_path, dataset_label) tuples, deduplicated by real path."""
    seen, files = set(), []
    for folder, label in folders.items():
        if not os.path.isdir(folder):
            continue
        for fname in sorted(os.listdir(folder)):
            if not fname.endswith('.txt'):
                continue
            path = os.path.join(folder, fname)
            real = os.path.realpath(path)
            if real not in seen:
                seen.add(real)
                files.append((path, label))
    return files


# ===========================================================================
# Main loop
# ===========================================================================

all_files = collect_files(SCAN_FOLDERS)
print(f"Found {len(all_files)} .txt files across {len(SCAN_FOLDERS)} folders.")

# Load existing CSV so already-computed networks are not reprocessed.
# Networks where ternary synergy was previously skipped (syn3_computed=False)
# but now fall within MAX_ERCS_3SYN are excluded from the cache so they get
# reprocessed to fill in the ternary columns.
if os.path.exists(CSV_FILE):
    _existing = pd.read_csv(CSV_FILE)
    _existing = _existing[_existing['n_ercs'] >= MIN_ERCS]
    if 'syn3_computed' in _existing.columns:
        needs_ternary = (
            (~_existing['syn3_computed'].astype(bool)) &
            (_existing['n_ercs'] <= MAX_ERCS_3SYN)
        )
        n_requeue = int(needs_ternary.sum())
        if n_requeue:
            print(f"Re-queuing {n_requeue} networks where ternary was previously skipped "
                  f"but n_ercs <= {MAX_ERCS_3SYN} now.")
        _existing = _existing[~needs_ternary]
    already_done = set(_existing['file'].astype(str))
    print(f"Resuming: {len(already_done)} networks already cached in {CSV_FILE}")
else:
    _existing = pd.DataFrame()
    already_done = set()

records = []
skipped = []

for idx, (fpath, dataset_label) in enumerate(all_files):
    fname = os.path.basename(fpath)
    if fname in already_done:
        print(f"[{idx+1}/{len(all_files)}] {fname}  — CACHED, skip")
        continue
    print(f"\n[{idx+1}/{len(all_files)}] {fname}")

    try:
        RN = read_txt(fpath)
        n_sp = len(RN.species())
        n_rx = len(RN.reactions())
        print(f"  {n_sp} species, {n_rx} reactions")

        if n_rx > MAX_REACTIONS:
            print(f"  SKIP: too many reactions ({n_rx} > {MAX_REACTIONS})")
            skipped.append((fname, f'too many reactions: {n_rx}'))
            continue

        t0 = time.time()
        ercs, _cached = load_ercs(fpath, RN, ERC)
        t_ercs = time.time() - t0

        # Filter out E_∅ (empty-closure ERC).
        ercs = [e for e in ercs if len(e.get_closure_names(RN)) > 0]
        n_ercs = len(ercs)
        cache_tag = ' [cache]' if _cached else ''
        print(f"  {n_ercs} ERCs (E_∅ excluded)  ({t_ercs:.1f}s){cache_tag}")

        if n_ercs > MAX_ERCS:
            print(f"  SKIP: too many ERCs ({n_ercs} > {MAX_ERCS})")
            skipped.append((fname, f'too many ERCs: {n_ercs}'))
            continue
        if n_ercs < MIN_ERCS:
            print(f"  SKIP: too few ERCs ({n_ercs} < {MIN_ERCS})")
            skipped.append((fname, f'too few ERCs: {n_ercs}'))
            continue
        if t_ercs > MAX_TIME_ERCS:
            print(f"  SKIP: ERC computation too slow ({t_ercs:.1f}s)")
            skipped.append((fname, f'ERC timeout: {t_ercs:.1f}s'))
            continue

        hierarchy = ERC_Hierarchy(RN, ercs)

        t0  = time.time()
        syn = compute_synergy_stats(ercs, hierarchy, RN)
        t_syn = time.time() - t0
        print(f"  Pairs   — Basic={syn['n_basic_pairs']}  "
              f"Maximal={syn['n_maximal_pairs']}  "
              f"Fundamental={syn['n_fundamental_pairs']}  ({t_syn:.1f}s)")
        print(f"  Targets — Basic={syn['n_basic_targets']}  "
              f"Maximal={syn['n_maximal_targets']}  "
              f"Fundamental={syn['n_fundamental_targets']}")
        if syn['syn3_computed']:
            print(f"  3-Syn triples — Basic={syn['n_basic_3syn_triples']}  "
                  f"Maximal={syn['n_maximal_3syn_triples']}  "
                  f"Fundamental={syn['n_fundamental_3syn_triples']}")
        else:
            print(f"  3-Syn — skipped (n_ercs={n_ercs} > MAX_ERCS_3SYN={MAX_ERCS_3SYN})")

        n_pairs   = n_ercs * (n_ercs - 1) // 2          # C(n,2)
        n_triples = n_ercs * (n_ercs - 1) * (n_ercs - 2) // 6  # C(n,3)
        records.append({
            'file':         fname,
            'dataset':      dataset_label,
            'n_species':    n_sp,
            'n_reactions':  n_rx,
            'n_ercs':       n_ercs,
            'n_pairs_max':   n_pairs,
            'n_triples_max': n_triples,
            # ── Binary (degree-2) synergy ──────────────────────────────────────
            'n_basic_pairs':         syn['n_basic_pairs'],
            'n_maximal_pairs':       syn['n_maximal_pairs'],
            'n_fundamental_pairs':   syn['n_fundamental_pairs'],
            'n_basic_targets':       syn['n_basic_targets'],
            'n_maximal_targets':     syn['n_maximal_targets'],
            'n_fundamental_targets': syn['n_fundamental_targets'],
            'ratio_basic':       syn['n_basic_pairs']       / n_pairs if n_pairs > 0 else 0,
            'ratio_maximal':     syn['n_maximal_pairs']     / n_pairs if n_pairs > 0 else 0,
            'ratio_fundamental': syn['n_fundamental_pairs'] / n_pairs if n_pairs > 0 else 0,
            'target_ratio_basic':       syn['n_basic_targets']       / n_ercs if n_ercs > 0 else 0,
            'target_ratio_maximal':     syn['n_maximal_targets']     / n_ercs if n_ercs > 0 else 0,
            'target_ratio_fundamental': syn['n_fundamental_targets'] / n_ercs if n_ercs > 0 else 0,
            # ── Ternary (degree-3) synergy (-1 = not computed) ────────────────
            'syn3_computed':              syn['syn3_computed'],
            'n_basic_3syn_triples':       syn['n_basic_3syn_triples'],
            'n_maximal_3syn_triples':     syn['n_maximal_3syn_triples'],
            'n_fundamental_3syn_triples': syn['n_fundamental_3syn_triples'],
            'n_basic_3syn_targets':       syn['n_basic_3syn_targets'],
            'n_maximal_3syn_targets':     syn['n_maximal_3syn_targets'],
            'n_fundamental_3syn_targets': syn['n_fundamental_3syn_targets'],
            'ratio_3syn_basic':
                syn['n_basic_3syn_triples'] / n_triples
                if (n_triples > 0 and syn['syn3_computed']) else -1,
            'ratio_3syn_fundamental':
                syn['n_fundamental_3syn_triples'] / n_triples
                if (n_triples > 0 and syn['syn3_computed']) else -1,
            # ── Timing ────────────────────────────────────────────────────────
            'time_ercs':      round(t_ercs, 2),
            'time_synergies': round(t_syn,  2),
        })

        # Write after every network so progress survives interruptions.
        _cur = pd.concat([_existing, pd.DataFrame(records)], ignore_index=True) \
               if not _existing.empty else pd.DataFrame(records)
        _cur.to_csv(CSV_FILE, index=False)
        print(f"  → CSV updated ({len(_cur)} total records)")

    except Exception as exc:
        print(f"  ERROR: {exc}")
        skipped.append((fname, str(exc)))

# ===========================================================================
# Save CSV
# ===========================================================================

new_df = pd.DataFrame(records)
if not _existing.empty and not new_df.empty:
    df = pd.concat([_existing, new_df], ignore_index=True)
elif not new_df.empty:
    df = new_df
else:
    df = _existing
df.to_csv(CSV_FILE, index=False)
print(f"\nSaved {len(df)} records to {CSV_FILE} "
      f"({len(new_df)} new, {len(_existing)} previously cached)")
if skipped:
    print(f"Skipped {len(skipped)} networks:")
    for fn, reason in skipped:
        print(f"  {fn}: {reason}")

if df.empty:
    print("No data to plot.")
    sys.exit(0)

# ===========================================================================
# Plots  (two focused panels)
# ===========================================================================

SYN_STYLES = {
    'basic':       {'color': '#A8D8EA', 'label': 'Basic',       'marker': 'o', 'zorder': 2},
    'maximal':     {'color': '#E67E22', 'label': 'Maximal',     'marker': 's', 'zorder': 3},
    'fundamental': {'color': '#27AE60', 'label': 'Fundamental', 'marker': '^', 'zorder': 4},
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# ---------------------------------------------------------------------------
# Plot 1 — Synergic pairs vs theoretical maximum C(n_ercs, 2)
#
# Log-log axes: networks range from a handful of ERCs to ~200, so C(n,2)
# spans ~1 – 20 000.  On a linear scale the data clusters at the bottom;
# log-log spreads the full dynamic range and keeps y = x as a 45° diagonal.
# Networks whose count is 0 for a given type are omitted (log(0) undefined);
# a note in the title flags this.
# ---------------------------------------------------------------------------
x_max = df['n_pairs_max'].values

for key, sty in SYN_STYLES.items():
    y    = df[f'n_{key}_pairs'].values
    mask = (x_max > 0) & (y > 0)          # log scale needs strictly positive values
    ax1.scatter(x_max[mask], y[mask],
                color=sty['color'], marker=sty['marker'],
                s=55, alpha=0.80, edgecolors='white', lw=0.5,
                label=sty['label'], zorder=sty['zorder'])

# Diagonal y = x is a straight 45° line in log-log space
x_pos = x_max[x_max > 0]
if len(x_pos):
    d_lo = x_pos.min() * 0.7
    d_hi = x_pos.max() * 1.5
    ax1.plot([d_lo, d_hi], [d_lo, d_hi], 'k--', lw=0.9, alpha=0.35,
             label='y = x  (all pairs synergic)')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('C(|ERCs|, 2)  —  theoretical max unordered pairs  [log]', fontsize=11)
ax1.set_ylabel('Synergic pairs  (unique (E1, E2) with ≥ 1 synergy)  [log]', fontsize=11)
ax1.set_title(
    'Synergic pairs vs theoretical maximum  (log–log)\n'
    'Distance below the diagonal reflects rarity  ·  zero counts omitted',
    fontsize=11, fontweight='bold')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.25, which='both')

# ---------------------------------------------------------------------------
# Plot 2 — Synergetic target fraction vs ERC pairs reactants fraction
#
# X: ratio_* = n_*_pairs / C(n_ercs,2)  ∈ [0,1]
#    Fraction of all ERC pairs that ARE reactants in ≥1 synergy.
# Y: target_ratio_* = n_*_targets / n_ercs  ∈ [0,1]
#    Fraction of ERCs that appear as a synergy TARGET.
# Node size ∝ n_ercs so the network scale is visible.
# ---------------------------------------------------------------------------
n_ercs_vals = df['n_ercs'].values
# Map n_ercs → scatter size: min 30 pts², max 400 pts²
_n_min, _n_max = n_ercs_vals.min(), n_ercs_vals.max()
_n_range = max(_n_max - _n_min, 1)
node_s = 30 + 370 * (n_ercs_vals - _n_min) / _n_range

for key, sty in SYN_STYLES.items():
    xv = df[f'ratio_{key}'].values
    yv = df[f'target_ratio_{key}'].values
    ax2.scatter(xv, yv, s=node_s, color=sty['color'], marker=sty['marker'],
                alpha=0.75, edgecolors='#333333', lw=0.4,
                label=sty['label'], zorder=sty['zorder'])

ax2.axvline(0.5, color='gray', lw=0.8, ls=':', alpha=0.5)
ax2.axhline(0.5, color='gray', lw=0.8, ls=':', alpha=0.5)
ax2.set_xscale('symlog', linthresh=0.01)
ax2.set_xlim(0, 1.05)
ax2.set_ylim(-0.02, 1.02)
ax2.set_xlabel(
    'ERC pairs reactants fraction  =  synergic pairs / C(|ERCs|, 2)  [symlog]',
    fontsize=11)
ax2.set_ylabel(
    'Synergetic target fraction  =  unique targets / |ERCs|',
    fontsize=11)
ax2.set_title(
    'Synergetic target fraction vs ERC pairs reactants fraction\n'
    'Node size ∝ |ERCs|  ·  x: symlog  ·  y: linear ∈ [0, 1]',
    fontsize=11, fontweight='bold')

# Size legend: show 3 representative n_ercs values
from matplotlib.lines import Line2D as _L2D
_ticks = [_n_min, (_n_min + _n_max) // 2, _n_max]
size_handles = [
    ax2.scatter([], [], s=30 + 370 * (n - _n_min) / _n_range,
                color='#888888', marker='o', alpha=0.7,
                edgecolors='#333333', lw=0.4, label=f'|ERCs| = {n}')
    for n in _ticks
]
# Merge synergy-type legend + size legend
type_handles = [ax2.scatter([], [], s=80, color=sty['color'],
                             marker=sty['marker'], alpha=0.9, label=sty['label'])
                for sty in SYN_STYLES.values()]
ax2.legend(handles=type_handles + size_handles,
           fontsize=8, loc='upper left',
           title='Type  /  scale', title_fontsize=8)
ax2.grid(True, alpha=0.25)

fig.suptitle(
    f'ERC synergy statistics  —  {len(df)} networks  '
    f'(E_∅ excluded,  max_ercs = {MAX_ERCS})',
    fontsize=12, fontweight='bold')
plt.tight_layout()

fig_path = os.path.join(OUT_DIR, 'synergy_stats.png')
fig.savefig(fig_path, dpi=150, bbox_inches='tight')
print(f"Plot saved to {fig_path}")
plt.show()

# ===========================================================================
# Summary table
# ===========================================================================
print("\n=== Summary ===")
cols = ['file', 'n_species', 'n_reactions', 'n_ercs',
        'n_basic_pairs', 'n_maximal_pairs', 'n_fundamental_pairs',
        'ratio_basic', 'ratio_fundamental',
        'n_basic_3syn_triples', 'n_fundamental_3syn_triples']
print(df[[c for c in cols if c in df.columns]].to_string(index=False))
