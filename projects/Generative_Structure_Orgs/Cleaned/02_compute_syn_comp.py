#!/usr/bin/env python3
"""
02_compute_syn_comp.py
=======================
Batch-compute ERC synergy and complementarity statistics for every reaction
network in SCAN_FOLDERS.  Both are computed in a single pass per network
(ERC loading and hierarchy building happen once per network).

Outputs (incremental — already-computed networks are skipped):
  outputs/synergy_stats.csv
  outputs/complementarity_stats.csv

Key configuration (edit in config.py):
  MAX_REACTIONS    — networks with more reactions are skipped
  MAX_ERCS         — networks with more ERCs are skipped
  MAX_ERCS_TERNARY — ternary synergies only for networks with <= this many ERCs
  MIN_ERCS         — networks with fewer ERCs are skipped

synergy_stats.csv columns
--------------------------
  file, dataset, n_species, n_reactions, n_ercs, n_pairs_max, n_triples_max,
  n_basic_pairs, n_maximal_pairs, n_fundamental_pairs,
  n_basic_targets, n_maximal_targets, n_fundamental_targets,
  ratio_basic, ratio_maximal, ratio_fundamental,
  target_ratio_basic, target_ratio_maximal, target_ratio_fundamental,
  syn3_computed,
  n_basic_3syn_triples, n_maximal_3syn_triples, n_fundamental_3syn_triples,
  n_basic_3syn_targets, n_maximal_3syn_targets, n_fundamental_3syn_targets,
  ratio_3syn_basic, ratio_3syn_fundamental,
  time_ercs, time_synergies, time_complementarity

complementarity_stats.csv columns
-----------------------------------
  file, dataset, n_species, n_reactions, n_ercs, n_pairs_max,
  n_complementary_pairs, n_pure_complementary_pairs, n_fundamental_edges,
  n_producer_ERCs, n_consumer_ERCs,
  ratio_complementary, ratio_pure, ratio_fundamental,
  ratio_producers, ratio_consumers,
  time_ercs, time_complementarity
"""

import os
import sys
import time
from itertools import combinations
from collections import defaultdict

import pandas as pd
import networkx as nx

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_SCRIPT_DIR, '..', '..', '..', 'src'))
sys.path.insert(0, _SCRIPT_DIR)

from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from utils_ercs import load_ercs
from config import (
    SCAN_FOLDERS, MAX_REACTIONS, MAX_ERCS, MAX_ERCS_TERNARY, MIN_ERCS,
    OUT_DIR,
)

os.makedirs(OUT_DIR, exist_ok=True)

SYN_CSV  = os.path.join(OUT_DIR, 'synergy_stats.csv')
COMP_CSV = os.path.join(OUT_DIR, 'complementarity_stats.csv')

# =============================================================================
# Synergy functions  (binary degree-2)
# =============================================================================

def _compute_basic(erc1, erc2, hierarchy, RN):
    if (erc1 in hierarchy.get_contain(erc2) or
            erc2 in hierarchy.get_contain(erc1)):
        return []
    cl1   = erc1.get_closure_names(RN)
    cl2   = erc2.get_closure_names(RN)
    joint = cl1 | cl2
    sub1  = {e.label for e in hierarchy.get_contain(erc1)}
    sub2  = {e.label for e in hierarchy.get_contain(erc2)}
    result = []
    for target in hierarchy.ercs:
        if target is erc1 or target is erc2:
            continue
        if target.label in sub1 or target.label in sub2:
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
    if not basics:
        return []
    by_pair = defaultdict(list)
    for e1, e2, target in basics:
        by_pair[tuple(sorted([e1.label, e2.label]))].append((e1, e2, target))
    result = []
    for syns in by_pair.values():
        closures = {s[2].label: s[2].get_closure_names(RN) for s in syns}
        for e1, e2, T in syns:
            dominated = any(
                closures[T.label] < closures[T2.label]
                for _, _, T2 in syns if T2.label != T.label
            )
            if not dominated:
                result.append((e1, e2, T))
    return result


def _filter_fundamental(maximals, hierarchy, RN):
    if not maximals:
        return []
    label_to_erc = {e.label: e for e in hierarchy.ercs}

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
        orig    = tuple(sorted([e1.label, e2.label]))
        for l1 in d1:
            if not is_fund:
                break
            for l2 in d2:
                if tuple(sorted([l1, l2])) == orig:
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


# =============================================================================
# Ternary (degree-3) synergy functions
# =============================================================================

def _compute_basic_3(erc1, erc2, erc3, hierarchy, RN):
    for a, b in [(erc1, erc2), (erc1, erc3), (erc2, erc3)]:
        if (a in hierarchy.get_contain(b) or b in hierarchy.get_contain(a)):
            return []
    cl1 = erc1.get_closure_names(RN)
    cl2 = erc2.get_closure_names(RN)
    cl3 = erc3.get_closure_names(RN)
    triple = cl1 | cl2 | cl3
    pair12, pair13, pair23 = cl1 | cl2, cl1 | cl3, cl2 | cl3
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
            if (gen_sp.issubset(pair12) or gen_sp.issubset(pair13) or
                    gen_sp.issubset(pair23)):
                continue
            result.append((erc1, erc2, erc3, target))
            break
    return result


def _filter_maximal_3(basics_3, RN):
    if not basics_3:
        return []
    by_triple = defaultdict(list)
    for e1, e2, e3, target in basics_3:
        by_triple[tuple(sorted([e1.label, e2.label, e3.label]))].append((e1, e2, e3, target))
    result = []
    for syns in by_triple.values():
        closures = {s[3].label: s[3].get_closure_names(RN) for s in syns}
        for e1, e2, e3, T in syns:
            dominated = any(
                closures[T.label] < closures[T2.label]
                for _, _, _, T2 in syns if T2.label != T.label
            )
            if not dominated:
                result.append((e1, e2, e3, T))
    return result


def _filter_fundamental_3(maximals_3, hierarchy, RN):
    if not maximals_3:
        return []
    label_to_erc = {e.label: e for e in hierarchy.ercs}

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
        d1, d2, d3 = desc_cache[e1.label], desc_cache[e2.label], desc_cache[e3.label]
        orig_key = tuple(sorted([e1.label, e2.label, e3.label]))

        for l1 in d1:
            if not is_fund: break
            for l2 in d2:
                if not is_fund: break
                for l3 in d3:
                    sub_key = tuple(sorted([l1, l2, l3]))
                    if sub_key == orig_key:
                        continue
                    s1, s2, s3 = label_to_erc.get(l1), label_to_erc.get(l2), label_to_erc.get(l3)
                    if None in (s1, s2, s3):
                        continue
                    sb = _compute_basic_3(s1, s2, s3, hierarchy, RN)
                    sm = _filter_maximal_3(sb, RN)
                    if any(s[3].label == target.label for s in sm):
                        is_fund = False
                        break

        if is_fund:
            for da, db in [(d1, d2), (d1, d3), (d2, d3)]:
                if not is_fund: break
                for la in da:
                    if not is_fund: break
                    for lb in db:
                        if la == lb: continue
                        sa, sb = label_to_erc.get(la), label_to_erc.get(lb)
                        if None in (sa, sb): continue
                        bas2 = _compute_basic(sa, sb, hierarchy, RN)
                        max2 = _filter_maximal(bas2, RN)
                        if any(s[2].label == target.label for s in max2):
                            is_fund = False
                            break

        if is_fund:
            result.append((e1, e2, e3, target))
    return result


def compute_synergy_stats(ercs, hierarchy, RN):
    b_pairs, m_pairs, f_pairs       = set(), set(), set()
    b_targets, m_targets, f_targets = set(), set(), set()

    for e1, e2 in combinations(ercs, 2):
        b = _compute_basic(e1, e2, hierarchy, RN)
        if not b:
            continue
        key = tuple(sorted([e1.label, e2.label]))
        b_pairs.add(key)
        b_targets.update(t.label for _, _, t in b)
        m = _filter_maximal(b, RN)
        if m:
            m_pairs.add(key)
            m_targets.update(t.label for _, _, t in m)
        f = _filter_fundamental(m, hierarchy, RN)
        if f:
            f_pairs.add(key)
            f_targets.update(t.label for _, _, t in f)

    stats = {
        'n_basic_pairs': len(b_pairs), 'n_maximal_pairs': len(m_pairs),
        'n_fundamental_pairs': len(f_pairs), 'n_basic_targets': len(b_targets),
        'n_maximal_targets': len(m_targets), 'n_fundamental_targets': len(f_targets),
    }

    if len(ercs) <= MAX_ERCS_TERNARY:
        b3, m3, f3     = set(), set(), set()
        bt3, mt3, ft3  = set(), set(), set()
        for e1, e2, e3 in combinations(ercs, 3):
            b3s = _compute_basic_3(e1, e2, e3, hierarchy, RN)
            if not b3s: continue
            key = tuple(sorted([e1.label, e2.label, e3.label]))
            b3.add(key); bt3.update(t.label for _, _, _, t in b3s)
            m3s = _filter_maximal_3(b3s, RN)
            if m3s:
                m3.add(key); mt3.update(t.label for _, _, _, t in m3s)
            f3s = _filter_fundamental_3(m3s, hierarchy, RN)
            if f3s:
                f3.add(key); ft3.update(t.label for _, _, _, t in f3s)
        stats.update({
            'n_basic_3syn_triples': len(b3), 'n_maximal_3syn_triples': len(m3),
            'n_fundamental_3syn_triples': len(f3), 'n_basic_3syn_targets': len(bt3),
            'n_maximal_3syn_targets': len(mt3), 'n_fundamental_3syn_targets': len(ft3),
            'syn3_computed': True,
        })
    else:
        stats.update({
            'n_basic_3syn_triples': -1, 'n_maximal_3syn_triples': -1,
            'n_fundamental_3syn_triples': -1, 'n_basic_3syn_targets': -1,
            'n_maximal_3syn_targets': -1, 'n_fundamental_3syn_targets': -1,
            'syn3_computed': False,
        })
    return stats


# =============================================================================
# Complementarity functions
# =============================================================================

def _is_incomparable(erc1, erc2, hierarchy):
    return (erc1 not in hierarchy.get_contain(erc2) and
            erc2 not in hierarchy.get_contain(erc1))


def _supply(erc_prod, erc_cons, RN):
    return erc_prod.get_produced_species(RN) & erc_cons.get_required_species(RN)


def _compute_minprod_mincons(ercs, hierarchy, RN):
    producers = defaultdict(set)
    consumers = defaultdict(set)
    for erc in ercs:
        for s in erc.get_produced_species(RN):
            producers[s].add(erc.label)
        for s in erc.get_required_species(RN):
            consumers[s].add(erc.label)

    label_to_erc = {e.label: e for e in ercs}

    def _minimal_ercs(label_set):
        minimal = []
        for lbl in label_set:
            erc = label_to_erc[lbl]
            descendants = {e.label for e in hierarchy.get_contain(erc)}
            if not descendants.intersection(label_set):
                minimal.append(erc)
        return minimal

    minprod = {s: _minimal_ercs(lbls) for s, lbls in producers.items()}
    mincons = {s: _minimal_ercs(lbls) for s, lbls in consumers.items()}
    return minprod, mincons


def _has_basic_synergy(erc1, erc2, hierarchy, RN):
    if (erc1 in hierarchy.get_contain(erc2) or
            erc2 in hierarchy.get_contain(erc1)):
        return False
    cl1   = erc1.get_closure_names(RN)
    cl2   = erc2.get_closure_names(RN)
    joint = cl1 | cl2
    sub1  = {e.label for e in hierarchy.get_contain(erc1)}
    sub2  = {e.label for e in hierarchy.get_contain(erc2)}
    for target in hierarchy.ercs:
        if target is erc1 or target is erc2:
            continue
        if target.label in sub1 or target.label in sub2:
            continue
        for gen in target.min_generators:
            gen_sp = set(species_list_to_names(gen))
            if (gen_sp.issubset(joint) and
                    not gen_sp.issubset(cl1) and
                    not gen_sp.issubset(cl2)):
                return True
    return False


def compute_complementarity_stats(ercs, hierarchy, RN):
    minprod, mincons = _compute_minprod_mincons(ercs, hierarchy, RN)
    comp_pairs = set(); pure_pairs = set()
    fund_pairs = set(); prod_ercs  = set(); cons_ercs = set()

    for e1, e2 in combinations(ercs, 2):
        if not _is_incomparable(e1, e2, hierarchy):
            continue
        s12 = _supply(e1, e2, RN)
        s21 = _supply(e2, e1, RN)
        if not s12 and not s21:
            continue
        pair_key = frozenset([e1.label, e2.label])
        comp_pairs.add(pair_key)
        if not _has_basic_synergy(e1, e2, hierarchy, RN):
            pure_pairs.add(pair_key)
        for s in s12:
            if e1 in minprod.get(s, []) and e2 in mincons.get(s, []):
                fund_pairs.add(pair_key)
                prod_ercs.add(e1.label); cons_ercs.add(e2.label)
        for s in s21:
            if e2 in minprod.get(s, []) and e1 in mincons.get(s, []):
                fund_pairs.add(pair_key)
                prod_ercs.add(e2.label); cons_ercs.add(e1.label)

    return {
        'n_complementary_pairs':      len(comp_pairs),
        'n_pure_complementary_pairs': len(pure_pairs),
        'n_fundamental_edges':        len(fund_pairs),
        'n_producer_ERCs':            len(prod_ercs),
        'n_consumer_ERCs':            len(cons_ercs),
    }


# =============================================================================
# File discovery
# =============================================================================

def collect_files(folders):
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


# =============================================================================
# Load existing CSVs for incremental resumption
# =============================================================================

def _load_existing(csv_path):
    if os.path.exists(csv_path):
        try:
            df = pd.read_csv(csv_path)
            return df, set(df['file'].astype(str))
        except Exception:
            pass
    return pd.DataFrame(), set()


syn_existing,  syn_done  = _load_existing(SYN_CSV)
comp_existing, comp_done = _load_existing(COMP_CSV)
already_done = syn_done & comp_done

all_files = collect_files(SCAN_FOLDERS)
print(f"Found {len(all_files)} .txt files across {len(SCAN_FOLDERS)} folders.")
print(f"synergy_stats.csv:        {len(syn_done)} cached")
print(f"complementarity_stats.csv: {len(comp_done)} cached")
print(f"Will skip {len(already_done)} networks present in both.\n")

# =============================================================================
# Main loop
# =============================================================================

syn_rows  = []
comp_rows = []
skipped   = []

for idx, (fpath, dataset_label) in enumerate(all_files):
    fname = os.path.basename(fpath)
    if fname in already_done:
        print(f"[{idx+1}/{len(all_files)}] {fname}  — CACHED, skip")
        continue
    print(f"\n[{idx+1}/{len(all_files)}] {fname}")

    try:
        RN   = read_txt(fpath)
        n_sp = len(RN.species())
        n_rx = len(RN.reactions())
        print(f"  {n_sp} species, {n_rx} reactions")

        if n_rx > MAX_REACTIONS:
            print(f"  SKIP: too many reactions ({n_rx} > {MAX_REACTIONS})")
            skipped.append((fname, f'too many reactions: {n_rx}'))
            continue

        t0 = time.time()
        ercs, _cached = load_ercs(fpath, RN, ERC)
        ercs = [e for e in ercs if len(e.get_closure_names(RN)) > 0]
        t_ercs = time.time() - t0

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

        hierarchy = ERC_Hierarchy(RN, ercs)

        n_pairs   = n_ercs * (n_ercs - 1) // 2
        n_triples = n_ercs * (n_ercs - 1) * (n_ercs - 2) // 6

        # ── Synergy ───────────────────────────────────────────────────────────
        t_s0 = time.time()
        syn  = compute_synergy_stats(ercs, hierarchy, RN)
        t_syn = time.time() - t_s0
        print(f"  Syn pairs   Basic={syn['n_basic_pairs']}  Max={syn['n_maximal_pairs']}"
              f"  Fund={syn['n_fundamental_pairs']}  ({t_syn:.1f}s)")
        if syn['syn3_computed']:
            print(f"  3-Syn triples  Basic={syn['n_basic_3syn_triples']}"
                  f"  Fund={syn['n_fundamental_3syn_triples']}")

        syn_rows.append({
            'file': fname, 'dataset': dataset_label,
            'n_species': n_sp, 'n_reactions': n_rx, 'n_ercs': n_ercs,
            'n_pairs_max': n_pairs, 'n_triples_max': n_triples,
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
            'time_ercs':      round(t_ercs, 2),
            'time_synergies': round(t_syn,  2),
        })

        # ── Complementarity ───────────────────────────────────────────────────
        t_c0 = time.time()
        cst  = compute_complementarity_stats(ercs, hierarchy, RN)
        t_comp = time.time() - t_c0
        print(f"  Comp pairs  All={cst['n_complementary_pairs']}"
              f"  Pure={cst['n_pure_complementary_pairs']}"
              f"  Fund={cst['n_fundamental_edges']}  ({t_comp:.1f}s)")

        comp_rows.append({
            'file': fname, 'dataset': dataset_label,
            'n_species': n_sp, 'n_reactions': n_rx, 'n_ercs': n_ercs,
            'n_pairs_max': n_pairs,
            'n_complementary_pairs':      cst['n_complementary_pairs'],
            'n_pure_complementary_pairs': cst['n_pure_complementary_pairs'],
            'n_fundamental_edges':        cst['n_fundamental_edges'],
            'n_producer_ERCs':            cst['n_producer_ERCs'],
            'n_consumer_ERCs':            cst['n_consumer_ERCs'],
            'ratio_complementary': cst['n_complementary_pairs']      / n_pairs if n_pairs > 0 else 0,
            'ratio_pure':          cst['n_pure_complementary_pairs'] / n_pairs if n_pairs > 0 else 0,
            'ratio_fundamental':   cst['n_fundamental_edges']        / n_pairs if n_pairs > 0 else 0,
            'ratio_producers':     cst['n_producer_ERCs'] / n_ercs if n_ercs > 0 else 0,
            'ratio_consumers':     cst['n_consumer_ERCs'] / n_ercs if n_ercs > 0 else 0,
            'time_ercs':            round(t_ercs,  2),
            'time_complementarity': round(t_comp,  2),
        })

        # Flush both CSVs after every network
        _sf = pd.concat([syn_existing,  pd.DataFrame(syn_rows)],  ignore_index=True) \
              if not syn_existing.empty  else pd.DataFrame(syn_rows)
        _cf = pd.concat([comp_existing, pd.DataFrame(comp_rows)], ignore_index=True) \
              if not comp_existing.empty else pd.DataFrame(comp_rows)
        _sf.to_csv(SYN_CSV,  index=False)
        _cf.to_csv(COMP_CSV, index=False)
        print(f"  → CSVs updated  (syn={len(_sf)}, comp={len(_cf)} total records)")

    except Exception as exc:
        import traceback
        print(f"  ERROR: {exc}")
        traceback.print_exc()
        skipped.append((fname, str(exc)))

# =============================================================================
# Final save
# =============================================================================

new_syn  = pd.DataFrame(syn_rows)
new_comp = pd.DataFrame(comp_rows)

df_syn  = pd.concat([syn_existing,  new_syn],  ignore_index=True) \
          if not syn_existing.empty  and not new_syn.empty  \
          else (new_syn  if not new_syn.empty  else syn_existing)
df_comp = pd.concat([comp_existing, new_comp], ignore_index=True) \
          if not comp_existing.empty and not new_comp.empty \
          else (new_comp if not new_comp.empty else comp_existing)

df_syn.to_csv(SYN_CSV,  index=False)
df_comp.to_csv(COMP_CSV, index=False)

print(f"\n=== Done ===")
print(f"synergy_stats.csv:        {len(df_syn)} records  ({len(new_syn)} new)")
print(f"complementarity_stats.csv: {len(df_comp)} records  ({len(new_comp)} new)")
if skipped:
    print(f"Skipped {len(skipped)} networks:")
    for fn, reason in skipped:
        print(f"  {fn}: {reason}")
