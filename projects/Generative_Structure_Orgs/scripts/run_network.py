#!/usr/bin/env python3
"""
run_network.py
===============
Compute the ERC hierarchy of a single reaction network and overlay
synergy and/or complementarity relations, restricted to whichever
classification levels are enabled in the CONFIGURATION block below.

Each relation type (synergy, complementarity) has its own three-tier
classification in this codebase; both are exposed here under the SAME
three toggles for a uniform mental model:

  BASIC       — the loosest, broadest relation.
                Synergy:         Def. 18 "basic synergy" (joint closure
                                  covers a minimal generator of the target
                                  that neither reactant covers alone).
                Complementarity: any directional supply prod->req between
                                  two incomparable ERCs, INCLUDING pairs
                                  that also happen to synergize (the
                                  'complementary' tier in
                                  script_erc_complementarity_viz.py).

  MAXMIN      — the refined, minimal/maximal-generator-based tier.
                Synergy:         Def. 19 "maximal synergy" (no other basic
                                  synergy for the same reactant pair has a
                                  strictly bigger target).
                Complementarity: "purely complementary" -- a basic supply
                                  between a pair that does NOT also
                                  synergize (the 'pure' tier).

  FUNDAMENTAL — the strictest tier.
                Synergy:         Def. 20 -- no strictly smaller reactant
                                  pair already achieves a maximal synergy
                                  to the same target.
                Complementarity: the supply species has the producer in
                                  minprod(s) AND the consumer in
                                  mincons(s) -- the two ERCs are each
                                  individually minimal for that species.

By default only the two FUNDAMENTAL toggles are True; every other toggle
is False, matching each relation's own "no shortcut / no redundant path"
notion of minimality.

Rendering (identical conventions to script_erc_synergy_viz.py /
script_erc_complementarity_viz.py, combined on one shared ERC hierarchy
layout):
  - ERC nodes: size ~ closure size; colour = maintenance class
    (green = self-maintaining, orange = semi-self-maintaining only,
    steel = neither / not reactive).
  - Synergy: diamond junction, converging arrows from both reactants,
    diverging arrow to the target. Colour by level (blue/orange/green
    for basic/maxmin/fundamental).
  - Complementarity: square junction labelled with the supplied species,
    arrow producer -> junction -> consumer. Colour by level
    (blue/orange/green for basic/maxmin/fundamental).

Both relation types are computed only if at least one of their three
toggles is True (skips the O(|ERCs|^2)-plus computation entirely
otherwise).

Correctness note
-----------------
The synergy/complementarity computations below are copied verbatim from
script_erc_synergy_viz.py / script_erc_complementarity_viz.py, which
already carry corrected implementations (the library's
get_maximal_synergies() has an inverted containment check).

HOW TO RUN
----------
Edit the CONFIGURATION block, then:
  python projects/Generative_Structure_Orgs/scripts/run_network.py
"""

import os
import sys
import time
from itertools import combinations
from collections import defaultdict

import numpy as np
from scipy.optimize import linprog
from matplotlib.patches import FancyArrowPatch
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import networkx as nx

# -- Path setup ----------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, '..', '..', '..'))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))
sys.path.insert(0, _SCRIPT_DIR)

from pyCOT.io.functions import read_txt
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from pyCOT.analysis.SORN_Generators import is_semi_self_maintaining
from utils_ercs import load_ercs

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  CONFIGURATION — edit this block, then run                               ║
# ╚══════════════════════════════════════════════════════════════════════════╝

RN_FILE = os.path.join(_PYCOT_ROOT, 'data', 'Examples_tests', 'FarmVariants',
                        'Farm_agro_stages.txt')
# Alternatives:
# RN_FILE = os.path.join(_PYCOT_ROOT, 'data', 'biomodels', 'BiGG', 'bigg_iAF692.txt')

# Synergy tiers to compute and render.
SHOW_SYN_BASIC       = False
SHOW_SYN_MAXMIN      = False   # "maximal" synergy
SHOW_SYN_FUNDAMENTAL = True

# Complementarity tiers to compute and render.
SHOW_COMP_BASIC       = False   # "complementary" (co-occurs with synergy)
SHOW_COMP_MAXMIN      = False   # "pure" (not synergetic)
SHOW_COMP_FUNDAMENTAL = True

COMPUTE_SYN  = SHOW_SYN_BASIC or SHOW_SYN_MAXMIN or SHOW_SYN_FUNDAMENTAL
COMPUTE_COMP = SHOW_COMP_BASIC or SHOW_COMP_MAXMIN or SHOW_COMP_FUNDAMENTAL

# -- Visual parameters ---------------------------------------------------------
NODE_SIZE_BASE  = 400
NODE_SIZE_SCALE = 120
COL_SM   = '#27AE60'
COL_SSM  = '#E67E22'
COL_DEF  = '#AED6F1'
SM_EPS   = 1e-6

COL_SYN = {'basic': '#1C5FB7', 'maximal': '#D76E0B', 'fundamental': "#10ABBC"}
COL_COMP = {'complementary': '#3498DB', 'pure': '#E67E22', 'fundamental': "#AE27AE"}

# Level -> whether it's enabled, keyed by each script's own internal name.
SYN_ENABLED = {
    'basic': SHOW_SYN_BASIC,
    'maximal': SHOW_SYN_MAXMIN,
    'fundamental': SHOW_SYN_FUNDAMENTAL,
}
COMP_ENABLED = {
    'complementary': SHOW_COMP_BASIC,
    'pure': SHOW_COMP_MAXMIN,
    'fundamental': SHOW_COMP_FUNDAMENTAL,
}

JUNC_SIZE_SYN  = 80
JUNC_SIZE_COMP = 90
ARROW_LW_SYN    = {'basic': 0.7, 'maximal': 0.7, 'fundamental': 0.7}
ARROW_ALPHA_SYN = {'basic': 0.7, 'maximal': 0.7, 'fundamental': 0.7}
ARROW_LW_COMP    = {'complementary': 0.7, 'pure': 0.9, 'fundamental': 1.1}
ARROW_ALPHA_COMP = {'complementary': 0.65, 'pure': 0.75, 'fundamental': 0.90}

SHRINK_ERC  = 12
SHRINK_JUNC = 5
PUSH_OFF_LEVEL_SYN  = 0.40
PUSH_OFF_LEVEL_COMP = 0.8


# ===========================================================================
# Node classification helpers (identical across the source scripts)
# ===========================================================================

def _is_reactive(RN, species_set):
    if not species_set:
        return False
    return any(r.support_indices() for r in RN.get_reactions_from_species(species_set))


def _is_sm(RN, species_set):
    sub = RN.sub_reaction_network(species_set)
    S = np.asarray(sub.stoichiometry_matrix(), dtype=float)
    _, n_rx = S.shape
    if n_rx == 0:
        return False
    res = linprog(np.zeros(n_rx),
                   A_ub=-S, b_ub=SM_EPS * S.sum(axis=1),
                   bounds=[(0, None)] * n_rx, method='highs')
    return res.status == 0


# ===========================================================================
# Corrected synergy detection (copied verbatim from script_erc_synergy_viz.py)
# ===========================================================================

def _compute_basic(erc1, erc2, hierarchy, RN):
    """Basic synergies (Def. 18)."""
    if (erc1 in hierarchy.get_contain(erc2) or
            erc2 in hierarchy.get_contain(erc1)):
        return []
    cl1 = erc1.get_closure_names(RN)
    cl2 = erc2.get_closure_names(RN)
    joint = cl1 | cl2
    sub1 = {e.label for e in hierarchy.get_contain(erc1)}
    sub2 = {e.label for e in hierarchy.get_contain(erc2)}
    result = []
    for target in hierarchy.ercs:
        if target is erc1 or target is erc2:
            continue
        if target.label in sub1 or target.label in sub2:
            continue
        for gen in target.min_generators:
            gen_sp = set(species_list_to_names(gen))
            if gen_sp.issubset(joint) and not gen_sp.issubset(cl1) and not gen_sp.issubset(cl2):
                result.append((erc1, erc2, target))
                break
    return result


def _filter_maximal_v2(basics, RN):
    """Definition 19: maximal synergy (biggest target per reactant pair)."""
    if not basics:
        return []
    by_pair = defaultdict(list)
    for e1, e2, t in basics:
        by_pair[tuple(sorted([e1.label, e2.label]))].append((e1, e2, t))
    result = []
    for syns in by_pair.values():
        cls = {s[2].label: s[2].get_closure_names(RN) for s in syns}
        for e1, e2, T in syns:
            if not any(cls[T.label] < cls[T2.label]
                       for _, _, T2 in syns if T2.label != T.label):
                result.append((e1, e2, T))
    return result


def _filter_fundamental(maximals, hierarchy, RN):
    """Definition 20: fundamental synergy (no strictly smaller pair suffices)."""
    if not maximals:
        return []

    def _desc(erc):
        if erc.label not in hierarchy.graph:
            return {erc.label}
        return {erc.label} | set(nx.descendants(hierarchy.graph, erc.label))

    desc = {}
    l2erc = {e.label: e for e in hierarchy.ercs}
    for e1, e2, _ in maximals:
        for e in (e1, e2):
            if e.label not in desc:
                desc[e.label] = _desc(e)

    result = []
    for e1, e2, target in maximals:
        orig = tuple(sorted([e1.label, e2.label]))
        is_fund = True
        for l1 in desc[e1.label]:
            if not is_fund:
                break
            for l2 in desc[e2.label]:
                if l1 == l2 or tuple(sorted([l1, l2])) == orig:
                    continue
                s1, s2 = l2erc.get(l1), l2erc.get(l2)
                if s1 is None or s2 is None:
                    continue
                sub_b = _compute_basic(s1, s2, hierarchy, RN)
                sub_m = _filter_maximal_v2(sub_b, RN)
                if any(s[2].label == target.label for s in sub_m):
                    is_fund = False
                    break
        if is_fund:
            result.append((e1, e2, target))
    return result


def compute_all_synergies(ercs, hierarchy, RN):
    """Returns (all_basics, all_maximals, all_fundamentals) as lists of (e1,e2,target)."""
    all_b, all_m, all_f = [], [], []
    pairs = list(combinations(ercs, 2))
    print(f"  Checking {len(pairs)} ERC pairs for synergy...")
    for e1, e2 in pairs:
        b = _compute_basic(e1, e2, hierarchy, RN)
        if not b:
            continue
        all_b.extend(b)
        m = _filter_maximal_v2(b, RN)
        all_m.extend(m)
        f = _filter_fundamental(m, hierarchy, RN)
        all_f.extend(f)
    return all_b, all_m, all_f


# ===========================================================================
# Complementarity detection (copied verbatim from
# script_erc_complementarity_viz.py)
# ===========================================================================

def _has_basic_synergy(erc1, erc2, hierarchy, RN):
    """Return True iff (erc1,erc2) has at least one basic synergy."""
    if (erc1 in hierarchy.get_contain(erc2) or
            erc2 in hierarchy.get_contain(erc1)):
        return False
    cl1 = erc1.get_closure_names(RN)
    cl2 = erc2.get_closure_names(RN)
    joint = cl1 | cl2
    sub1 = {e.label for e in hierarchy.get_contain(erc1)}
    sub2 = {e.label for e in hierarchy.get_contain(erc2)}
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


def _is_incomparable(erc1, erc2, hierarchy):
    return (erc1 not in hierarchy.get_contain(erc2) and
            erc2 not in hierarchy.get_contain(erc1))


def _supply(erc_prod, erc_cons, RN):
    """supl(erc_prod, erc_cons) = prod(R_{erc_prod}) ∩ req(erc_cons)."""
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

    def _minimal(label_set):
        minimal = []
        for lbl in label_set:
            erc = label_to_erc[lbl]
            descendants = {e.label for e in hierarchy.get_contain(erc)}
            if not descendants.intersection(label_set):
                minimal.append(erc)
        return minimal

    minprod = {s: _minimal(lbls) for s, lbls in producers.items()}
    mincons = {s: _minimal(lbls) for s, lbls in consumers.items()}
    return minprod, mincons


def _classify_supply(e_prod, e_cons, supply_species, pair_is_syn, minprod, mincons):
    """Highest classification level for a directional supply e_prod -> e_cons."""
    for s in supply_species:
        if e_prod in minprod.get(s, []) and e_cons in mincons.get(s, []):
            return 'fundamental'
    if not pair_is_syn:
        return 'pure'
    return 'complementary'


def compute_complementarity_entries(ercs, hierarchy, RN):
    """
    Return a list of dicts, one per directional supply (producer -> consumer):
      {'producer': ERC, 'consumer': ERC, 'species': set[str], 'level': str}
    """
    minprod, mincons = _compute_minprod_mincons(ercs, hierarchy, RN)
    entries = []
    print(f"  Checking {len(list(combinations(ercs, 2)))} ERC pairs for complementarity...")
    for e1, e2 in combinations(ercs, 2):
        if not _is_incomparable(e1, e2, hierarchy):
            continue
        s12 = _supply(e1, e2, RN)
        s21 = _supply(e2, e1, RN)
        if not s12 and not s21:
            continue
        pair_syn = _has_basic_synergy(e1, e2, hierarchy, RN)
        if s12:
            level = _classify_supply(e1, e2, s12, pair_syn, minprod, mincons)
            entries.append({'producer': e1, 'consumer': e2, 'species': s12, 'level': level})
        if s21:
            level = _classify_supply(e2, e1, s21, pair_syn, minprod, mincons)
            entries.append({'producer': e2, 'consumer': e1, 'species': s21, 'level': level})
    return entries


# ===========================================================================
# Load network, compute ERCs, build hierarchy
# ===========================================================================
print(f"Loading: {RN_FILE}")
RN = read_txt(RN_FILE)
print(f"  {len(RN.species())} species, {len(RN.reactions())} reactions")

print("Computing ERCs...")
t0 = time.time()
ercs, from_cache = load_ercs(RN_FILE, RN, ERC)
print(f"  {len(ercs)} ERCs  ({time.time() - t0:.1f}s, {'cached' if from_cache else 'fresh'})")

print("Building hierarchy...")
hierarchy = ERC_Hierarchy(RN, ercs)
G = hierarchy.graph

# -- Classify nodes -------------------------------------------------------------
print("Classifying nodes (SSM -> SM)...")
node_sizes, node_colors = {}, {}
for erc in ercs:
    cl = erc.get_closure(RN)
    n_sp = len(cl)
    node_sizes[erc.label] = NODE_SIZE_BASE + NODE_SIZE_SCALE * n_sp
    if not _is_reactive(RN, cl):
        node_colors[erc.label] = COL_DEF
    elif not is_semi_self_maintaining(RN, cl):
        node_colors[erc.label] = COL_DEF
    elif _is_sm(RN, cl):
        node_colors[erc.label] = COL_SM
    else:
        node_colors[erc.label] = COL_SSM

# -- Compute synergies (only the selected tiers get rendered, but all three
#    tiers are always computed together since maximal/fundamental are
#    filters over basic -- there is no cheaper way to get just one tier) ------
all_b = all_m = all_f = []
if COMPUTE_SYN:
    print("Computing synergies...")
    t0 = time.time()
    all_b, all_m, all_f = compute_all_synergies(ercs, hierarchy, RN)
    print(f"  Basic={len(all_b)}  Maximal={len(all_m)}  Fundamental={len(all_f)}"
          f"  ({time.time() - t0:.1f}s)")

syn_level = {}
syn_reactants = {}
for level_name, triples in (('basic', all_b), ('maximal', all_m), ('fundamental', all_f)):
    if not SYN_ENABLED[level_name]:
        continue
    for e1, e2, t in triples:
        k = (tuple(sorted([e1.label, e2.label])), t.label)
        syn_level[k] = level_name          # later tiers overwrite with their own level
        syn_reactants[k] = (e1.label, e2.label)

# -- Compute complementarities (same reasoning: fundamental/pure are
#    sub-classifications of the same entries, computed together) -------------
comp_entries = []
if COMPUTE_COMP:
    print("Computing complementarities...")
    t0 = time.time()
    comp_entries = compute_complementarity_entries(ercs, hierarchy, RN)
    t_comp = time.time() - t0
    n_fund = sum(1 for e in comp_entries if e['level'] == 'fundamental')
    n_pure = sum(1 for e in comp_entries if e['level'] == 'pure')
    n_comp = sum(1 for e in comp_entries if e['level'] == 'complementary')
    print(f"  Fundamental={n_fund}  Pure={n_pure}  Complementary(+synergy)={n_comp}"
          f"  ({t_comp:.1f}s)")

comp_entries_shown = [e for e in comp_entries if COMP_ENABLED[e['level']]]

# ===========================================================================
# Layout: level-based (identical to script_erc_hierarchy.py and both viz
# scripts, so this reads the same as the single-relation plots)
# ===========================================================================
levels = ERC.get_node_levels(G)
level_nodes = defaultdict(list)
for node, lvl in levels.items():
    level_nodes[lvl].append(node)

pos = {}
for lvl, nodes in level_nodes.items():
    nodes.sort(key=lambda n: len(nx.ancestors(G, n)), reverse=True)
    for i, node in enumerate(nodes):
        pos[node] = np.array([(i - (len(nodes) - 1) / 2) * 2.0, lvl * 2.0])

level_ys = sorted(set(float(p[1]) for p in pos.values()))


def _clear_of_levels(y, y_target, step, max_iter=8):
    direction = 1 if y_target >= y else -1
    for _ in range(max_iter):
        if not any(abs(y - ly) < 0.22 for ly in level_ys):
            break
        y += direction * step
    return y


# -- Synergy junction positions (3-arrow diamond topology) --------------------
syn_junc_pos = {}
syn_pair_junc_count = defaultdict(int)
for key, (l1, l2) in syn_reactants.items():
    _, t_label = key
    if l1 not in pos or l2 not in pos or t_label not in pos:
        continue
    p1, p2, pt = pos[l1], pos[l2], pos[t_label]
    x_junc = (p1[0] + p2[0]) / 2.0
    y_hi = max(p1[1], p2[1])
    y_T = float(pt[1])
    if abs(y_T - y_hi) > 0.5:
        y_junc = (y_hi + y_T) / 2.0
    else:
        y_sign = 1 if y_T >= y_hi else -1
        y_junc = y_hi + y_sign * 0.7
    y_junc = _clear_of_levels(y_junc, y_T, PUSH_OFF_LEVEL_SYN)
    base = np.array([x_junc, y_junc])
    pair_key = tuple(sorted([l1, l2]))
    count = syn_pair_junc_count[pair_key]
    if count > 0:
        sign = 1 if count % 2 == 1 else -1
        offset = ((count + 1) // 2) * 0.42
        base = base + np.array([sign * offset, 0.0])
    syn_pair_junc_count[pair_key] += 1
    syn_junc_pos[key] = base

# -- Complementarity junction positions (2-arrow square topology) -------------
comp_junc_pos = {}
comp_pair_count = defaultdict(int)
for idx, entry in enumerate(comp_entries_shown):
    lp, lc = entry['producer'].label, entry['consumer'].label
    if lp not in pos or lc not in pos:
        continue
    pp, pc = pos[lp], pos[lc]
    mid = (pp + pc) * 0.5
    direction = pc - pp
    dist = np.linalg.norm(direction)
    perp = np.array([-direction[1], direction[0]]) / dist if dist > 1e-9 else np.array([0.0, 1.0])
    mid_y = float(mid[1])
    if any(abs(mid_y - ly) < 0.25 for ly in level_ys):
        sign = 1 if hash(frozenset([lp, lc])) % 2 == 0 else -1
        mid = mid + sign * PUSH_OFF_LEVEL_COMP * perp
    pair_key = frozenset([lp, lc])
    count = comp_pair_count[pair_key]
    if count > 0:
        sign = 1 if count % 2 == 1 else -1
        mid = mid + sign * (count + 1) * 0.32 * perp
    comp_pair_count[pair_key] += 1
    comp_junc_pos[idx] = mid


# ===========================================================================
# Plot
# ===========================================================================
ordered_nodes = list(G.nodes())
sizes_list = [node_sizes.get(n, NODE_SIZE_BASE) for n in ordered_nodes]
colors_list = [node_colors.get(n, COL_DEF) for n in ordered_nodes]

fig, ax = plt.subplots(figsize=(14, 9))

scatter = nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=ordered_nodes,
                                  node_color=colors_list, node_size=sizes_list, alpha=0.92)
scatter.set_zorder(3)
nx.draw_networkx_edges(G, pos, ax=ax, edge_color='#555555',
                        arrows=True, arrowsize=14, width=1.1, alpha=0.65)
nx.draw_networkx_labels(G, pos, ax=ax, font_size=12, font_weight='bold')


def _farrow(src, dst, col, lw, alpha, shrinkA=SHRINK_ERC, shrinkB=SHRINK_ERC):
    patch = FancyArrowPatch(posA=tuple(src), posB=tuple(dst), arrowstyle='->',
                             lw=lw, color=col, alpha=alpha,
                             shrinkA=shrinkA, shrinkB=shrinkB,
                             mutation_scale=10, zorder=1)
    ax.add_patch(patch)


# -- Draw synergies (diamond junctions) ----------------------------------------
for key, stype in syn_level.items():
    if key not in syn_junc_pos:
        continue
    l1, l2 = syn_reactants[key]
    _, t_label = key
    if l1 not in pos or l2 not in pos or t_label not in pos:
        continue
    col, lw, alpha = COL_SYN[stype], ARROW_LW_SYN[stype], ARROW_ALPHA_SYN[stype]
    junc = syn_junc_pos[key]
    ax.scatter(*junc, s=JUNC_SIZE_SYN, c=col, marker='D',
               edgecolors='white', linewidths=0.8, alpha=0.95, zorder=3)
    _farrow(pos[l1], junc, col, lw, alpha, shrinkA=SHRINK_ERC, shrinkB=SHRINK_JUNC)
    _farrow(pos[l2], junc, col, lw, alpha, shrinkA=SHRINK_ERC, shrinkB=SHRINK_JUNC)
    _farrow(junc, pos[t_label], col, lw, alpha, shrinkA=SHRINK_JUNC, shrinkB=SHRINK_ERC)

# -- Draw complementarities (square junctions) ---------------------------------
for idx, entry in enumerate(comp_entries_shown):
    if idx not in comp_junc_pos:
        continue
    lp, lc, level = entry['producer'].label, entry['consumer'].label, entry['level']
    if lp not in pos or lc not in pos:
        continue
    col, lw, alpha = COL_COMP[level], ARROW_LW_COMP[level], ARROW_ALPHA_COMP[level]
    jpos = comp_junc_pos[idx]
    ax.scatter(*jpos, s=JUNC_SIZE_COMP, c=col, marker='s',
               edgecolors='white', linewidths=0.8, alpha=0.95, zorder=4)
    species_label = ', '.join(sorted(entry['species']))
    if len(species_label) > 20:
        species_label = species_label[:18] + '…'
    ax.text(jpos[0], jpos[1] + 0.18, species_label, ha='center', va='bottom',
            fontsize=6, color=col, zorder=5, alpha=min(alpha + 0.1, 1.0))
    _farrow(pos[lp], jpos, col, lw, alpha, shrinkA=SHRINK_ERC, shrinkB=SHRINK_JUNC)
    _farrow(jpos, pos[lc], col, lw, alpha, shrinkA=SHRINK_JUNC, shrinkB=SHRINK_ERC)

# -- Hover tooltip (ERC nodes) --------------------------------------------------
annot = ax.annotate("", xy=(0, 0), xytext=(15, 15), textcoords="offset points",
                     bbox=dict(boxstyle="round", fc="white", ec="0.5", alpha=0.9),
                     visible=False, zorder=10)
erc_by_label = {erc.label: erc for erc in ercs}


def _on_hover(event):
    if event.inaxes != ax:
        return
    cont, ind = scatter.contains(event)
    if cont:
        node = ordered_nodes[ind["ind"][0]]
        erc = erc_by_label.get(node)
        cl = species_list_to_names(erc.get_closure(RN)) if erc else []
        req = sorted(erc.get_required_species(RN)) if erc else []
        prod = sorted(erc.get_produced_species(RN)) if erc else []
        annot.xy = tuple(pos[node])
        annot.set_text(f"{node}\nclosure: {cl}\nreq: {req}\nprod: {prod}")
        annot.set_visible(True)
    else:
        annot.set_visible(False)
    fig.canvas.draw_idle()


fig.canvas.mpl_connect("motion_notify_event", _on_hover)

# -- Legends ---------------------------------------------------------------
maintenance_legend = ax.legend(
    handles=[
        mpatches.Patch(facecolor=COL_SM, label='Self-maintaining (SM)'),
        mpatches.Patch(facecolor=COL_SSM, label='Semi-self-maintaining (SSM)'),
        mpatches.Patch(facecolor=COL_DEF, label='Neither / not reactive'),
    ],
    loc='upper right', fontsize=10, framealpha=0.9,
    title='Maintenance class', title_fontsize=10,
)
ax.add_artist(maintenance_legend)

if COMPUTE_SYN:
    syn_handles = []
    for level_name, count in (('fundamental', len(all_f)), ('maximal', len(all_m)), ('basic', len(all_b))):
        if SYN_ENABLED[level_name]:
            syn_handles.append(Line2D([0], [0], color=COL_SYN[level_name], lw=2.0,
                                       label=f'{level_name.capitalize()} synergy ({count})'))
    if syn_handles:
        syn_legend = ax.legend(handles=syn_handles, loc='lower right', fontsize=10,
                                framealpha=0.9, title='Synergy type (◆)', title_fontsize=10)
        ax.add_artist(syn_legend)

if COMPUTE_COMP:
    comp_handles = []
    for level_name, label in (('fundamental', 'Fundamental'), ('pure', 'Pure (not synergetic)'),
                               ('complementary', 'Complementary (+synergy)')):
        if COMP_ENABLED[level_name]:
            n = sum(1 for e in comp_entries_shown if e['level'] == level_name)
            comp_handles.append(Line2D([0], [0], color=COL_COMP[level_name], lw=2.0,
                                        label=f'{label} ({n})'))
    if comp_handles:
        comp_legend = ax.legend(handles=comp_handles, loc='upper left', fontsize=10,
                                 framealpha=0.9, title='Complementarity type (■)',
                                 title_fontsize=10)
        ax.add_artist(comp_legend)

all_counts = sorted({len(erc.get_closure(RN)) for erc in ercs})
if len(all_counts) > 4:
    tick_idx = [0, len(all_counts) // 3, 2 * len(all_counts) // 3, -1]
    size_ticks = sorted({all_counts[i] for i in tick_idx})
else:
    size_ticks = all_counts
size_handles = [
    ax.scatter([], [], s=NODE_SIZE_BASE + NODE_SIZE_SCALE * n, color='#888888',
               alpha=0.85, label=f'{n} species')
    for n in size_ticks
]
ax.legend(handles=size_handles, loc='lower left', fontsize=10, framealpha=0.9,
          title='Closure size', title_fontsize=10, labelspacing=1.2, handletextpad=1.0)

# -- Title -----------------------------------------------------------------
net_name = os.path.splitext(os.path.basename(RN_FILE))[0]
syn_tiers = [k for k, v in SYN_ENABLED.items() if v]
comp_tiers = [k for k, v in COMP_ENABLED.items() if v]
ax.set_title(
    f"ERC Hierarchy — {net_name}\n"
    f"Synergy tiers: {', '.join(syn_tiers) or 'none'}   |   "
    f"Complementarity tiers: {', '.join(comp_tiers) or 'none'}",
    fontsize=12)
ax.axis('off')
plt.tight_layout()

_out_dir = os.path.join(_SCRIPT_DIR, '..', 'outputs', 'run_network')
os.makedirs(_out_dir, exist_ok=True)
_out_path = os.path.join(_out_dir, f'run_network_{net_name}.png')
plt.savefig(_out_path, dpi=150, bbox_inches='tight')
print(f"Saved: {os.path.abspath(_out_path)}")
plt.show()
