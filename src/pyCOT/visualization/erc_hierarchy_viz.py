"""
erc_hierarchy_viz.py — ERC hierarchy + fundamental relations visualization.

Renders the generative structure of a reaction network — the object the
Persistent_Modules_Generator.py module docstring points to as worth
studying on its own, before any organization/EPM/ESPM computation — as an
interactive graph:

  - gray edges   : ERC containment (the Hasse diagram of the hierarchy)
  - orange nodes : fundamental synergies (one diamond per triad: lines from
                    the two contributing ERCs, arrow to the target ERC
                    they jointly unlock)
  - blue dashed  : fundamental complementarities (producer ERC -> consumer
                    ERC, labeled with the species being supplied)

Ported from projects/COT_Fundamental_Generators_Exploration/cot_gen/deep_report.py and
deep_report_viz.py (compute_hierarchy_stats + plot_hierarchy_overview),
trimmed to just the hierarchy/relations visualization (the degeneracy-
tracking and SO-lattice plotting in the original stay in that project,
which is a research/instrumentation tool rather than a library-visualization
concern).

Public API
----------
compute_hierarchy_stats(ercs, hier, syn_result, comp_result) -> HierarchyStats
plot_erc_hierarchy(ercs, hier, syn_result, comp_result, out_path, *,
                    stats=None, highlight_ercs=None, highlight_label="highlighted",
                    max_nodes=250, show_synergy=True, show_complementarity=True,
                    title="ERC hierarchy") -> str (path to the written HTML file)

Quick start
-----------
    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations import (
        build_rndata, compute_ercs, build_hierarchy,
        compute_synergies_basis_first, compute_complementarities,
    )
    from pyCOT.visualization.erc_hierarchy_viz import plot_erc_hierarchy

    rn = read_txt("data/Examples_tests/Centler2006_EcoliSugar/centler_glucose.txt")
    rn_data = build_rndata(rn, network_id="centler_glucose")
    ercs = compute_ercs(rn_data)
    hier = build_hierarchy(ercs)
    syn = compute_synergies_basis_first(ercs, hier)
    comp = compute_complementarities(ercs, hier, syn)
    plot_erc_hierarchy(ercs, hier, syn, comp, "erc_hierarchy.html")
"""
from __future__ import annotations

import os
from dataclasses import dataclass

# Shared palette -- consistent with cot_gen/explorer.py's relation colors,
# so plots from either tool read the same way.
COLOR_CONTAINMENT = "#7f8c8d"
COLOR_SYNERGY = "#e67e22"
COLOR_COMPLEMENTARITY = "#2980b9"
COLOR_HIGHLIGHT = "#27ae60"
COLOR_NEUTRAL = "#d0d3d4"


@dataclass
class HierarchyStats:
    """Per-ERC structural statistics, used both for the plot and for reports."""
    n_ercs: int
    n_persistent: int
    sizes: list          # species count per ERC
    levels: list          # containment level per ERC (0 = subset-minimal)
    n_levels: int
    syn_degree: list      # fundamental synergies touching each ERC
    comp_out_degree: list  # as fundamental producer
    comp_in_degree: list   # as fundamental consumer
    n_fundamental_syn: int
    n_fundamental_comp: int
    hasse_edges: int
    n_comparable_pairs: int
    n_incomparable_pairs: int


def e0_cover_targets(hier, n: int, e0_index: int) -> list:
    """Subset-minimal non-E0 ERCs: the ERCs that directly cover E0."""
    return [i for i in range(n) if i != e0_index and not hier.children[i]]


def relation_level_counts(stats: HierarchyStats, syn_result, comp_result) -> dict:
    """
    Count fundamental relations by the containment levels of the ERCs involved.

    Synergy (E_i, E_j) -> E_k is a triad: it counts under
    (a, b, c) = (sorted levels of the two contributing ERCs, level of the
    target). Complementarity producer -> consumer counts under the unordered
    pair (a, b) of the two ERCs' levels, so (1, 2) and (2, 1) share a bin.

    Returns {"syn": {(a, b, c): int}, "comp": {(a, b): int}}, keys sorted.
    """
    lv = stats.levels
    syn: dict = {}
    for (i, j, k) in synergy_triads(syn_result):   # same triads the graph draws
        a, b = sorted((lv[i], lv[j]))
        syn[(a, b, lv[k])] = syn.get((a, b, lv[k]), 0) + 1
    comp: dict = {}
    for fc in comp_result.fundamental:
        key = tuple(sorted((lv[fc.prod_idx], lv[fc.cons_idx])))
        comp[key] = comp.get(key, 0) + 1
    return {"syn": dict(sorted(syn.items())), "comp": dict(sorted(comp.items()))}


def synergy_triads(syn_result, node_set=None) -> list:
    """Distinct fundamental synergies as (i, j, k) with i < j, all in node_set."""
    triads = {(min(st.i, st.j), max(st.i, st.j), st.k) for st in syn_result.fundamental}
    if node_set is not None:
        triads = {t for t in triads if all(x in node_set for x in t)}
    return sorted(triads)


def syn_triad_label(key, arrow: str = "->") -> str:
    """(a, b, c) -> 'a+b->c' (levels of the two contributors, then the target).
    ASCII arrow by default so it prints on any console; HTML passes '→'."""
    return f"{key[0]}+{key[1]}{arrow}{key[2]}"


def _compute_containment_levels(hier, n: int) -> list:
    """Level 0 = subset-minimal ERCs; level(i) = 1 + max(level(child))."""
    order = sorted(range(n), key=lambda i: len(hier.descendants[i]))
    level = [0] * n
    children_of = [[] for _ in range(n)]
    for i in range(n):
        for p in hier.parents[i]:
            children_of[p].append(i)
    for i in order:
        if children_of[i]:
            level[i] = 1 + max(level[c] for c in children_of[i])
    return level


def compute_hierarchy_stats(ercs, hier, syn_result, comp_result, *, e0_index=None) -> HierarchyStats:
    """
    Compute per-ERC level/degree statistics used by plot_erc_hierarchy.

    e0_index : index of the inflow ERC E0, if the network has inflow.
        build_hierarchy works on E0-stripped (quotiented) species masks, so
        E0 comes out isolated, with no parents. In the full network,
        E0 is contained in every ERC's closure, i.e. it is the bottom of the
        hierarchy. When e0_index is given, E0 is put at level 0, every other
        ERC is shifted up one level, and the implicit E0 -> (subset-minimal
        ERC) cover edges are counted in hasse_edges.
    """
    n = len(ercs)
    sizes = [e.size() for e in ercs]
    levels = _compute_containment_levels(hier, n)
    if e0_index is not None:
        levels = [0 if i == e0_index else lvl + 1 for i, lvl in enumerate(levels)]
    syn_degree = [0] * n
    for st in syn_result.fundamental:
        syn_degree[st.i] += 1
        syn_degree[st.j] += 1
    comp_out = [0] * n
    comp_in = [0] * n
    for fc in comp_result.fundamental:
        comp_out[fc.prod_idx] += 1
        comp_in[fc.cons_idx] += 1
    n_pairs = n * (n - 1) // 2
    n_comp_pairs = sum(len(hier.ancestors[i]) for i in range(n))
    hasse_edges = sum(len(hier.parents[i]) for i in range(n))
    if e0_index is not None:
        hasse_edges += len(e0_cover_targets(hier, n, e0_index))
        n_comp_pairs += n - 1   # E0 is below every other ERC
    return HierarchyStats(
        n_ercs=n,
        n_persistent=sum(1 for e in ercs if e.is_persistent()),
        sizes=sizes, levels=levels, n_levels=(max(levels) + 1 if levels else 0),
        syn_degree=syn_degree, comp_out_degree=comp_out, comp_in_degree=comp_in,
        n_fundamental_syn=len(syn_result.fundamental),
        n_fundamental_comp=len(comp_result.fundamental),
        hasse_edges=hasse_edges, n_comparable_pairs=n_comp_pairs,
        n_incomparable_pairs=n_pairs - n_comp_pairs,
    )


def _new_pyvis_network(height="800px", *, level_separation=140, node_spacing=110):
    from pyvis.network import Network
    net = Network(height=height, width="100%", directed=True, notebook=False, bgcolor="#ffffff")
    net.set_options(f"""
    {{
      "layout": {{"hierarchical": {{"enabled": true, "direction": "DU",
                                    "sortMethod": "hubsize",
                                    "levelSeparation": {level_separation}, "nodeSpacing": {node_spacing}}}}},
      "physics": {{"enabled": false}},
      "interaction": {{"hover": true, "navigationButtons": true, "keyboard": true}}
    }}
    """)
    return net


def _inject_title(html_path: str, title: str) -> None:
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    banner = (f'<div style="font-family:sans-serif;padding:10px 16px;'
              f'background:#f4f4f4;border-bottom:1px solid #ddd;font-size:16px;">'
              f'<b>{title}</b></div>')
    content = content.replace("<body>", "<body>" + banner, 1)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(content)


def _histogram_svg(counts: dict, color: str, noun: str, xlabel: str) -> str:
    """Single-series bar chart of {label: count}, inline SVG."""
    from html import escape
    keys = list(counts.keys())
    if not keys:
        return f"<p style='color:#666'>No fundamental {escape(noun)}.</p>"
    vmax = max(counts.values()) or 1
    bar_w, slot = 18, 40
    left, top, plot_h, bottom = 34, 18, 150, 34
    width = left + len(keys) * slot + 10
    height = top + plot_h + bottom

    def _bar(x, v, color, tip):
        if v == 0:
            return ""
        h = max(2.0, plot_h * v / vmax)
        y0, y1 = top + plot_h, top + plot_h - h
        r = min(4.0, h / 2, bar_w / 2)
        # Rounded top corners only; flat on the baseline.
        d = (f"M{x},{y0} L{x},{y1 + r} Q{x},{y1} {x + r},{y1} L{x + bar_w - r},{y1} "
             f"Q{x + bar_w},{y1} {x + bar_w},{y1 + r} L{x + bar_w},{y0} Z")
        return (f'<path d="{d}" fill="{color}"><title>{escape(tip)}</title></path>'
                f'<text x="{x + bar_w / 2}" y="{y1 - 3}" font-size="10" fill="#333" '
                f'text-anchor="middle">{v}</text>')

    parts = [f'<svg width="{width}" height="{height}" font-family="sans-serif" role="img" '
             f'aria-label="Fundamental {escape(noun)} per {escape(xlabel)}">']
    for frac in (0.0, 0.5, 1.0):   # recessive gridlines
        y = top + plot_h - plot_h * frac
        parts.append(f'<line x1="{left}" x2="{width - 6}" y1="{y}" y2="{y}" stroke="#e5e5e5"/>'
                     f'<text x="{left - 5}" y="{y + 3}" font-size="10" fill="#666" '
                     f'text-anchor="end">{round(vmax * frac)}</text>')
    for g, lbl in enumerate(keys):
        x = left + (slot - bar_w) / 2 + g * slot
        parts.append(_bar(x, counts[lbl], color, f"{lbl}: {counts[lbl]} fundamental {noun}"))
        parts.append(f'<text x="{x + bar_w / 2}" y="{top + plot_h + 14}" font-size="10" '
                     f'fill="#333" text-anchor="middle">{escape(lbl)}</text>')
    parts.append(f'<text x="{left + (width - left) / 2}" y="{height - 4}" font-size="11" '
                 f'fill="#666" text-anchor="middle">{escape(xlabel)}</text>')
    parts.append("</svg>")
    return "".join(parts)


def _inject_stats_panel(html_path: str, stats: HierarchyStats, level_counts: dict) -> None:
    """Insert a summary table + per-level histograms/tables above the graph."""
    from html import escape
    td = 'style="padding:2px 10px;border-bottom:1px solid #eee"'
    tdr = 'style="padding:2px 10px;border-bottom:1px solid #eee;text-align:right"'
    summary_rows = [
        ("ERCs", stats.n_ercs),
        ("persistent ERCs", stats.n_persistent),
        ("hierarchy levels", stats.n_levels),
        ("containment (Hasse) edges", stats.hasse_edges),
        ("comparable ERC pairs", stats.n_comparable_pairs),
        ("incomparable ERC pairs", stats.n_incomparable_pairs),
        ("fundamental synergies", stats.n_fundamental_syn),
        ("fundamental complementarities", stats.n_fundamental_comp),
    ]
    summary = "".join(f"<tr><td {td}>{k}</td><td {tdr}>{v}</td></tr>" for k, v in summary_rows)
    syn_counts = {syn_triad_label(k, "→"): v for k, v in level_counts["syn"].items()}
    comp_counts = {f"({a},{b})": v for (a, b), v in level_counts["comp"].items()}

    def _count_table(counts, head):
        rows = "".join(f"<tr><td {td}>{escape(k)}</td><td {tdr}>{v}</td></tr>"
                       for k, v in counts.items())
        return ('<div style="max-height:230px;overflow-y:auto"><table style="border-collapse:collapse">'
                f'<tr><th {td}>{head}</th><th {tdr}>n</th></tr>{rows}</table></div>')

    def _block(heading, color, counts, noun, xlabel, head):
        swatch = (f'<span style="display:inline-block;width:10px;height:10px;border-radius:2px;'
                  f'background:{color};margin-right:6px"></span>')
        return ('<div style="display:flex;gap:14px;align-items:flex-start;max-width:100%">'
                '<div style="max-width:100%;overflow-x:auto">'
                f'<div style="font-weight:bold;margin-bottom:4px">{swatch}{heading}</div>'
                f'{_histogram_svg(counts, color, noun, xlabel)}</div>'
                f'{_count_table(counts, head)}</div>')

    panel = (
        '<details open style="font-family:sans-serif;font-size:13px;padding:8px 16px;'
        'border-bottom:1px solid #ddd;background:#fafafa">'
        '<summary style="cursor:pointer;font-weight:bold">Fundamental hierarchy statistics</summary>'
        '<div style="display:flex;flex-wrap:wrap;gap:28px;align-items:flex-start;margin-top:8px">'
        f'<table style="border-collapse:collapse">{summary}</table>'
        + _block("Fundamental synergies per level triad", COLOR_SYNERGY, syn_counts,
                 "synergies", "levels: contributor + contributor → target", "triad")
        + _block("Fundamental complementarities per level pair", COLOR_COMPLEMENTARITY,
                 comp_counts, "complementarities", "level pair (unordered)", "pair")
        + '</div>'
        '<div style="color:#666;font-size:11px;margin-top:6px">Level 0 = bottom of the hierarchy '
        '(the inflow ERC E0, when present). A synergy E<sub>i</sub> + E<sub>j</sub> &rarr; '
        'E<sub>k</sub> is binned as level(i)+level(j)&rarr;level(k), contributors sorted; '
        'complementarities by the unordered producer/consumer level pair.</div>'
        '</details>'
    )
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    # Goes right after the title banner, which _inject_title placed after <body>.
    marker = "</b></div>"
    idx = content.find(marker)
    if idx >= 0:
        idx += len(marker)
        content = content[:idx] + panel + content[idx:]
    else:
        content = content.replace("<body>", "<body>" + panel, 1)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(content)


# ---------------------------------------------------------------------------
# Click-to-inspect details (ERC contents, synergy generators, comp. reactions)
# ---------------------------------------------------------------------------

_MAX_LIST = 60   # cap long reaction lists in the inspector


def _closure_reactions(rn_data, mask: int) -> list:
    """R_E: reactions whose (E0-stripped) support lies inside `mask`."""
    return [r for r in range(rn_data.n_reactions) if rn_data.supp_q[r] & ~mask == 0]


def _erc_detail_data(rn_data, erc) -> dict:
    """Species-level contents of one ERC, all as bitmasks over raw species."""
    rxns = _closure_reactions(rn_data, erc.species_mask)
    supp = prod = 0
    for r in rxns:
        supp |= rn_data.supp_raw[r]
        prod |= rn_data.prod_raw[r]
    return {"rxns": rxns, "own": list(erc.reaction_indices),
            "supp": supp, "prod": prod, "req": erc.req_mask}


def build_inspector_details(ercs, hier, syn_result, comp_result, rn_data, stats,
                            *, e0_index=None, node_set=None) -> dict:
    """
    HTML snippets, keyed by pyvis node/edge id, shown when that element is
    clicked:

      ERC node "i"        : species, minimal generators, reactions (own and
                            all of R_E), supp / req / prod
      synergy "syn_i_j_k" : which minimal generator of the target E_k is
                            assembled from E_i's and E_j's species, and
                            which of E_k's reactions that unlocks
      complementarity edge "comp_p_c_s" : the species, the producer's
                            reactions making it, the consumer's needing it
      containment edge "cont_c_p" : child/parent species, what the parent
                            adds (species, reactions, prod) and how req changes

    Species sets include the inflow closure E0 where it is physically present
    (species, supp, prod); req and generators are E0-stripped, since inflow
    species are never required and never need to be generated.
    """
    from html import escape
    if node_set is None:
        node_set = set(range(len(ercs)))
    e0 = rn_data.E0_mask
    e0_names = ", ".join(rn_data.bitset_to_names(e0))
    data = {i: _erc_detail_data(rn_data, ercs[i]) for i in range(len(ercs))}

    def sp(mask):
        names = rn_data.bitset_to_names(mask)
        return "{" + ", ".join(escape(s) for s in names) + "}" if names else "&empty;"

    def rxn(r):
        return (f"<b>{escape(rn_data.reaction_names[r])}</b>: "
                f"{' + '.join(escape(s) for s in rn_data.bitset_to_names(rn_data.supp_raw[r])) or '&empty;'}"
                f" &rarr; {' + '.join(escape(s) for s in rn_data.bitset_to_names(rn_data.prod_raw[r])) or '&empty;'}")

    def rxn_list(rs):
        items = "".join(f"<li>{rxn(r)}</li>" for r in rs[:_MAX_LIST])
        more = f"<li><i>... {len(rs) - _MAX_LIST} more</i></li>" if len(rs) > _MAX_LIST else ""
        return f"<ul>{items}{more}</ul>"

    def erc_name(i):
        return f"E{i}" + (" (inflow)" if i == e0_index else "") + f" [L{stats.levels[i]}]"

    out: dict = {}
    for i in node_set:
        e, d = ercs[i], data[i]
        gens = ("&empty; (inflow: needs no reactant)" if not e.min_bases
                else "<ul>" + "".join(f"<li>{sp(b)}</li>" for b in e.min_bases) + "</ul>")
        parents = ", ".join(f"E{p}" for p in hier.parents[i]) or "-"
        children = ", ".join(f"E{c}" for c in hier.children[i]) or "-"
        out[str(i)] = (
            f"<h3>{erc_name(i)}</h3>"
            f"<p>{'persistent (req = &empty;)' if e.is_persistent() else 'not persistent'}"
            f" &middot; parents: {parents} &middot; children: {children}</p>"
            f"<p><b>species</b> ({bin(e.species_mask | e0).count('1')}): {sp(e.species_mask | e0)}</p>"
            f"<p><b>minimal generators</b> (smallest species sets whose closure is this ERC; "
            f"inflow {{{escape(e0_names)}}} always implicit):</p>{gens}"
            f"<p><b>supp</b> (consumed by R<sub>E</sub>): {sp(d['supp'])}</p>"
            f"<p><b>req</b> (consumed, not produced &mdash; needed from outside): {sp(d['req'])}</p>"
            f"<p><b>prod</b> (produced by R<sub>E</sub>): {sp(d['prod'])}</p>"
            f"<p><b>own reactions</b> (closure of support = this ERC; {len(d['own'])}):</p>{rxn_list(d['own'])}"
            f"<p><b>all reactions R<sub>E</sub></b> (incl. sub-ERCs; {len(d['rxns'])}):</p>{rxn_list(d['rxns'])}"
        )

    cover_edges = [(c, p) for c in node_set for p in hier.parents[c] if p in node_set]
    if e0_index is not None and e0_index in node_set:
        cover_edges += [(e0_index, p) for p in e0_cover_targets(hier, len(ercs), e0_index)
                        if p in node_set]
    for c, p in cover_edges:
        dc, dp = data[c], data[p]
        mc, mp = ercs[c].species_mask | e0, ercs[p].species_mask | e0
        gained = sorted(set(dp["rxns"]) - set(dc["rxns"]))
        out[f"cont_{c}_{p}"] = (
            f"<h3>Containment {erc_name(c)} &sub; {erc_name(p)}</h3>"
            f"<p>{'Implicit edge: the inflow closure is contained in every ERC. ' if c == e0_index else ''}"
            f"Direct cover (no ERC strictly in between).</p>"
            f"<p><b>E{c} species</b>: {sp(mc)}<br><b>E{p} species</b>: {sp(mp)}</p>"
            f"<p><b>species added</b> (E{p} &setminus; E{c}): {sp(mp & ~mc)}</p>"
            f"<p><b>req</b>: E{c} {sp(ercs[c].req_mask)} &rarr; E{p} {sp(ercs[p].req_mask)}<br>"
            f"&nbsp;&nbsp;satisfied inside E{p}: {sp(ercs[c].req_mask & ~ercs[p].req_mask)}<br>"
            f"&nbsp;&nbsp;newly required: {sp(ercs[p].req_mask & ~ercs[c].req_mask)}</p>"
            f"<p><b>prod added</b>: {sp(dp['prod'] & ~dc['prod'])}</p>"
            f"<p><b>E{p} minimal generators</b>: "
            + (", ".join(sp(b) for b in ercs[p].min_bases) or "&empty;") + "</p>"
            f"<p><b>reactions added</b> (R<sub>E{p}</sub> &setminus; R<sub>E{c}</sub>; {len(gained)}):</p>"
            f"{rxn_list(gained)}"
        )

    for (i, j, k) in synergy_triads(syn_result, node_set):
        mi, mj = ercs[i].species_mask, ercs[j].species_mask
        both = mi | mj
        gen_rows = []
        for b in ercs[k].min_bases:
            if b & ~both == 0 and b & ~mi and b & ~mj:   # needs both ERCs
                rs = [r for r in ercs[k].reaction_indices
                      if rn_data.supp_q[r] & ~both == 0 and rn_data.supp_q[r] & ~b == 0]
                gen_rows.append(
                    f"<li>generator {sp(b)}<br>"
                    f"&nbsp;&nbsp;from E{i}: {sp(b & mi)}<br>"
                    f"&nbsp;&nbsp;from E{j}: {sp(b & mj)}"
                    + (f"{rxn_list(rs)}" if rs else "") + "</li>")
        lv = stats.levels
        out[f"syn_{i}_{j}_{k}"] = (
            f"<h3>Fundamental synergy {erc_name(i)} + {erc_name(j)} &rarr; {erc_name(k)}</h3>"
            f"<p>level triad: {escape(syn_triad_label((*sorted((lv[i], lv[j])), lv[k]), chr(0x2192)))}</p>"
            f"<p>E{i}: {sp(mi | e0)}<br>E{j}: {sp(mj | e0)}<br>"
            f"<b>target E{k}</b>: {sp(ercs[k].species_mask | e0)}</p>"
            f"<p>Neither ERC alone contains a minimal generator of E{k}; together they do. "
            f"Each generator is split by which ERC supplies each species (a species in both "
            f"appears on both sides), followed by the target reactions it fires.</p>"
            f"<ul>{''.join(gen_rows) or '<li><i>no single generator split found</i></li>'}</ul>")

    for fc in comp_result.fundamental:
        p, c, s = fc.prod_idx, fc.cons_idx, fc.species
        if p not in node_set or c not in node_set:
            continue
        bit = 1 << s
        making = [r for r in data[p]["rxns"] if rn_data.prod_raw[r] & bit]
        needing = [r for r in data[c]["rxns"] if rn_data.supp_raw[r] & bit]
        out[f"comp_{p}_{c}_{s}"] = (
            f"<h3>Fundamental complementarity</h3>"
            f"<p>{erc_name(p)} supplies <b>{escape(rn_data.species_names[s])}</b> to {erc_name(c)}"
            f"{' (chain: hierarchy-comparable pair)' if fc.chain else ''}</p>"
            f"<p>E{c} req: {sp(ercs[c].req_mask)}<br>E{p} prod: {sp(data[p]['prod'])}</p>"
            f"<p><b>reactions of E{p} producing it</b> ({len(making)}):</p>{rxn_list(making)}"
            f"<p><b>reactions of E{c} consuming it</b> ({len(needing)}):</p>{rxn_list(needing)}"
        )
    return out


def _inject_inspector(html_path: str, details: dict) -> None:
    """Add a right-hand panel filled on node/edge click from `details`."""
    import json
    payload = json.dumps(details).replace("</", "<\\/")
    block = (
        '<style>'
        '#inspector{position:fixed;top:60px;right:12px;width:420px;max-width:calc(100vw - 32px);'
        'max-height:calc(100vh - 80px);overflow:auto;background:#fff;border:1px solid #ccc;'
        'border-radius:6px;box-shadow:0 2px 10px rgba(0,0,0,.15);padding:10px 14px;'
        'font-family:sans-serif;font-size:12.5px;line-height:1.4;z-index:1000}'
        '#inspector h3{margin:4px 0 8px;font-size:14px}#inspector p{margin:6px 0}'
        '#inspector ul{margin:2px 0 6px;padding-left:18px}'
        '#inspector .close{float:right;cursor:pointer;color:#888;font-size:16px}'
        '.vis-tooltip{white-space:pre-line}'
        '</style>'
        '<div id="inspector" hidden><span class="close" title="close">&times;</span>'
        '<div id="inspector-body"></div></div>'
        '<script>'
        f'var INSPECTOR_DETAILS = {payload};'
        '(function(){'
        ' var panel=document.getElementById("inspector"), body=document.getElementById("inspector-body");'
        ' panel.querySelector(".close").onclick=function(){panel.hidden=true;};'
        ' function show(id){var h=INSPECTOR_DETAILS[String(id)]; if(!h) return false;'
        '   body.innerHTML=h; panel.hidden=false; panel.scrollTop=0; return true;}'
        ' function hook(){ if(typeof network==="undefined"||!network){setTimeout(hook,100);return;}'
        '   network.on("click",function(p){'
        '     if(p.nodes.length&&show(p.nodes[0])) return;'
        '     if(p.edges.length){for(var k=0;k<p.edges.length;k++){if(show(p.edges[k])) return;}}'
        '     if(!p.nodes.length&&!p.edges.length) panel.hidden=true;'
        '   });}'
        ' hook();'
        '})();'
        '</script>'
    )
    with open(html_path, "r", encoding="utf-8") as f:
        content = f.read()
    content = content.replace("</body>", block + "</body>", 1)
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(content)


def plot_erc_hierarchy(
    ercs, hier, syn_result, comp_result, out_path: str,
    *,
    stats: HierarchyStats | None = None,
    highlight_ercs=None,
    highlight_label: str = "highlighted",
    max_nodes: int = 250,
    show_synergy: bool = True,
    show_complementarity: bool = True,
    title: str = "ERC hierarchy",
    e0_index: int | None = None,
    species_names=None,
    stats_panel: bool = False,
    rn_data=None,
) -> str:
    """
    Render the ERC containment hierarchy with fundamental synergies and
    complementarities overlaid, as an interactive pyvis HTML graph.

    Nodes are leveled by containment depth (subset-minimal ERCs at the
    bottom). Gray edges are containment (child -> parent). Orange diamond
    nodes mark fundamental synergies (edges to both contributing ERCs).
    Blue dashed edges are fundamental complementarities (producer -> consumer).

    Parameters
    ----------
    ercs, hier, syn_result, comp_result : outputs of
        pyCOT.analysis.organizations.compute_ercs / build_hierarchy /
        compute_synergies_basis_first / compute_complementarities
    out_path : where to write the HTML file (parent dirs created if needed)
    stats : precomputed HierarchyStats; computed automatically if omitted
    highlight_ercs : iterable of ERC indices to color distinctly (e.g. the
        members of one particular EPM/ESPM), so their position in the full
        hierarchy is immediately visible
    max_nodes : caps total nodes shown for readability on genome-scale
        networks -- keeps the highlighted set plus a level-balanced sample
    show_synergy, show_complementarity : toggle each relation overlay
    title : banner text written into the HTML page
    e0_index : index of the inflow ERC E0, if any -- drawn as the bottom of
        the hierarchy with containment edges to every subset-minimal ERC
        (see compute_hierarchy_stats)
    species_names : list mapping species bit index -> name, used to label
        complementarity tooltips (falls back to the bit index)
    stats_panel : add a collapsible panel above the graph with summary
        counts and a histogram/table of fundamental relations per level pair
    rn_data : the compiled RNData the ERCs came from. When given, clicking
        an ERC / synergy diamond / complementarity edge opens a side panel
        with its reactions, minimal generators and supp/req/prod sets, and
        for synergies how each target generator splits across the pair
        (see build_inspector_details). Also supplies species_names.

    Returns
    -------
    str : out_path, for chaining
    """
    if stats is None:
        stats = compute_hierarchy_stats(ercs, hier, syn_result, comp_result, e0_index=e0_index)
    if rn_data is not None and species_names is None:
        species_names = rn_data.species_names
    click_hint = "\n(click for details)" if rn_data is not None else ""

    n = len(ercs)
    highlight = set(highlight_ercs or [])

    if n > max_nodes:
        by_level: dict = {}
        for i in range(n):
            by_level.setdefault(stats.levels[i], []).append(i)
        keep = set(highlight)
        budget = max_nodes - len(keep)
        levels_sorted = sorted(by_level.keys())
        per_level = max(1, budget // max(1, len(levels_sorted)))
        for lvl in levels_sorted:
            for i in by_level[lvl]:
                if len(keep) >= max_nodes:
                    break
                if i not in keep and len([x for x in keep if stats.levels[x] == lvl]) < per_level:
                    keep.add(i)
        node_set = keep
        capped_note = f" (showing {len(node_set)} of {n} ERCs, level-balanced sample)"
    else:
        node_set = set(range(n))
        capped_note = ""

    net = _new_pyvis_network()
    for i in node_set:
        is_hl = i in highlight
        color = COLOR_HIGHLIGHT if is_hl else COLOR_NEUTRAL
        size = 22 if is_hl else 14
        # pyvis 0.3.2 / vis-network 9 render `title` as plain text, so use
        # newlines (shown via .vis-tooltip{white-space:pre-line}), not <br>.
        title_txt = (f"E{i}  |  {stats.sizes[i]} species  |  level {stats.levels[i]}\n"
                     f"synergy degree: {stats.syn_degree[i]}\n"
                     f"complementarity out/in: {stats.comp_out_degree[i]}/{stats.comp_in_degree[i]}")
        if is_hl:
            title_txt += f"\n{highlight_label}"
        title_txt += click_hint
        label = f"E{i} (inflow)" if i == e0_index else f"E{i}"
        net.add_node(i, label=label, level=2 * stats.levels[i], color=color, size=size,
                     title=title_txt, shape="dot")

    for i in node_set:
        for p in hier.parents[i]:
            if p in node_set:
                net.add_edge(i, p, color=COLOR_CONTAINMENT, width=1.2, arrows="to",
                             id=f"cont_{i}_{p}", title=f"E{i} ⊂ E{p}" + click_hint)
    if e0_index is not None and e0_index in node_set:
        for p in e0_cover_targets(hier, n, e0_index):
            if p in node_set:
                net.add_edge(e0_index, p, color=COLOR_CONTAINMENT, width=1.2, arrows="to",
                             id=f"cont_{e0_index}_{p}", title=f"E{e0_index} ⊂ E{p}" + click_hint)

    if show_synergy:
        # One diamond per synergy triad: E_i, E_j -- diamond --> E_k. ERCs sit
        # on even layout rows (2*level); a diamond sits on the odd row just
        # above its higher contributor, so it never shares a row with an ERC.
        for (i, j, k) in synergy_triads(syn_result, node_set):
            d_id = f"syn_{i}_{j}_{k}"
            net.add_node(d_id, label="+", shape="diamond", color=COLOR_SYNERGY, size=10,
                         level=2 * max(stats.levels[i], stats.levels[j]) + 1,
                         title=f"fundamental synergy E{i} + E{j} -> E{k}" + click_hint)
            net.add_edge(i, d_id, color=COLOR_SYNERGY, width=1, arrows="")
            net.add_edge(j, d_id, color=COLOR_SYNERGY, width=1, arrows="")
            net.add_edge(d_id, k, color=COLOR_SYNERGY, width=1.6, arrows="to")

    if show_complementarity:
        for fc in comp_result.fundamental:
            if fc.prod_idx in node_set and fc.cons_idx in node_set:
                net.add_edge(fc.prod_idx, fc.cons_idx, color=COLOR_COMPLEMENTARITY, width=1,
                             dashes=True, arrows="to",
                             id=f"comp_{fc.prod_idx}_{fc.cons_idx}_{fc.species}",
                             title=f"E{fc.prod_idx} supplies "
                                   f"{species_names[fc.species] if species_names else 'species ' + str(fc.species)}"
                                   f" to E{fc.cons_idx}" + click_hint)

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    net.write_html(out_path, open_browser=False, notebook=False)
    _inject_title(out_path, title + capped_note)
    if stats_panel:
        _inject_stats_panel(out_path, stats,
                            relation_level_counts(stats, syn_result, comp_result))
    if rn_data is not None:
        _inject_inspector(out_path, build_inspector_details(
            ercs, hier, syn_result, comp_result, rn_data, stats,
            e0_index=e0_index, node_set=node_set))
    return out_path
