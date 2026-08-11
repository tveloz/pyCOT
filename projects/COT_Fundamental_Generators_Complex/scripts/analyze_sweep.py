#!/usr/bin/env python
"""
analyze_sweep.py

Post-processing / visualization for sweep_results.csv (produced by
run_bigg_sweep.py) and its deadend_economics/ sidecar files. Read-only:
never touches the sweep itself, safe to run at any point mid-sweep (it
just analyzes whatever rows exist so far).

Produces, in analysis_output/:
  summary_table.csv          -- one row per network, the key comparison
                                 columns side by side (baseline vs adaptive
                                 vs original-ESPM), sorted by syn_density
  speedup_vs_density.png     -- does synergy density predict the
                                 propagation speedup/slowdown?
  walltime_vs_size.png       -- does network size (n_ercs) predict cost?
  deadend_gini_vs_density.png-- does density predict bottleneck concentration
                                 ("critical needs" hypothesis)?
  scc_vs_density.png         -- structural cross-check: does density track
                                 the earlier entanglement (SCC) finding?
  epm_espm_by_order.png      -- EPM (order 0) vs ESPM (order >=1) counts per
                                 network, stacked by order -- only for
                                 networks where the epm_espm_original stage
                                 completed without timing out
  top_blocking_species.png   -- species that recur as dead-end bottlenecks
                                 across the MOST networks (aggregated from
                                 every deadend_economics/*.csv sidecar) --
                                 the cross-organism "critical needs" view

Usage (run from the pyCOT repository root):
    python projects/COT_Fundamental_Generators_Complex/scripts/analyze_sweep.py
"""
import os
import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # no display available; always save to file
import matplotlib.pyplot as plt

_here = os.path.dirname(os.path.abspath(__file__))
SWEEP_CSV = os.path.join(_here, 'sweep_results.csv')
SIDECAR_DIR = os.path.join(_here, 'deadend_economics')
OUT_DIR = os.path.join(_here, 'analysis_output')


def load_sweep():
    if not os.path.exists(SWEEP_CSV):
        raise SystemExit(f'No {SWEEP_CSV} found -- run run_bigg_sweep.py first.')
    df = pd.read_csv(SWEEP_CSV)
    if df.empty:
        raise SystemExit(f'{SWEEP_CSV} is empty -- no results to analyze yet.')

    # Defensive: if this CSV predates a schema change (e.g. the ESPM columns
    # added later), older rows just won't have those columns at all rather
    # than having them blank -- add them as all-NaN so every plot/column
    # reference below works uniformly regardless of which schema vintage
    # produced the file.
    RISK_COLS = [
        'riskepm_n_epms',
        'riskepm_syn_safe', 'riskepm_syn_mid', 'riskepm_syn_risky',
        'riskepm_comp_safe', 'riskepm_comp_mid', 'riskepm_comp_risky',
        'riskepm_states_syn_safe', 'riskepm_states_syn_mid', 'riskepm_states_syn_risky',
        'riskepm_states_comp_safe', 'riskepm_states_comp_mid', 'riskepm_states_comp_risky',
        'riskepm_ssm_syn_safe', 'riskepm_ssm_syn_mid', 'riskepm_ssm_syn_risky',
        'riskepm_ssm_comp_safe', 'riskepm_ssm_comp_mid', 'riskepm_ssm_comp_risky',
        'riskepm_deadend_syn_safe', 'riskepm_deadend_syn_mid', 'riskepm_deadend_syn_risky',
        'riskepm_deadend_comp_safe', 'riskepm_deadend_comp_mid', 'riskepm_deadend_comp_risky',
        'riskespm_total', 'riskespm_max_order',
        'riskespm_syn_safe', 'riskespm_syn_mid', 'riskespm_syn_risky',
        'riskespm_comp_safe', 'riskespm_comp_mid', 'riskespm_comp_risky',
        'riskespm_hierarchy_unclassified',
        'riskespm_states_syn_safe', 'riskespm_states_syn_mid', 'riskespm_states_syn_risky',
        'riskespm_states_comp_safe', 'riskespm_states_comp_mid', 'riskespm_states_comp_risky',
        'riskespm_ssm_syn_safe', 'riskespm_ssm_syn_mid', 'riskespm_ssm_syn_risky',
        'riskespm_ssm_comp_safe', 'riskespm_ssm_comp_mid', 'riskespm_ssm_comp_risky',
    ]

    expected_cols = [
        'espm_n_epm', 'espm_total', 'espm_max_order', 'espm_by_order', 'espm_complete',
        'dead_end_req_gini', 'dead_end_req_top1_share', 'dead_end_req_top10pct_share',
        'n_species_ever_blocking', 'progress_fraction', 'max_scc_fraction',
        'syn_density', 'comp_density', 'syn_to_comp_ratio',
    ] + RISK_COLS
    for col in expected_cols:
        if col not in df.columns:
            df[col] = pd.NA

    # numeric coercion: blank cells (e.g. baseline-only columns on an
    # epm_espm_original row) become NaN rather than staying strings
    numeric_cols = [
        'n_species', 'n_reactions', 'n_ercs', 'n_fund_synergies', 'n_fund_complementarities',
        'syn_density', 'comp_density', 'syn_to_comp_ratio', 'max_syn_outdegree',
        'max_scc_fraction', 'n_nontrivial_sccs', 'wall_time_s', 'progress_fraction',
        'states_explored', 'ssm_found', 'leaves_found', 'n_epms',
        'dead_end_req_gini', 'dead_end_req_top1_share', 'dead_end_req_top10pct_share',
        'n_species_ever_blocking', 'espm_n_epm', 'espm_total', 'espm_max_order',
    ] + RISK_COLS
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')

    # 'timed_out'/'error' round-trip through CSV as strings ("True"/"False"/
    # blank) -- normalize explicitly rather than relying on pandas' type
    # inference, which can vary with how much NaN is mixed in per column.
    if 'timed_out' in df.columns:
        df['timed_out'] = df['timed_out'].map({'True': True, 'False': False, True: True, False: False})
    if 'espm_complete' in df.columns:
        df['espm_complete'] = df['espm_complete'].map({'True': True, 'False': False, True: True, False: False})
    if 'error' in df.columns:
        df['error'] = df['error'].replace('', pd.NA)

    # Resuming means the SAME (network, config) can now appear multiple times
    # across separate sweep invocations (a paused/timed-out attempt, then a
    # later one that picks up from its checkpoint and gets further). Each
    # later row is a strict continuation -- same or more progress -- of the
    # one before it, so keeping only the chronologically last row per
    # (network, config) always reflects the most-advanced known state.
    # Without this, downstream merges on 'network' would silently row-explode
    # (cartesian product) once duplicates exist.
    if 'timestamp' in df.columns:
        df = df.sort_values('timestamp').drop_duplicates(['network', 'config'], keep='last')
        df = df.reset_index(drop=True)

    return df


def build_summary_table(df):
    """One row per network: structure + baseline vs adaptive timing + ESPM summary."""
    struct_cols = ['network', 'n_species', 'n_reactions', 'n_ercs', 'syn_density',
                   'comp_density', 'syn_to_comp_ratio', 'max_scc_fraction']
    struct = df[df.config == 'baseline'][struct_cols].drop_duplicates('network')

    base = df[df.config == 'baseline'][['network', 'wall_time_s', 'timed_out',
                                         'progress_fraction', 'dead_end_req_gini',
                                         'n_species_ever_blocking']]
    base = base.rename(columns={'wall_time_s': 'baseline_wall_s', 'timed_out': 'baseline_timed_out',
                                 'progress_fraction': 'baseline_progress'})

    adap = df[df.config == 'adaptive'][['network', 'wall_time_s', 'timed_out', 'progress_fraction']]
    adap = adap.rename(columns={'wall_time_s': 'adaptive_wall_s', 'timed_out': 'adaptive_timed_out',
                                 'progress_fraction': 'adaptive_progress'})

    espm = df[df.config == 'epm_espm_original'][['network', 'timed_out', 'espm_complete', 'espm_n_epm',
                                                   'espm_total', 'espm_max_order']]
    espm = espm.rename(columns={'timed_out': 'espm_timed_out'})

    out = struct.merge(base, on='network', how='left') \
                .merge(adap, on='network', how='left') \
                .merge(espm, on='network', how='left')

    # speedup only meaningful when NEITHER side timed out (a timed-out run's
    # wall_time_s is just the budget, not a real completion time)
    both_finished = (~out['baseline_timed_out'].fillna(True)) & (~out['adaptive_timed_out'].fillna(True))
    out['speedup'] = pd.NA
    out.loc[both_finished, 'speedup'] = (
        out.loc[both_finished, 'baseline_wall_s'] / out.loc[both_finished, 'adaptive_wall_s'].replace(0, pd.NA)
    )

    return out.sort_values('syn_density')


def plot_speedup_vs_density(summary, out_dir):
    """
    Compares two search methods on the same networks: 'baseline' vs.
    'adaptive' (adaptive adds an extra shortcut -- forced-move detection --
    on top of baseline). Question: does that shortcut help more on networks
    where reaction-classes are more richly interconnected?
    """
    d = summary.dropna(subset=['speedup'])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(d['syn_density'], d['speedup'], s=40)
    for _, r in d.iterrows():
        ax.annotate(r['network'], (r['syn_density'], r['speedup']), fontsize=7,
                    xytext=(3, 3), textcoords='offset points')
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=1, label='same speed (no benefit)')
    ax.set_xlabel('network interconnection density\n(avg. number of "synergy" partners per reaction-class)')
    ax.set_ylabel('speed-up from the extra shortcut\n(> 1 = faster,  < 1 = slower,  1 = no difference)')
    ax.set_title('Does the extra search shortcut pay off more on densely-connected networks?\n'
                  '(each dot = one network; only networks where both methods finished in time)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'speedup_vs_density.png'), dpi=150)
    plt.close(fig)


def plot_walltime_vs_size(summary, out_dir):
    """How fast does search time grow as networks get bigger?"""
    d = summary.dropna(subset=['baseline_wall_s'])
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = ['tab:red' if t else 'tab:blue' for t in d['baseline_timed_out'].fillna(False)]
    ax.scatter(d['n_ercs'], d['baseline_wall_s'], s=40, c=colors)
    for _, r in d.iterrows():
        ax.annotate(r['network'], (r['n_ercs'], r['baseline_wall_s']), fontsize=7,
                    xytext=(3, 3), textcoords='offset points')
    ax.set_yscale('log')
    ax.set_xlabel('network size (number of reaction-classes)')
    ax.set_ylabel('search time in seconds (log scale)')
    ax.set_title('How much does search time grow with network size?')
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:blue', markersize=8,
                           label='finished in time'),
               plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:red', markersize=8,
                           label='hit the time limit (true cost is even higher than shown)')]
    ax.legend(handles=handles, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'walltime_vs_size.png'), dpi=150)
    plt.close(fig)


def plot_deadend_gini_vs_density(summary, out_dir):
    """
    A 'dead end' is a search path that got stuck missing an ingredient it
    could never produce. Question: in each network, are dead ends caused by
    a few chronic "bottleneck" ingredients, or spread evenly across many
    different ones?
    """
    d = summary.dropna(subset=['dead_end_req_gini'])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(d['syn_density'], d['dead_end_req_gini'], s=40)
    for _, r in d.iterrows():
        ax.annotate(r['network'], (r['syn_density'], r['dead_end_req_gini']), fontsize=7,
                    xytext=(3, 3), textcoords='offset points')
    ax.set_xlabel('network interconnection density')
    ax.set_ylabel('how concentrated the dead-end causes are\n(0 = spread evenly across many ingredients,\n'
                  '1 = almost always the same one or two "bottleneck" ingredients)')
    ax.set_title('Do dead ends have a few chronic "bottleneck" causes, or many different ones?')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'deadend_gini_vs_density.png'), dpi=150)
    plt.close(fig)


def plot_scc_vs_density(summary, out_dir):
    """
    A tightly-knit cluster is a group of reaction-classes that all depend on
    each other in a loop (none can be resolved without the others). Cross-
    checks the density metric against this independent structural signal.
    """
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(summary['syn_density'], summary['max_scc_fraction'], s=40)
    for _, r in summary.iterrows():
        ax.annotate(r['network'], (r['syn_density'], r['max_scc_fraction']), fontsize=7,
                    xytext=(3, 3), textcoords='offset points')
    ax.set_xlabel('network interconnection density')
    ax.set_ylabel('fraction of the network in its largest\nmutually-dependent (looped) cluster')
    ax.set_title('Cross-check: does interconnection density line up with\nhaving one big tangle of mutually-dependent reaction-classes?')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'scc_vs_density.png'), dpi=150)
    plt.close(fig)


def parse_by_order(s):
    """'1:45;2:12' -> {1: 45, 2: 12}"""
    if not isinstance(s, str) or not s:
        return {}
    out = {}
    for part in s.split(';'):
        if ':' not in part:
            continue
        k, v = part.split(':')
        out[int(k)] = int(v)
    return out


def plot_epm_espm_by_order(df, out_dir):
    """
    Shows every network with any ESPM data so far -- both fully-completed
    BFS runs AND ones still mid-progress across checkpointed sweep
    invocations (espm_complete=False, e.g. a multi-night run still going).
    Partial networks are marked with a trailing '*' on their x-axis label
    and a lighter bar edge, since their true final counts are still growing.
    """
    espm_rows = df[(df.config == 'epm_espm_original') & (df.error.isna()) & df['espm_total'].notna()]
    if espm_rows.empty:
        print('  [skip] epm_espm_by_order.png -- no epm_espm_original runs with any data yet')
        return

    espm_rows = espm_rows.sort_values('espm_max_order' if 'espm_max_order' in espm_rows else 'network')
    is_complete = {r['network']: bool(r.get('espm_complete')) for _, r in espm_rows.iterrows()}
    labels = [n + ('' if is_complete[n] else '*') for n in espm_rows['network']]
    networks = espm_rows['network'].tolist()
    orders_seen = set()
    per_network_by_order = {}
    for _, r in espm_rows.iterrows():
        by_order = parse_by_order(r.get('espm_by_order', ''))
        by_order[0] = int(r['espm_n_epm']) if pd.notna(r.get('espm_n_epm')) else 0
        per_network_by_order[r['network']] = by_order
        orders_seen.update(by_order.keys())

    max_order_to_show = min(max(orders_seen, default=0), 12)   # cap for plot legibility
    orders = list(range(0, max_order_to_show + 1))

    fig, ax = plt.subplots(figsize=(max(8, len(networks) * 0.6), 6))
    bottom = [0] * len(networks)
    cmap = plt.get_cmap('viridis', len(orders))
    for oi, order in enumerate(orders):
        heights = [per_network_by_order[n].get(order, 0) for n in networks]
        ax.bar(labels, heights, bottom=bottom, color=cmap(oi),
               label=f'depth {order}' + (' (smallest, indivisible)' if order == 0 else ''))
        bottom = [b + h for b, h in zip(bottom, heights)]
    for label, n, total in zip(labels, networks, bottom):
        if not is_complete[n]:
            ax.annotate('still growing', (label, total), fontsize=7, color='tab:red',
                        ha='center', xytext=(0, 3), textcoords='offset points')
    ax.set_ylabel('number of self-sufficient modules found')
    ax.set_title('Self-sufficient modules found per network, by construction "depth"\n'
                  '(depth 0 = smallest, indivisible modules; depth 1+ = built by combining smaller ones -- '
                  'darker = more combining steps)\n'
                  '* = still in progress (this network has more to find; numbers will keep growing)')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8, title='construction depth')
    plt.setp(ax.get_xticklabels(), rotation=60, ha='right', fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'epm_espm_by_order.png'), dpi=150)
    plt.close(fig)


def plot_top_blocking_species(out_dir, top_n=25):
    """
    Aggregate every deadend_economics/*.csv sidecar: which species recur as
    dead-end bottlenecks across the MOST different networks? (Not just
    highest total block_count in one network -- a species that's the #1
    bottleneck in 10 different organisms is a more interesting "universal
    critical need" than one that's extremely common in only one.)
    """
    files = glob.glob(os.path.join(SIDECAR_DIR, '*_deadend_req_species.csv'))
    if not files:
        print('  [skip] top_blocking_species.png -- no deadend_economics/ sidecar files found')
        return

    frames = []
    for path in files:
        base = os.path.basename(path).replace('_deadend_req_species.csv', '')
        # base is '{network}_{config}' -- split off the trailing config token
        network = base.rsplit('_baseline', 1)[0] if base.endswith('_baseline') else \
                  base.rsplit('_adaptive', 1)[0] if base.endswith('_adaptive') else base
        d = pd.read_csv(path)
        d['network'] = network
        frames.append(d)
    all_rows = pd.concat(frames, ignore_index=True)

    agg = all_rows.groupby('species_name').agg(
        n_networks=('network', 'nunique'),
        total_block_count=('block_count', 'sum'),
    ).sort_values(['n_networks', 'total_block_count'], ascending=False).head(top_n)

    fig, ax = plt.subplots(figsize=(9, max(6, top_n * 0.3)))
    y_pos = range(len(agg))
    ax.barh(y_pos, agg['n_networks'])
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(agg.index, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel('number of different organisms (networks) where this ingredient caused a dead end')
    ax.set_title(f'Top {top_n} ingredients that repeatedly cause dead ends across organisms\n'
                  '(a "universal bottleneck": missing this ingredient stalls the search in many different networks)')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'top_blocking_species.png'), dpi=150)
    plt.close(fig)


# The search grows a module one piece at a time. Each addition is either:
#   safe  -- a sure bet, it can never create a new unmet need
#   mid   -- partly justified -- only half of what justifies the move is
#            already backed by what's been built so far
#   risky -- a gamble -- neither half is backed yet; it might still pay off,
#            or it might lead nowhere (a dead end)
CATS = ['safe', 'mid', 'risky']
CAT_COLORS = {'safe': 'tab:green', 'mid': 'tab:orange', 'risky': 'tab:red'}
CAT_LABELS = {'safe': 'safe (sure bet)', 'mid': 'mid (partly justified)', 'risky': 'risky (a gamble)'}


def _risk_rows(df):
    """Rows from the risk_analysis config that actually finished (no timeout/error)."""
    r = df[(df.config == 'risk_analysis') & df['timed_out'].fillna(True).eq(False) & df['error'].isna()]
    return r[r['riskepm_n_epms'].notna()]


def plot_risk_move_classification(df, out_dir):
    """
    Every step the search takes to grow a module is one of two kinds:
      - a "synergy" step: two existing pieces team up to unlock a third
      - a "complementarity" step: one piece supplies an ingredient another
        piece is still missing
    This shows, for each kind of step, what fraction were safe/mid/risky
    (see CATS comment above) -- separately for the two search stages:
      - "building blocks" = finding the smallest, indivisible modules
      - "combining blocks" = building bigger modules out of smaller ones
    """
    r = _risk_rows(df)
    if r.empty:
        print('  [skip] risk_move_classification.png -- no completed risk_analysis runs yet')
        return

    r = r.sort_values('network')
    networks = r['network'].tolist()
    panels = [
        ('Stage 1 (building blocks): synergy steps', 'riskepm_syn'),
        ('Stage 1 (building blocks): complementarity steps', 'riskepm_comp'),
        ('Stage 2 (combining blocks): synergy steps', 'riskespm_syn'),
        ('Stage 2 (combining blocks): complementarity steps', 'riskespm_comp'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(max(9, len(networks) * 1.3), 10), sharex=True)
    for ax, (title, prefix) in zip(axes.flat, panels):
        bottom = [0.0] * len(networks)
        for cat in CATS:
            col = f'{prefix}_{cat}'
            counts = r[col].fillna(0).to_numpy()
            totals = sum(r[f'{prefix}_{c}'].fillna(0).to_numpy() for c in CATS)
            totals = [t if t > 0 else 1 for t in totals]
            fracs = [c / t for c, t in zip(counts, totals)]
            ax.bar(networks, fracs, bottom=bottom, color=CAT_COLORS[cat], label=CAT_LABELS[cat])
            bottom = [b + fr for b, fr in zip(bottom, fracs)]
        ax.set_title(title, fontsize=10)
        ax.set_ylabel('share of steps taken')
        ax.set_ylim(0, 1.02)
        plt.setp(ax.get_xticklabels(), rotation=45, ha='right', fontsize=8)
    axes[0, 0].legend(loc='upper right', fontsize=8)
    fig.suptitle('What fraction of search steps were safe bets vs. gambles?\n(one bar per network; green+orange+red always add up to 100%)')
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'risk_move_classification.png'), dpi=150)
    plt.close(fig)


STAGES = [
    ('Stage 1: building the smallest modules', 'riskepm'),
    ('Stage 2: combining modules into bigger ones', 'riskespm'),
]
MOVE_TYPES = [('synergy steps', 'syn'), ('complementarity steps', 'comp')]


def plot_risk_contribution_vs_effort(df, out_dir):
    """
    THE key question this whole instrumentation exists to answer: does each
    category (safe/mid/risky) find its FAIR share of modules given how much
    search effort it costs? Solid bar = % of ALL modules found (in that
    stage) that came from this category AND this kind of step. Hatched bar =
    % of the total search work that category+kind used up. If the solid bar
    is taller than its hatched twin, that combination is "punching above its
    weight" -- finding more than its share of the effort would predict.
    One row per stage, one column per kind of step (synergy/complementarity)
    -- so this connects directly to the move-classification plot, which
    shows the same 2x2 breakdown for raw step counts.
    """
    r = _risk_rows(df)
    if r.empty:
        print('  [skip] risk_contribution_vs_effort.png -- no completed risk_analysis runs yet')
        return

    r = r.sort_values('network')
    networks = r['network'].tolist()

    fig, axes = plt.subplots(2, 2, figsize=(max(11, len(networks) * 2.0), 10))
    for row_i, (stage_title, stage_prefix) in enumerate(STAGES):
        # shared denominator across BOTH move types, so "share" always means
        # "share of this whole stage's modules/effort", not just within one type
        ssm_tot = sum(r[f'{stage_prefix}_ssm_{t}_{c}'].fillna(0) for _, t in MOVE_TYPES for c in CATS).replace(0, 1)
        states_tot = sum(r[f'{stage_prefix}_states_{t}_{c}'].fillna(0) for _, t in MOVE_TYPES for c in CATS).replace(0, 1)
        for col_i, (type_title, mtype) in enumerate(MOVE_TYPES):
            ax = axes[row_i, col_i]
            x = range(len(networks))
            width = 0.2
            for ci, cat in enumerate(CATS):
                contrib = (r[f'{stage_prefix}_ssm_{mtype}_{cat}'].fillna(0) / ssm_tot).to_numpy()
                effort = (r[f'{stage_prefix}_states_{mtype}_{cat}'].fillna(0) / states_tot).to_numpy()
                offset = (ci - 1) * width
                xs = [xi + offset for xi in x]
                ax.bar([xi - width / 4 for xi in xs], contrib, width=width / 2.2,
                       color=CAT_COLORS[cat], label=CAT_LABELS[cat] if (row_i == 0 and col_i == 0) else None)
                ax.bar([xi + width / 4 for xi in xs], effort, width=width / 2.2,
                       color=CAT_COLORS[cat], alpha=0.4, hatch='//')
            ax.set_xticks(list(x))
            ax.set_xticklabels(networks, rotation=45, ha='right', fontsize=8)
            ax.set_ylabel('share of this stage\'s total')
            ax.set_title(f'{stage_title}\n{type_title}', fontsize=9)
    fig.suptitle('Does each category+kind of move find its fair share of modules,\ngiven how much search effort it costs?\n'
                 '(solid bar = % of modules found, hatched twin = % of search effort spent)', fontsize=12)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), fontsize=8, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 0.92))
    fig.tight_layout(rect=(0, 0, 1, 0.84))
    fig.savefig(os.path.join(out_dir, 'risk_contribution_vs_effort.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_risk_efficiency_ratio(df, out_dir):
    """
    Same comparison as risk_contribution_vs_effort.png, collapsed into one
    number per category+kind: (% of modules found) divided by (% of search
    effort spent). Above the dashed line at 1.0 = this combination finds
    more than its "fair share" of modules for the effort it costs
    (efficient). Below the line = costs more effort than the modules it
    finds are worth (inefficient, relatively speaking).
    """
    r = _risk_rows(df)
    if r.empty:
        print('  [skip] risk_efficiency_ratio.png -- no completed risk_analysis runs yet')
        return

    r = r.sort_values('network')
    networks = r['network'].tolist()

    fig, axes = plt.subplots(2, 2, figsize=(max(10, len(networks) * 1.8), 9), sharey=True)
    for row_i, (stage_title, stage_prefix) in enumerate(STAGES):
        ssm_tot = sum(r[f'{stage_prefix}_ssm_{t}_{c}'].fillna(0) for _, t in MOVE_TYPES for c in CATS).replace(0, 1)
        states_tot = sum(r[f'{stage_prefix}_states_{t}_{c}'].fillna(0) for _, t in MOVE_TYPES for c in CATS).replace(0, 1)
        for col_i, (type_title, mtype) in enumerate(MOVE_TYPES):
            ax = axes[row_i, col_i]
            x = range(len(networks))
            width = 0.25
            for ci, cat in enumerate(CATS):
                contrib = r[f'{stage_prefix}_ssm_{mtype}_{cat}'].fillna(0) / ssm_tot
                effort = r[f'{stage_prefix}_states_{mtype}_{cat}'].fillna(0) / states_tot
                ratio = (contrib / effort.replace(0, float('nan'))).astype(float).to_numpy()
                offset = (ci - 1) * width
                ax.bar([xi + offset for xi in x], ratio, width=width, color=CAT_COLORS[cat],
                       label=CAT_LABELS[cat] if (row_i == 0 and col_i == 0) else None)
            ax.axhline(1.0, color='gray', linestyle='--', linewidth=1,
                       label='breaks even' if (row_i == 0 and col_i == 0) else None)
            ax.set_xticks(list(x))
            ax.set_xticklabels(networks, rotation=45, ha='right', fontsize=8)
            ax.set_title(f'{stage_title}\n{type_title}', fontsize=9)
    axes[0, 0].set_ylabel('modules found per unit\nof search effort spent')
    axes[1, 0].set_ylabel('modules found per unit\nof search effort spent')
    fig.suptitle('Which category+kind of move gives the most modules for the search effort it costs?\n'
                 '(1.0 = exactly proportional; higher = more efficient)', fontsize=12)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), fontsize=8, ncol=4, loc='upper center', bbox_to_anchor=(0.5, 0.94))
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(os.path.join(out_dir, 'risk_efficiency_ratio.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_risk_success_rate(df, out_dir):
    """
    Stage 1 (building the smallest modules) only: of every search path that
    used a move of this category+kind, what fraction actually finished as a
    complete module rather than getting stuck (a dead end)? This is a purely
    per-attempt view -- it does NOT account for how many attempts of each
    kind were made (that's what the contribution-vs-effort and efficiency
    plots are for). Stage 2's dead-end counts weren't recorded separately in
    this data collection run, so this view is Stage-1-only. One panel per
    kind of step (synergy/complementarity), matching the other risk plots.
    """
    r = _risk_rows(df)
    if r.empty:
        print('  [skip] risk_success_rate.png -- no completed risk_analysis runs yet')
        return

    r = r.sort_values('network')
    networks = r['network'].tolist()
    fig, axes = plt.subplots(1, 2, figsize=(max(9, len(networks) * 1.6), 5), sharey=True)
    for col_i, (type_title, mtype) in enumerate(MOVE_TYPES):
        ax = axes[col_i]
        x = range(len(networks))
        width = 0.25
        for ci, cat in enumerate(CATS):
            ssm = r[f'riskepm_ssm_{mtype}_{cat}'].fillna(0)
            dead = r[f'riskepm_deadend_{mtype}_{cat}'].fillna(0)
            tot = (ssm + dead).replace(0, float('nan'))
            rate = (ssm / tot).astype(float).to_numpy()
            offset = (ci - 1) * width
            ax.bar([xi + offset for xi in x], rate, width=width, color=CAT_COLORS[cat],
                   label=CAT_LABELS[cat] if col_i == 0 else None)
        ax.set_xticks(list(x))
        ax.set_xticklabels(networks, rotation=45, ha='right', fontsize=8)
        ax.set_title(type_title, fontsize=10)
    axes[0].set_ylabel('% of attempts that succeeded\n(finished as a complete module, not a dead end)')
    fig.suptitle('Stage 1 (building the smallest modules):\nhow often does each category+kind of move actually pay off?')
    fig.legend(*axes[0].get_legend_handles_labels(), fontsize=8, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 0.90))
    fig.tight_layout(rect=(0, 0, 1, 0.82))
    fig.savefig(os.path.join(out_dir, 'risk_success_rate.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    df = load_sweep()
    n_networks = df['network'].nunique()
    print(f'Loaded {len(df)} rows covering {n_networks} networks from {SWEEP_CSV}')

    summary = build_summary_table(df)
    summary_path = os.path.join(OUT_DIR, 'summary_table.csv')
    summary.to_csv(summary_path, index=False)
    print(f'  -> {summary_path}')

    pd.set_option('display.width', 160)
    pd.set_option('display.max_columns', 20)
    cols_to_print = ['network', 'n_ercs', 'syn_density', 'max_scc_fraction',
                      'baseline_wall_s', 'adaptive_wall_s', 'speedup',
                      'dead_end_req_gini', 'espm_total', 'espm_max_order']
    print('\n' + summary[[c for c in cols_to_print if c in summary.columns]].to_string(index=False))

    print('\nGenerating plots...')
    plot_speedup_vs_density(summary, OUT_DIR)
    plot_walltime_vs_size(summary, OUT_DIR)
    plot_deadend_gini_vs_density(summary, OUT_DIR)
    plot_scc_vs_density(summary, OUT_DIR)
    plot_epm_espm_by_order(df, OUT_DIR)
    plot_top_blocking_species(OUT_DIR)
    plot_risk_move_classification(df, OUT_DIR)
    plot_risk_contribution_vs_effort(df, OUT_DIR)
    plot_risk_efficiency_ratio(df, OUT_DIR)
    plot_risk_success_rate(df, OUT_DIR)

    n_risk_done = int(len(_risk_rows(df)))
    n_risk_total = int((df.config == 'risk_analysis').sum())
    print(f'\nRisk-analysis stage coverage: {n_risk_done}/{n_risk_total} networks completed '
          f'(the rest timed out and will retry on the next sweep invocation).')

    espm_all = df[df.config == 'epm_espm_original']
    n_espm_done = int((espm_all['espm_complete'] == True).sum())
    n_espm_partial = int(((espm_all['espm_complete'] == False) & espm_all['espm_total'].notna()).sum())
    n_espm_total = int(len(espm_all))
    print(f'\nESPM stage coverage: {n_espm_done}/{n_espm_total} networks fully complete, '
          f'{n_espm_partial} still resuming (partial data available).')
    print(f'\nAll output in {OUT_DIR}')


if __name__ == '__main__':
    main()
