#!/usr/bin/env python
"""
run_bigg_sweep.py

Sweeps every BiGG network up to MAX_REACTIONS reactions, comparing:
  - "baseline": reach(k) pruning only (validated identical results to the
    unmodified cot_gen.epm.compute_epms — see epm_reach_propagate.py's
    module docstring for the correctness argument)
  - "adaptive": reach(k) pruning + density-gated branch-scoped unit
    propagation (epm_reach_propagate.compute_epms_adaptive)

Results are appended to sweep_results.csv, flushed AND fsynced to disk
after every single run — safe to interrupt (Ctrl+C) at any point;
everything completed so far is already saved when you stop it.

Two kinds of extra diagnostics, beyond wall-time/states comparisons:

1. Progress-on-timeout. A run that hits its time budget still reports
   states_explored/ssm_found/leaves_found so far, plus progress_fraction
   (states_explored / (states_explored + still-pending) -- 1.0 for a
   normal completion) and how deep the still-pending frontier had gotten
   (stack_erc_size_mean/max_at_cutoff) -- distinguishes "barely started"
   from "got deep, ran out of time near the end".

2. Dead-end "economics". Every dead end is characterized by what it still
   needed (req) and what it had already built (prod, erc_set size) when it
   died. req is tracked per-species (cheap -- dead-end req sets are small,
   typically 1-10 species) and summarized per network via a Gini
   coefficient + top-1/top-10% share: a network where a handful of species
   account for most dead-end failures has "bottleneck" resources -- a few
   critical, scarce needs that repeatedly block completion; a flat
   (low-Gini) distribution means failure is spread evenly, no clear
   bottleneck. prod/erc_set are tracked only by size (their popcounts can
   be large; per-species tallying would be expensive and less diagnostic,
   since prod at a dead end is "whatever got built", not a scarce
   resource) -- their mean/max size is the "how much was already invested
   before this branch died" signal. Full per-species req histograms
   (species name + block count) are written to sidecar CSVs in
   deadend_economics/ for deeper offline analysis.

Usage (run from the pyCOT repository root):
    python run_bigg_sweep.py

Requires epm_reach_propagate.py in the same directory as this script.
"""
import sys, os, csv, time, glob, traceback, importlib.util, multiprocessing as mp

sys.path.insert(0, 'src')

from pyCOT.io.functions import read_txt
from pyCOT.analysis.organizations.io_pyCOT import build_rndata
from pyCOT.analysis.organizations.erc import compute_ercs
from pyCOT.analysis.organizations.hierarchy import build_hierarchy
from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
from pyCOT.analysis.organizations.complementarity import compute_complementarities

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    'epm_reach_propagate', os.path.join(_here, 'epm_reach_propagate.py'))
opt = importlib.util.module_from_spec(_spec)
sys.modules['epm_reach_propagate'] = opt
_spec.loader.exec_module(opt)

_espm_spec = importlib.util.spec_from_file_location(
    'espm_resumable', os.path.join(_here, 'espm_resumable.py'))
espm_resumable = importlib.util.module_from_spec(_espm_spec)
sys.modules['espm_resumable'] = espm_resumable
_espm_spec.loader.exec_module(espm_resumable)

_epmrisk_spec = importlib.util.spec_from_file_location(
    'epm_risk_analysis', os.path.join(_here, 'epm_risk_analysis.py'))
epm_risk_analysis = importlib.util.module_from_spec(_epmrisk_spec)
sys.modules['epm_risk_analysis'] = epm_risk_analysis
_epmrisk_spec.loader.exec_module(epm_risk_analysis)

_espmrisk_spec = importlib.util.spec_from_file_location(
    'espm_risk_analysis', os.path.join(_here, 'espm_risk_analysis.py'))
espm_risk_analysis = importlib.util.module_from_spec(_espmrisk_spec)
sys.modules['espm_risk_analysis'] = espm_risk_analysis
_espmrisk_spec.loader.exec_module(espm_risk_analysis)


# =============================================================================
# CONFIGURATION
# =============================================================================
DATA_DIR = 'data/biochemical_databases/biomodels_all_txt'  # BiGG networks (bigg_*.txt) -- see discover_networks()
MIN_REACTIONS = 50     # skip networks smaller than this (0 = no lower bound)
MAX_REACTIONS = 1500
PER_RUN_TIME_BUDGET_S = 600.0   # 5 minutes max, per (network, config) run
ESPM_TIME_BUDGET_S = 600.0      # separate budget for the ESPM stage (see below);
                                 # enforced by hard process kill, not internal checks
ESPM_MAX_ORDER = 20
ESPM_EXTERNAL_KILL_GRACE_S = 600.0  # extra time beyond ESPM_TIME_BUDGET_S before the
                                     # external hard-kill fires -- the internal deadline
                                     # (checked once per BFS round) needs room to finish
                                     # whatever round is in flight and checkpoint cleanly;
                                     # the external kill is only the true last-resort net
                                     # against a single round that itself runs away
RISK_TIME_BUDGET_S = 300.0          # budget for the risk-classification stage (EPM+ESPM
                                     # combined); no checkpointing here (unlike the ESPM
                                     # stage) -- a timeout produces no partial data, and
                                     # (like baseline/adaptive) is simply retried whole on
                                     # the next sweep invocation, not treated as permanent
RISK_EXTERNAL_KILL_GRACE_S = 60.0   # external hard-kill grace beyond RISK_TIME_BUDGET_S
OUTPUT_CSV = os.path.join(_here, 'sweep_results.csv')
SIDECAR_DIR = os.path.join(_here, 'deadend_economics')
ESPM_CKPT_DIR = os.path.join(_here, 'espm_checkpoints')

FIELDNAMES = [
    'timestamp', 'network', 'n_species', 'n_reactions', 'n_ercs',
    'n_fund_synergies', 'n_fund_complementarities',
    'syn_density', 'comp_density', 'syn_to_comp_ratio', 'max_syn_outdegree',
    'max_scc_fraction', 'n_nontrivial_sccs',
    'config', 'effective_cap', 'time_budget_s', 'wall_time_s', 'timed_out',
    'progress_fraction', 'stack_remaining_at_cutoff',
    'stack_erc_size_mean_at_cutoff', 'stack_erc_size_max_at_cutoff',
    'states_explored', 'ssm_found', 'leaves_found',
    'propagation_steps', 'propagation_gated',
    'canonical_pruned', 'reach_pruned_at_construction', 'n_epms',
    'dead_end_req_size_mean', 'dead_end_req_size_max',
    'dead_end_prod_size_mean', 'dead_end_prod_size_max',
    'dead_end_erc_size_mean', 'dead_end_erc_size_max',
    'dead_end_req_gini', 'dead_end_req_top1_share', 'dead_end_req_top10pct_share',
    'n_species_ever_blocking',
    'espm_n_epm', 'espm_total', 'espm_max_order', 'espm_by_order', 'espm_complete',
    # -- risk classification (config='risk_analysis'; see epm_risk_analysis.py /
    #    espm_risk_analysis.py module docstrings for the exact safe/mid/risky
    #    criteria). Raw counts only -- relative/contribution/effort-share views
    #    are derived downstream in analyze_sweep.py from these totals.
    'riskepm_n_epms',
    'riskepm_syn_safe', 'riskepm_syn_mid', 'riskepm_syn_risky',
    'riskepm_comp_safe', 'riskepm_comp_mid', 'riskepm_comp_risky',
    # states/ssm/deadend below: each attributed to whichever move type (syn/comp)
    # is currently responsible for that state's category -- see
    # epm_risk_analysis.py's module docstring. A handful of never-extended seed
    # states have no move type yet, so these per-type columns can sum to
    # slightly less than the untyped total would (there is no untyped total
    # column stored here -- see riskepm_states_syn_* + riskepm_states_comp_*
    # vs. states_explored in the underlying stats dict if that gap matters).
    'riskepm_states_syn_safe', 'riskepm_states_syn_mid', 'riskepm_states_syn_risky',
    'riskepm_states_comp_safe', 'riskepm_states_comp_mid', 'riskepm_states_comp_risky',
    'riskepm_ssm_syn_safe', 'riskepm_ssm_syn_mid', 'riskepm_ssm_syn_risky',
    'riskepm_ssm_comp_safe', 'riskepm_ssm_comp_mid', 'riskepm_ssm_comp_risky',
    'riskepm_deadend_syn_safe', 'riskepm_deadend_syn_mid', 'riskepm_deadend_syn_risky',
    'riskepm_deadend_comp_safe', 'riskepm_deadend_comp_mid', 'riskepm_deadend_comp_risky',
    'riskespm_total', 'riskespm_max_order', 'riskespm_by_order',
    'riskespm_syn_safe', 'riskespm_syn_mid', 'riskespm_syn_risky',
    'riskespm_comp_safe', 'riskespm_comp_mid', 'riskespm_comp_risky',
    'riskespm_hierarchy_unclassified',
    'riskespm_states_syn_safe', 'riskespm_states_syn_mid', 'riskespm_states_syn_risky',
    'riskespm_states_comp_safe', 'riskespm_states_comp_mid', 'riskespm_states_comp_risky',
    'riskespm_ssm_syn_safe', 'riskespm_ssm_syn_mid', 'riskespm_ssm_syn_risky',
    'riskespm_ssm_comp_safe', 'riskespm_ssm_comp_mid', 'riskespm_ssm_comp_risky',
    'error',
]

# config values that can appear in the 'config' column:
#   'baseline'          -- reach(k) pruning only, from epm_reach_propagate
#   'adaptive'           -- reach(k) + density-gated propagation, from epm_reach_propagate
#   'epm_espm_original'  -- original unoptimized cot_gen.epm (both EPM and ESPM
#                           stages), run in an isolated process; only the
#                           espm_* columns (and n_species/n_reactions/n_ercs/
#                           density metrics) are meaningful on this row --
#                           timing/state-count columns are left blank since
#                           the original module has no comparable instrumentation
#   'risk_analysis'      -- safe/mid/risky move classification for both EPM
#                           (Mode-1) and ESPM (Mode-2) discovery, isolated in a
#                           subprocess with a hard kill, no checkpointing (a
#                           timeout is simply retried whole next run, like
#                           baseline/adaptive -- see load_done_configs)
REQUIRED_CONFIGS = ('baseline', 'adaptive', 'epm_espm_original', 'risk_analysis')


def encode_by_order(by_order: dict) -> str:
    """{1: 45, 2: 12} -> '1:45;2:12' -- compact, single-CSV-cell encoding."""
    return ';'.join(f'{k}:{v}' for k, v in sorted(by_order.items()))


def gini(values):
    """Standard Gini coefficient (0 = perfectly even, 1 = maximally concentrated)."""
    if not values:
        return 0.0
    vals = sorted(values)
    n = len(vals)
    total = sum(vals)
    if total == 0:
        return 0.0
    cum = sum((i + 1) * v for i, v in enumerate(vals))
    return (2 * cum) / (n * total) - (n + 1) / n


def top_share(counts_desc, frac):
    """Fraction of total mass held by the top `frac` share of entries (already sorted descending)."""
    if not counts_desc:
        return 0.0
    n = len(counts_desc)
    k = max(1, round(n * frac))
    total = sum(counts_desc)
    return sum(counts_desc[:k]) / total if total else 0.0


def hist_mean_max(hist: dict):
    """(mean, max) of a {size: count} histogram, weighted by count."""
    if not hist:
        return 0.0, 0
    total_n = sum(hist.values())
    mean = sum(size * cnt for size, cnt in hist.items()) / total_n
    return mean, max(hist.keys())


def write_deadend_sidecar(network, config, req_species_counter, species_names):
    """Full per-species dead-end-blocking histogram, species name attached, sorted by count desc."""
    if not req_species_counter:
        return
    os.makedirs(SIDECAR_DIR, exist_ok=True)
    path = os.path.join(SIDECAR_DIR, f'{network}_{config}_deadend_req_species.csv')
    rows = sorted(req_species_counter.items(), key=lambda kv: -kv[1])
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['species_bit', 'species_name', 'block_count'])
        for bit, count in rows:
            name = species_names[bit] if bit < len(species_names) else f'bit{bit}'
            w.writerow([bit, name, count])


def _summarize_espm_result(result: dict, is_complete: bool) -> dict:
    """result: the dict returned by espm_resumable.compute_espm_resumable / a loaded checkpoint's 'result'."""
    by_order = result['espm_by_order']
    counts = {k: len(v) for k, v in by_order.items()}
    return {
        'ok': True,
        'is_complete': is_complete,
        'espm_n_epm': len(result['epm_masks']),
        'espm_by_order': counts,
        'espm_total': sum(counts.values()),
        'espm_max_order': max(counts.keys()) if counts else 0,
    }


def _precompute_cache_path(name):
    return os.path.join(ESPM_CKPT_DIR, f'{name}_precompute.pkl')


def _load_or_compute_precompute(path, name):
    """
    ERCs, hierarchy, fundamental synergies/complementarities, and the full
    order-0 EPM search (via the ORIGINAL, unoptimized cot_gen.epm.compute_epms
    -- the one both the ESPM stage and the risk-analysis stage need) are
    expensive and 100% deterministic from the network file. Both
    _espm_worker and _risk_worker used to recompute all five from scratch on
    EVERY invocation -- including every retry of a network that's still
    resuming its ESPM/risk BFS, and even across the two DIFFERENT stages for
    the same network -- even though only the downstream BFS-by-order state
    (the actual espm_checkpoints/{name}.pkl) was ever being resumed. This
    cache closes that gap: computed once, reused by both stages, invalidated
    only by deleting the file (e.g. after a code change to ERC/hierarchy/
    synergy/complementarity computation -- there is no automatic staleness
    check, same as the ESPM BFS checkpoint itself).
    """
    import pickle
    import sys as _sys
    _sys.path.insert(0, 'src')
    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations.io_pyCOT import build_rndata
    from pyCOT.analysis.organizations.erc import compute_ercs
    from pyCOT.analysis.organizations.hierarchy import build_hierarchy
    from pyCOT.analysis.organizations.synergy import compute_synergies_basis_first
    from pyCOT.analysis.organizations.complementarity import compute_complementarities
    from pyCOT.analysis.organizations.epm import compute_epms as compute_epms_original

    cache_path = _precompute_cache_path(name)
    if os.path.exists(cache_path):
        t0 = time.perf_counter()
        print(f'  [{name}] loading cached ERCs/hierarchy/synergy/complementarity/EPMs '
              f'from {os.path.basename(cache_path)} ...')
        with open(cache_path, 'rb') as f:
            data = pickle.load(f)
        print(f'  [{name}] loaded from cache in {time.perf_counter()-t0:.2f}s '
              f'(n_ercs={len(data["ercs"])}, n_epms={len(data["epm_result"].all_epm_masks)}) '
              f'-- skipped recomputation')
        return data['rn_data'], data['ercs'], data['hier'], data['syn'], data['comp'], data['epm_result']

    print(f'  [{name}] no precompute cache found -- computing ERCs, hierarchy, '
          f'synergy, complementarity, EPMs from scratch...')
    t0 = time.perf_counter()
    rn = read_txt(path, exact_names=True)
    rn_data = build_rndata(rn, network_id=name)
    ercs = compute_ercs(rn_data, verify=False)
    print(f'  [{name}]   ERCs: {len(ercs)}  ({time.perf_counter()-t0:.2f}s elapsed)')
    hier = build_hierarchy(ercs)
    print(f'  [{name}]   hierarchy built  ({time.perf_counter()-t0:.2f}s elapsed)')
    syn = compute_synergies_basis_first(ercs, hier)
    print(f'  [{name}]   fundamental synergies: {len(syn.fundamental)}  ({time.perf_counter()-t0:.2f}s elapsed)')
    comp = compute_complementarities(ercs, hier, syn)
    print(f'  [{name}]   fundamental complementarities: {len(comp.fundamental)}  '
          f'({time.perf_counter()-t0:.2f}s elapsed)')
    epm_result = compute_epms_original(rn_data, ercs, hier, syn_result=syn, comp_result=comp, verbose=False)
    print(f'  [{name}]   order-0 EPMs: {len(epm_result.all_epm_masks)}  '
          f'({time.perf_counter()-t0:.2f}s elapsed, total)')

    os.makedirs(ESPM_CKPT_DIR, exist_ok=True)
    tmp_path = cache_path + '.tmp'
    with open(tmp_path, 'wb') as f:
        pickle.dump({'rn_data': rn_data, 'ercs': ercs, 'hier': hier, 'syn': syn,
                     'comp': comp, 'epm_result': epm_result}, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, cache_path)
    print(f'  [{name}] cached to {os.path.basename(cache_path)} for future resumes/stages')
    return rn_data, ercs, hier, syn, comp, epm_result


def _espm_worker(queue, path, name, max_order, checkpoint_path, progress_csv_path, deadline_ts):
    """
    Runs in a SEPARATE PROCESS (see run_espm_stage). Computes order-level
    EPM/ESPM counts using the original, unmodified, fully-validated
    cot_gen.epm module for the EPM (order-0) stage, and the checkpointed
    espm_resumable.compute_espm_resumable wrapper (a deliberately-unchanged
    copy of cot_gen.epm.compute_espm's per-round logic, validated
    exact-match against it) for the ESPM (order>=1) stage.

    Why a separate process:
    - Even with the internal deadline, a single BFS round can itself run
      long (dense networks) -- the external hard kill below is the true
      safety net against a runaway round, independent of the graceful
      internal checkpoint-and-return path.

    Why resumable: compute_espm_resumable checkpoints to checkpoint_path
    after every completed round (round-boundary pause point), so even a
    hard-killed run's progress up to its last completed round survives on
    disk -- run_espm_stage reads it back directly if the external kill
    fires before the internal deadline gets a chance to return cleanly.
    """
    try:
        rn_data, ercs, hier, syn, comp, epm_result = _load_or_compute_precompute(path, name)
        result, is_complete = espm_resumable.compute_espm_resumable(
            rn_data, ercs, hier, syn, comp, epm_result,
            max_order=max_order, checkpoint_path=checkpoint_path,
            progress_csv_path=progress_csv_path, deadline_ts=deadline_ts, verbose=False)
        queue.put(_summarize_espm_result(result, is_complete))
    except Exception as e:
        queue.put({'ok': False, 'error': str(e)})


def run_espm_stage(path, name, time_budget_s, max_order, checkpoint_path, progress_csv_path):
    """
    Run _espm_worker in an isolated process with a hard wall-clock kill.

    The worker itself stops gracefully (checkpointing) at deadline_ts =
    now + time_budget_s; the external process join uses a longer timeout
    (+ESPM_EXTERNAL_KILL_GRACE_S) so that graceful path gets a real chance
    to fire before the hard kill does. If the hard kill *does* fire first
    (a single round ran long), the checkpoint from the last completed
    round -- if any -- is read back directly so progress isn't lost.
    """
    queue = mp.Queue()
    deadline_ts = time.time() + time_budget_s
    proc = mp.Process(target=_espm_worker,
                       args=(queue, path, name, max_order, checkpoint_path, progress_csv_path, deadline_ts))
    proc.start()
    proc.join(timeout=time_budget_s + ESPM_EXTERNAL_KILL_GRACE_S)
    if proc.is_alive():
        proc.terminate()
        proc.join(5)
        if proc.is_alive():
            proc.kill()
            proc.join(5)
        checkpoint = espm_resumable.load_checkpoint(checkpoint_path)
        if checkpoint is not None:
            result = _summarize_espm_result(checkpoint['result'], checkpoint['is_complete'])
            result['timed_out'] = True   # hard-killed, even though partial data survived
            return result
        return {'ok': False, 'timed_out': True}
    if not queue.empty():
        result = queue.get()
        result['timed_out'] = False
        return result
    return {'ok': False, 'timed_out': False, 'error': 'process exited without a result'}


def discover_networks(data_dir, max_reactions, min_reactions=0):
    """Cheap pre-check (parse only, no ERC computation) to filter by size."""
    paths = sorted(glob.glob(os.path.join(data_dir, '*.txt')))
    found = []
    for path in paths:
        name = os.path.splitext(os.path.basename(path))[0]
        if name.startswith('bigg_'):
            name = name[len('bigg_'):]
        try:
            rn = read_txt(path, exact_names=True)
            rn_data = build_rndata(rn, network_id=name)
        except Exception as e:
            print(f'  [skip] {name}: failed to parse ({e})')
            continue
        if rn_data.n_reactions > max_reactions or rn_data.n_reactions < min_reactions:
            continue
        found.append((path, name, rn_data.n_reactions))
    found.sort(key=lambda t: t[2])   # smallest first: broad, representative
    return found                     # partial results even if stopped early


def write_row(writer, f, row):
    full_row = {k: row.get(k, '') for k in FIELDNAMES}
    writer.writerow(full_row)
    f.flush()
    os.fsync(f.fileno())


def stats_to_row_fields(stats, result):
    """Extract every diagnostic field (progress + dead-end economics) from a run's stats dict."""
    req_size_mean, req_size_max = hist_mean_max(stats.get('dead_end_req_size_hist', {}))
    prod_size_mean, prod_size_max = hist_mean_max(stats.get('dead_end_prod_size_hist', {}))
    erc_size_mean, erc_size_max = hist_mean_max(stats.get('dead_end_erc_size_hist', {}))

    req_species = stats.get('dead_end_req_species', {})
    counts = list(req_species.values())
    counts_desc = sorted(counts, reverse=True)

    return {
        'timed_out': stats.get('time_budget_exceeded', False),
        'progress_fraction': f'{stats.get("progress_fraction", 1.0):.4f}',
        'stack_remaining_at_cutoff': stats.get('stack_remaining_at_cutoff', 0),
        'stack_erc_size_mean_at_cutoff': (f'{stats["stack_erc_size_mean_at_cutoff"]:.1f}'
                                           if 'stack_erc_size_mean_at_cutoff' in stats else ''),
        'stack_erc_size_max_at_cutoff': stats.get('stack_erc_size_max_at_cutoff', ''),
        'states_explored': stats.get('states_explored', 0),
        'ssm_found': stats.get('ssm_found', 0),
        'leaves_found': stats.get('leaves_found', 0),
        'propagation_steps': stats.get('propagation_steps', 0),
        'propagation_gated': stats.get('propagation_gated', 0),
        'canonical_pruned': stats.get('canonical_pruned', 0),
        'reach_pruned_at_construction': stats.get('reach_pruned_at_construction', 0),
        'n_epms': len(result.all_epm_masks),
        'dead_end_req_size_mean': f'{req_size_mean:.2f}', 'dead_end_req_size_max': req_size_max,
        'dead_end_prod_size_mean': f'{prod_size_mean:.2f}', 'dead_end_prod_size_max': prod_size_max,
        'dead_end_erc_size_mean': f'{erc_size_mean:.2f}', 'dead_end_erc_size_max': erc_size_max,
        'dead_end_req_gini': f'{gini(counts):.3f}',
        'dead_end_req_top1_share': f'{top_share(counts_desc, 1.0/max(len(counts_desc),1)):.3f}',
        'dead_end_req_top10pct_share': f'{top_share(counts_desc, 0.10):.3f}',
        'n_species_ever_blocking': len(req_species),
    }, req_species


def run_baseline_stage(base_row, rn_data, ercs, hier, syn, comp, reach_table, name):
    row = dict(base_row)
    row['config'] = 'baseline'
    row['time_budget_s'] = PER_RUN_TIME_BUDGET_S
    try:
        t0 = time.perf_counter()
        result = opt.compute_epms(rn_data, ercs, hier, syn_result=syn, comp_result=comp,
                                   time_budget_s=PER_RUN_TIME_BUDGET_S,
                                   reach_table=reach_table, propagation_cap=0)
        wall = time.perf_counter() - t0
        stats = result.stats
        row['effective_cap'] = 0
        row['wall_time_s'] = f'{wall:.2f}'
        fields, req_species = stats_to_row_fields(stats, result)
        row.update(fields)
        write_deadend_sidecar(name, 'baseline', req_species, rn_data.species_names)
        print(f'  baseline  wall={wall:9.2f}s  {"TIMEOUT" if stats.get("time_budget_exceeded") else "done   "}'
              f'  progress={fields["progress_fraction"]}  states={stats.get("states_explored",0):8d}  '
              f'ssm={stats.get("ssm_found",0)}  n_epms={len(result.all_epm_masks)}')
    except Exception as e:
        row['error'] = str(e)
        print(f'  baseline !!! ERROR: {e}')
        traceback.print_exc()
    return row


def run_adaptive_stage(base_row, rn_data, ercs, hier, syn, comp, reach_table, density, name):
    row = dict(base_row)
    row['config'] = 'adaptive'
    row['time_budget_s'] = PER_RUN_TIME_BUDGET_S
    try:
        t0 = time.perf_counter()
        effective_cap = 0 if density > opt.DENSITY_GUARD_THRESHOLD else opt.DEFAULT_PROPAGATION_CAP
        result = opt.compute_epms(rn_data, ercs, hier, syn_result=syn, comp_result=comp,
                                   time_budget_s=PER_RUN_TIME_BUDGET_S,
                                   reach_table=reach_table, propagation_cap=effective_cap)
        wall = time.perf_counter() - t0
        stats = result.stats
        row['effective_cap'] = effective_cap
        row['wall_time_s'] = f'{wall:.2f}'
        fields, req_species = stats_to_row_fields(stats, result)
        row.update(fields)
        write_deadend_sidecar(name, 'adaptive', req_species, rn_data.species_names)
        print(f'  adaptive  wall={wall:9.2f}s  {"TIMEOUT" if stats.get("time_budget_exceeded") else "done   "}'
              f'  progress={fields["progress_fraction"]}  states={stats.get("states_explored",0):8d}  '
              f'ssm={stats.get("ssm_found",0)}  n_epms={len(result.all_epm_masks)}  (cap={effective_cap})')
    except Exception as e:
        row['error'] = str(e)
        print(f'  adaptive !!! ERROR: {e}')
        traceback.print_exc()
    return row


def run_espm_original_stage(base_row, path, name):
    """
    Order-level EPM/ESPM counts via the original, unoptimized cot_gen.epm
    module for EPMs, and the checkpointed espm_resumable wrapper for ESPMs,
    isolated in a subprocess with a hard kill (see run_espm_stage /
    _espm_worker docstrings for why).

    Unlike the old one-shot version, this now checkpoints after every
    completed BFS round (see espm_resumable.py) to
    espm_checkpoints/{name}.pkl: a timed-out or hard-killed run still
    reports whatever ESPMs were found up through its last completed round,
    and a later invocation (e.g. the next sweep run) picks up exactly where
    it left off instead of restarting from order 1. row['espm_complete']
    distinguishes a genuinely-finished BFS (True) from a paused one that
    will resume next time (False) -- load_done_configs uses this to decide
    whether to re-invoke this stage.
    """
    row = dict(base_row)
    row['config'] = 'epm_espm_original'
    row['time_budget_s'] = ESPM_TIME_BUDGET_S
    os.makedirs(ESPM_CKPT_DIR, exist_ok=True)
    checkpoint_path = os.path.join(ESPM_CKPT_DIR, f'{name}.pkl')
    progress_csv_path = os.path.join(ESPM_CKPT_DIR, f'{name}_progress.csv')
    t0 = time.perf_counter()
    result = run_espm_stage(path, name, ESPM_TIME_BUDGET_S, ESPM_MAX_ORDER, checkpoint_path, progress_csv_path)
    wall = time.perf_counter() - t0
    row['wall_time_s'] = f'{wall:.2f}'
    if not result.get('ok'):
        row['timed_out'] = result.get('timed_out', False)
        if row['timed_out']:
            print(f'  espm      wall={wall:9.2f}s  TIMEOUT  (hard-killed, no checkpoint yet -- no partial data)')
        else:
            row['error'] = result.get('error', 'unknown error')
            print(f'  espm      wall={wall:9.2f}s  !!! ERROR: {row["error"]}')
    else:
        row['timed_out'] = result.get('timed_out', False)
        row['espm_complete'] = result['is_complete']
        row['espm_n_epm'] = result['espm_n_epm']
        row['espm_total'] = result['espm_total']
        row['espm_max_order'] = result['espm_max_order']
        row['espm_by_order'] = encode_by_order(result['espm_by_order'])
        status = 'done' if result['is_complete'] else ('killed, checkpointed' if row['timed_out'] else 'paused, checkpointed')
        print(f'  espm      wall={wall:9.2f}s  {status:22s} n_epm={result["espm_n_epm"]}  '
              f'n_espm_total={result["espm_total"]}  max_order={result["espm_max_order"]}')
    return row


def _risk_worker(queue, path, name, espm_max_order):
    """
    Runs in a SEPARATE PROCESS (see run_risk_stage). Computes the safe/mid/
    risky move classification for both EPM (Mode-1) and ESPM (Mode-2)
    discovery -- see epm_risk_analysis.py / espm_risk_analysis.py module
    docstrings for the exact criteria and the correctness argument (both
    validated exact-match against epm_reach_propagate / espm_resumable
    before being trusted). Isolated the same way as _espm_worker: neither
    the EPM nor ESPM risk engines carry their own internal time budget for
    the FULL traversal (epm_risk_analysis does accept one and stops
    cleanly; espm_risk_analysis's outer round loop also respects one) --
    but a single very slow round/branch could still run long, so the
    external hard-kill in run_risk_stage is the real safety net, same
    reasoning as the ESPM stage. No checkpointing here (unlike ESPM): if
    killed, nothing is reported for this attempt -- retried whole next run.
    """
    try:
        rn_data, ercs, hier, syn, comp, epm_result_for_espm = _load_or_compute_precompute(path, name)

        epm_r = epm_risk_analysis.compute_epms_risk(
            rn_data, ercs, hier, syn_result=syn, comp_result=comp)
        s = epm_r.stats

        # ESPM risk needs the ORIGINAL compute_epms's EPMResult (._graph etc.)
        # as its Mode-2 starting point -- same reason as the ESPM stage. Reused
        # from the shared precompute cache above rather than recomputed here.
        espm_r = espm_risk_analysis.compute_espm_risk(
            rn_data, ercs, hier, syn, comp, epm_result_for_espm, max_order=espm_max_order)
        es = espm_r.stats
        espm_counts = {k: len(v) for k, v in espm_r.espm_by_order.items()}

        queue.put({
            'ok': True,
            'riskepm_n_epms': len(epm_r.all_epm_masks),
            'riskepm_syn_safe': s.get('syn_moves_safe', 0),
            'riskepm_syn_mid': s.get('syn_moves_mid', 0),
            'riskepm_syn_risky': s.get('syn_moves_risky', 0),
            'riskepm_comp_safe': s.get('comp_moves_safe', 0),
            'riskepm_comp_mid': s.get('comp_moves_mid', 0),
            'riskepm_comp_risky': s.get('comp_moves_risky', 0),
            'riskepm_states_syn_safe': s.get('states_explored_syn_safe', 0),
            'riskepm_states_syn_mid': s.get('states_explored_syn_mid', 0),
            'riskepm_states_syn_risky': s.get('states_explored_syn_risky', 0),
            'riskepm_states_comp_safe': s.get('states_explored_comp_safe', 0),
            'riskepm_states_comp_mid': s.get('states_explored_comp_mid', 0),
            'riskepm_states_comp_risky': s.get('states_explored_comp_risky', 0),
            'riskepm_ssm_syn_safe': s.get('ssm_via_syn_safe', 0),
            'riskepm_ssm_syn_mid': s.get('ssm_via_syn_mid', 0),
            'riskepm_ssm_syn_risky': s.get('ssm_via_syn_risky', 0),
            'riskepm_ssm_comp_safe': s.get('ssm_via_comp_safe', 0),
            'riskepm_ssm_comp_mid': s.get('ssm_via_comp_mid', 0),
            'riskepm_ssm_comp_risky': s.get('ssm_via_comp_risky', 0),
            'riskepm_deadend_syn_safe': s.get('dead_end_via_syn_safe', 0),
            'riskepm_deadend_syn_mid': s.get('dead_end_via_syn_mid', 0),
            'riskepm_deadend_syn_risky': s.get('dead_end_via_syn_risky', 0),
            'riskepm_deadend_comp_safe': s.get('dead_end_via_comp_safe', 0),
            'riskepm_deadend_comp_mid': s.get('dead_end_via_comp_mid', 0),
            'riskepm_deadend_comp_risky': s.get('dead_end_via_comp_risky', 0),
            'riskespm_total': sum(espm_counts.values()),
            'riskespm_max_order': max(espm_counts.keys()) if espm_counts else 0,
            'riskespm_by_order': espm_counts,
            'riskespm_syn_safe': es.get('mode2_syn_moves_safe', 0),
            'riskespm_syn_mid': es.get('mode2_syn_moves_mid', 0),
            'riskespm_syn_risky': es.get('mode2_syn_moves_risky', 0),
            'riskespm_comp_safe': es.get('mode2_comp_moves_safe', 0),
            'riskespm_comp_mid': es.get('mode2_comp_moves_mid', 0),
            'riskespm_comp_risky': es.get('mode2_comp_moves_risky', 0),
            'riskespm_hierarchy_unclassified': es.get('mode2_hierarchy_moves', 0),
            'riskespm_states_syn_safe': es.get('states_explored_syn_safe', 0),
            'riskespm_states_syn_mid': es.get('states_explored_syn_mid', 0),
            'riskespm_states_syn_risky': es.get('states_explored_syn_risky', 0),
            'riskespm_states_comp_safe': es.get('states_explored_comp_safe', 0),
            'riskespm_states_comp_mid': es.get('states_explored_comp_mid', 0),
            'riskespm_states_comp_risky': es.get('states_explored_comp_risky', 0),
            'riskespm_ssm_syn_safe': es.get('ssm_via_syn_safe', 0),
            'riskespm_ssm_syn_mid': es.get('ssm_via_syn_mid', 0),
            'riskespm_ssm_syn_risky': es.get('ssm_via_syn_risky', 0),
            'riskespm_ssm_comp_safe': es.get('ssm_via_comp_safe', 0),
            'riskespm_ssm_comp_mid': es.get('ssm_via_comp_mid', 0),
            'riskespm_ssm_comp_risky': es.get('ssm_via_comp_risky', 0),
        })
    except Exception as e:
        queue.put({'ok': False, 'error': str(e)})


def run_risk_stage(path, name, time_budget_s, espm_max_order):
    """Run _risk_worker in an isolated process with a hard wall-clock kill (no checkpointing)."""
    queue = mp.Queue()
    proc = mp.Process(target=_risk_worker, args=(queue, path, name, espm_max_order))
    proc.start()
    proc.join(timeout=time_budget_s + RISK_EXTERNAL_KILL_GRACE_S)
    if proc.is_alive():
        proc.terminate()
        proc.join(5)
        if proc.is_alive():
            proc.kill()
            proc.join(5)
        return {'ok': False, 'timed_out': True}
    if not queue.empty():
        result = queue.get()
        result['timed_out'] = False
        return result
    return {'ok': False, 'timed_out': False, 'error': 'process exited without a result'}


def run_risk_analysis_stage(base_row, path, name):
    """
    Safe/mid/risky move classification for EPM + ESPM discovery -- see
    epm_risk_analysis.py / espm_risk_analysis.py. Raw counts only; relative
    (contribution share, exploration-effort share, efficiency ratio) views
    are computed downstream in analyze_sweep.py.
    """
    row = dict(base_row)
    row['config'] = 'risk_analysis'
    row['time_budget_s'] = RISK_TIME_BUDGET_S
    t0 = time.perf_counter()
    result = run_risk_stage(path, name, RISK_TIME_BUDGET_S, ESPM_MAX_ORDER)
    wall = time.perf_counter() - t0
    row['wall_time_s'] = f'{wall:.2f}'
    if not result.get('ok'):
        row['timed_out'] = result.get('timed_out', False)
        if row['timed_out']:
            print(f'  risk      wall={wall:9.2f}s  TIMEOUT  (hard-killed, no data -- will retry next run)')
        else:
            row['error'] = result.get('error', 'unknown error')
            print(f'  risk      wall={wall:9.2f}s  !!! ERROR: {row["error"]}')
    else:
        row['timed_out'] = False
        for k, v in result.items():
            if k in ('ok', 'timed_out'):
                continue
            row[k] = encode_by_order(v) if k == 'riskespm_by_order' else v
        print(f'  risk      wall={wall:9.2f}s  done      '
              f'n_epms={result["riskepm_n_epms"]}  '
              f'epm syn(safe/mid/risky)={result["riskepm_syn_safe"]}/{result["riskepm_syn_mid"]}/{result["riskepm_syn_risky"]}  '
              f'espm_total={result["riskespm_total"]}  '
              f'espm syn(safe/mid/risky)={result["riskespm_syn_safe"]}/{result["riskespm_syn_mid"]}/{result["riskespm_syn_risky"]}')
    return row


def run_network(writer, f, path, name, done_configs=frozenset()):
    """
    done_configs: set of config strings already computed for this network
    (from a prior run of this script) -- each one already present is
    skipped entirely, no recomputation.
    """
    needed = [c for c in REQUIRED_CONFIGS if c not in done_configs]
    if not needed:
        print('  [skip] all configs already computed')
        return

    base_row = {'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'), 'network': name}
    setup_ok = True
    try:
        rn = read_txt(path, exact_names=True)
        rn_data = build_rndata(rn, network_id=name)
        ercs = compute_ercs(rn_data, verify=False)
        hier = build_hierarchy(ercs)
        syn = compute_synergies_basis_first(ercs, hier)
        comp = compute_complementarities(ercs, hier, syn)
        n_ercs = len(ercs)
        metrics = opt.compute_structural_metrics(n_ercs, syn, comp)
        density = metrics['syn_density']   # kept as local name for the DENSITY_GUARD check below
        base_row.update({
            'n_species': rn_data.n_species, 'n_reactions': rn_data.n_reactions,
            'n_ercs': n_ercs, 'n_fund_synergies': len(syn.fundamental),
            'n_fund_complementarities': len(comp.fundamental),
            'syn_density': f'{metrics["syn_density"]:.2f}',
            'comp_density': f'{metrics["comp_density"]:.2f}',
            'syn_to_comp_ratio': f'{metrics["syn_to_comp_ratio"]:.2f}',
            'max_syn_outdegree': metrics['max_syn_outdegree'],
            'max_scc_fraction': f'{metrics["max_scc_fraction"]:.3f}',
            'n_nontrivial_sccs': metrics['n_nontrivial_sccs'],
        })
        print(f'  setup: n_reactions={rn_data.n_reactions}  n_ercs={n_ercs}  '
              f'syn_density={metrics["syn_density"]:.2f}  comp_density={metrics["comp_density"]:.2f}  '
              f'max_scc_fraction={metrics["max_scc_fraction"]:.3f}')
    except Exception as e:
        setup_ok = False
        base_row['error'] = f'setup failed: {e}'
        print(f'  !!! setup failed: {e}')
        traceback.print_exc()

    if not setup_ok:
        for cfg in needed:
            write_row(writer, f, {**base_row, 'config': cfg, 'time_budget_s': PER_RUN_TIME_BUDGET_S})
        return

    g = opt.FundamentalGraph(ercs, hier, syn, comp)
    reach_table = opt.build_reach_table(g, n_ercs)

    if 'baseline' in needed:
        write_row(writer, f, run_baseline_stage(base_row, rn_data, ercs, hier, syn, comp, reach_table, name))
    else:
        print('  [skip] baseline already computed')

    if 'adaptive' in needed:
        write_row(writer, f, run_adaptive_stage(base_row, rn_data, ercs, hier, syn, comp, reach_table, density, name))
    else:
        print('  [skip] adaptive already computed')

    if 'epm_espm_original' in needed:
        write_row(writer, f, run_espm_original_stage(base_row, path, name))
    else:
        print('  [skip] epm_espm_original already computed')

    if 'risk_analysis' in needed:
        write_row(writer, f, run_risk_analysis_stage(base_row, path, name))
    else:
        print('  [skip] risk_analysis already computed')


def load_done_configs(output_csv, expected_fieldnames):
    """
    Resume support: read any existing results file and return
    {network: {config, ...}} for rows already computed, so run_network can
    skip them. If the file's header doesn't match the CURRENT schema (e.g.
    new columns were added -- as just happened with the ESPM fields), the
    old data can't be safely appended to: it's archived under a timestamped
    name instead of silently mixed with the new schema, and a fresh sweep
    starts (this run will recompute everything, exactly once, since the
    archived file won't be seen as "existing" for the new OUTPUT_CSV path).

    A row is only "done" (skippable) if it actually finished:
    - epm_espm_original: only when espm_complete == 'True'. A paused/
      hard-killed row is left OUT of the done set so run_network invokes
      this stage again -- which transparently resumes from its checkpoint
      (see run_espm_original_stage) rather than starting over.
    - baseline/adaptive/risk_analysis: only when timed_out != 'True'. These
      have no checkpointing, so a timed-out run is simply retried from
      scratch; still better than permanently treating a partial result as
      final.
    """
    if not os.path.exists(output_csv):
        return {}

    with open(output_csv, newline='') as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            return {}   # empty file, nothing to preserve

    if header != expected_fieldnames:
        backup = output_csv.replace('.csv', f'.schema_{time.strftime("%Y%m%d_%H%M%S")}.bak.csv')
        os.rename(output_csv, backup)
        print(f'NOTE: {output_csv} used an older column schema (e.g. before ESPM support was added).')
        print(f'      It has been preserved as: {backup}')
        print(f'      Starting a fresh sweep on the current schema -- all networks will be recomputed once.\n')
        return {}

    done: dict[str, set] = {}
    with open(output_csv, newline='') as f:
        for row in csv.DictReader(f):
            if row.get('error'):
                continue   # a failed run isn't "done" -- worth retrying
            cfg = row['config']
            if cfg == 'epm_espm_original':
                if row.get('espm_complete') != 'True':
                    continue   # paused/killed -- re-invoke to resume from checkpoint
            elif row.get('timed_out') == 'True':
                continue       # baseline/adaptive hit their time budget -- worth retrying
            done.setdefault(row['network'], set()).add(cfg)
    return done


def main():
    bound_desc = f'{MIN_REACTIONS}-{MAX_REACTIONS}' if MIN_REACTIONS > 0 else f'up to {MAX_REACTIONS}'
    print(f'Discovering BiGG networks with {bound_desc} reactions in {DATA_DIR} ...')
    networks = discover_networks(DATA_DIR, MAX_REACTIONS, MIN_REACTIONS)
    print(f'Found {len(networks)} networks in range (smallest reactions first).')

    done_by_network = load_done_configs(OUTPUT_CSV, FIELDNAMES)
    n_fully_done = sum(1 for cfgs in done_by_network.values() if set(REQUIRED_CONFIGS) <= cfgs)
    if done_by_network:
        print(f'Resuming: {n_fully_done} network(s) already fully computed, will be skipped.')
    print(f'Results -> {OUTPUT_CSV}  (flushed after every run -- safe to Ctrl+C at any time)\n')

    write_header = not os.path.exists(OUTPUT_CSV)
    with open(OUTPUT_CSV, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        if write_header:
            writer.writeheader()
            f.flush()
            os.fsync(f.fileno())

        try:
            for i, (path, name, n_reactions) in enumerate(networks, 1):
                print(f'\n=== [{i}/{len(networks)}] {name}  (n_reactions={n_reactions}) ===')
                run_network(writer, f, path, name, done_configs=done_by_network.get(name, frozenset()))
        except KeyboardInterrupt:
            print(f'\n\nInterrupted by user -- all completed runs already saved to {OUTPUT_CSV}.')
            print('Exiting cleanly.')
            sys.exit(0)

    print(f'\nSweep complete. Results in {OUTPUT_CSV}')


if __name__ == '__main__':
    main()
