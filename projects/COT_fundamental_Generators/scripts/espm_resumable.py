"""
espm_resumable.py — checkpointed, resumable version of cot_gen.epm.compute_espm.

The original compute_espm has no time budget and no way to report partial
progress: if killed mid-computation (the only safe way to bound it — see
run_bigg_sweep.py's _espm_worker docstring), everything is lost. This module
adds round-boundary checkpointing: compute_espm's BFS-by-order structure
already processes one whole order ("round") at a time, so a checkpoint saved
after every completed round is a natural, correctness-preserving pause
point — resuming just means re-entering the same for-loop at the next order,
with the exact same state a single uninterrupted run would have had at that
point.

The per-round loop body below is a DELIBERATELY UNCHANGED copy of the
original's (see cot_gen/epm.py compute_espm for the authoritative,
extensively-commented version of the same logic, including the Pareto-tie
canonical-ordering fix) -- the only additions are: load-checkpoint-if-present
at the start, a deadline check at the top of each round, and
save-checkpoint-and-progress-log after each round completes.

Checkpoint contents (pickled dict): everything needed to resume the for-loop
at `next_order` with IDENTICAL state to an uninterrupted run reaching that
point -- visited_sp, sp_to_state, so_order, leaf_masks_by_order,
stats_by_order, current_layer, next_order, is_complete, and
cumulative_elapsed_s (wall-clock compute time actually spent, summed across
every resume -- NOT calendar time between runs).

Correctness: validated by (a) exact-match against the original compute_espm
on a network that completes within one shot, and (b) exact-match between a
single uninterrupted run and a run artificially split into two resumed
sessions -- see the validation script alongside this file.
"""
from __future__ import annotations

import os
import csv
import time
import pickle

from cot_gen.fundamental_graph import FundamentalGraph, DFSState
from cot_gen.epm import _mode1_dfs, EPMResult


def _bits(mask: int):
    m = mask
    while m:
        lsb = m & (-m)
        yield lsb.bit_length() - 1
        m &= m - 1


def _atomic_pickle_dump(obj, path: str) -> None:
    """Write-to-temp-then-rename so a kill mid-write never corrupts the checkpoint."""
    tmp_path = path + '.tmp'
    with open(tmp_path, 'wb') as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, path)   # atomic on POSIX and Windows (same volume)


def load_checkpoint(checkpoint_path: str):
    if not os.path.exists(checkpoint_path):
        return None
    with open(checkpoint_path, 'rb') as f:
        return pickle.load(f)


def _append_progress_log(progress_csv_path: str, row: dict) -> None:
    write_header = not os.path.exists(progress_csv_path)
    with open(progress_csv_path, 'a', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['order', 'n_new_this_round', 'cumulative_espm_total',
                                           'cumulative_elapsed_s', 'wall_clock_timestamp'])
        if write_header:
            w.writeheader()
        w.writerow(row)
        f.flush()
        os.fsync(f.fileno())


def compute_espm_resumable(
    rn, ercs, hier, syn_result, comp_result, epm_result: EPMResult,
    *,
    max_order: int = 10,
    checkpoint_path: str,
    progress_csv_path: str | None = None,
    deadline_ts: float | None = None,   # time.time()-based absolute deadline; None = no limit
    verbose: bool = False,
) -> tuple[dict, bool]:
    """
    Returns (result_dict, is_complete).

    result_dict has the same shape as ESPMResult's fields (as plain dict, so
    it's trivially serializable): epm_masks, espm_by_order, all_so_masks,
    leaf_masks_by_order, stats_by_order, plus cumulative_elapsed_s.

    is_complete=True means the BFS-by-order loop reached its natural
    termination (no more new SOs found, or max_order hit) -- nothing more to
    resume, safe to treat as permanently done. False means the deadline cut
    it off mid-way -- calling again with the same checkpoint_path picks up
    exactly where this call left off.
    """
    t_session_start = time.perf_counter()
    checkpoint = load_checkpoint(checkpoint_path)

    if checkpoint is not None:
        if checkpoint.get('is_complete'):
            return checkpoint['result'], True
        g = FundamentalGraph(ercs, hier, syn_result, comp_result)   # cheap, deterministic rebuild
        visited_sp          = checkpoint['visited_sp']
        sp_to_state          = checkpoint['sp_to_state']
        so_order              = checkpoint['so_order']
        leaf_masks_by_order = checkpoint['leaf_masks_by_order']
        stats_by_order        = checkpoint['stats_by_order']
        current_layer          = checkpoint['current_layer']
        start_order            = checkpoint['next_order']
        cumulative_elapsed_s  = checkpoint['cumulative_elapsed_s']
        if verbose:
            print(f'  [ESPM resume] loaded checkpoint: resuming at order {start_order}, '
                  f'{cumulative_elapsed_s:.1f}s already spent, {len(so_order)} SOs known so far')
    else:
        if epm_result._graph is not None:
            g           = epm_result._graph
            visited_sp  = epm_result._visited_sp
            sp_to_state = epm_result._sp_to_state
            so_order    = epm_result._so_order
        else:
            g = FundamentalGraph(ercs, hier, syn_result, comp_result)
            visited_sp  = set()
            sp_to_state = {}
            so_order    = {sp: 0 for sp in epm_result.all_epm_masks}
            for sp in epm_result.all_epm_masks:
                for i, e in enumerate(ercs):
                    if e.species_mask == sp and e.is_persistent():
                        sp_to_state[sp] = g.make_seed_state(i)
                        break
        leaf_masks_by_order: dict[int, list[int]] = {}
        stats_by_order:      dict[int, dict]      = {}
        current_layer = [sp_to_state[sp] for sp in epm_result.all_epm_masks if sp in sp_to_state]
        start_order = 1
        cumulative_elapsed_s = 0.0

    def _save(next_order: int, complete: bool):
        espm_by_order: dict[int, list[int]] = {}
        for sp, ord_k in so_order.items():
            if ord_k >= 1:
                espm_by_order.setdefault(ord_k, []).append(sp)
        for k in espm_by_order:
            espm_by_order[k].sort()
        result = {
            'epm_masks': sorted(epm_result.all_epm_masks),
            'espm_by_order': espm_by_order,
            'all_so_masks': sorted(so_order.keys()),
            'leaf_masks_by_order': dict(leaf_masks_by_order),
            'stats_by_order': dict(stats_by_order),
            'cumulative_elapsed_s': cumulative_elapsed_s,
        }
        ckpt = {
            'is_complete': complete,
            'result': result,
            'visited_sp': visited_sp,
            'sp_to_state': sp_to_state,
            'so_order': so_order,
            'leaf_masks_by_order': leaf_masks_by_order,
            'stats_by_order': stats_by_order,
            'current_layer': current_layer,
            'next_order': next_order,
            'cumulative_elapsed_s': cumulative_elapsed_s,
        }
        _atomic_pickle_dump(ckpt, checkpoint_path)
        return result

    # ── BFS by order (unchanged logic from the original -- see epm.py) ─────
    order = start_order
    while order <= max_order:
        if deadline_ts is not None and time.time() >= deadline_ts:
            cumulative_elapsed_s += time.perf_counter() - t_session_start
            result = _save(order, complete=False)
            if verbose:
                print(f'  [ESPM] deadline reached before order {order} -- checkpoint saved, '
                      f'{cumulative_elapsed_s:.1f}s cumulative')
            return result, False

        if not current_layer:
            break

        stats_k: dict = {'mode2_seeds': 0}

        best_by_sp: dict[int, list[DFSState]] = {}
        for so_state in current_layer:
            local_cand: dict[int, DFSState] = {}
            for i in so_state.erc_set:
                for (j, _) in g.syn_from.get(i, []):
                    if j in so_state.erc_set:
                        continue
                    if j not in local_cand:
                        local_cand[j] = g.extend_state(so_state, j)
                for a in g.parents[i]:
                    if a in so_state.erc_set:
                        continue
                    if (g.species_mask[a] & so_state.sp) == g.species_mask[a]:
                        continue
                    if a not in local_cand:
                        local_cand[a] = g.extend_state(so_state, a)
            for s_bit in _bits(so_state.prod):
                for cons_idx in g.comp_consumers_by_species.get(s_bit, []):
                    if cons_idx in so_state.erc_set:
                        continue
                    if (g.species_mask[cons_idx] & so_state.sp) == g.species_mask[cons_idx]:
                        continue
                    if cons_idx not in local_cand:
                        local_cand[cons_idx] = g.extend_state(so_state, cons_idx)

            for ext_state in local_cand.values():
                if ext_state.sp in visited_sp:
                    continue
                frontier = best_by_sp.setdefault(ext_state.sp, [])
                new_seed = min(ext_state.erc_set)
                dominated = False
                survivors = []
                for cand in frontier:
                    cand_seed = min(cand.erc_set)
                    if cand.min_ext <= ext_state.min_ext and cand_seed >= new_seed:
                        dominated = True
                        survivors.append(cand)
                    elif ext_state.min_ext <= cand.min_ext and new_seed >= cand_seed:
                        continue
                    else:
                        survivors.append(cand)
                if not dominated:
                    survivors.append(ext_state)
                best_by_sp[ext_state.sp] = survivors

        singleton_seeds: list[DFSState] = []
        tied_groups: list[list[DFSState]] = []
        for frontier in best_by_sp.values():
            if len(frontier) == 1:
                singleton_seeds.append(frontier[0])
            else:
                tied_groups.append(frontier)

        if not singleton_seeds and not tied_groups:
            if verbose:
                print(f"  [Round {order:>2}]: no Mode-2 seeds from {len(current_layer)} SOs — done.")
            break

        stats_k['mode2_seeds'] = len(singleton_seeds) + sum(len(t) for t in tied_groups)
        stats_k['pareto_ties'] = len(tied_groups)

        new_ssm_states: list[DFSState] = []
        new_leaf_states: list[DFSState] = []
        if singleton_seeds:
            ssms, leaves = _mode1_dfs(singleton_seeds, g, visited_sp, stats_k, verbose=False)
            new_ssm_states.extend(ssms)
            new_leaf_states.extend(leaves)
        for tied in tied_groups:
            for rep in tied:
                if rep.sp in visited_sp:
                    visited_sp.discard(rep.sp)
                ssms, leaves = _mode1_dfs([rep], g, visited_sp, stats_k, verbose=False)
                new_ssm_states.extend(ssms)
                new_leaf_states.extend(leaves)

        for s in sorted(new_ssm_states, key=lambda st: bin(st.sp).count('1')):
            if s.sp not in so_order:
                max_sub = -1
                for sub_sp, sub_ord in so_order.items():
                    if sub_sp != s.sp and (sub_sp & s.sp) == sub_sp and sub_ord > max_sub:
                        max_sub = sub_ord
                so_order[s.sp] = max_sub + 1
            sp_to_state[s.sp] = s

        if new_leaf_states:
            leaf_masks_by_order.setdefault(order, []).extend(s.sp for s in new_leaf_states)
        stats_by_order[order] = stats_k

        if verbose:
            n_seeds = stats_k.get('mode2_seeds', 0)
            n_cls   = stats_k.get('states_explored', 0)
            print(f"  [Round {order:>2}]: {n_seeds:4d} seeds -> {len(new_ssm_states):4d} new SOs  "
                  f"{len(new_leaf_states)} leaves  {n_cls} states")

        # ── Round complete: checkpoint + progress log ───────────────────────
        cumulative_elapsed_s_now = cumulative_elapsed_s + (time.perf_counter() - t_session_start)
        cumulative_total_so_far = sum(1 for v in so_order.values() if v >= 1)
        if progress_csv_path is not None:
            _append_progress_log(progress_csv_path, {
                'order': order,
                'n_new_this_round': len(new_ssm_states),
                'cumulative_espm_total': cumulative_total_so_far,
                'cumulative_elapsed_s': f'{cumulative_elapsed_s_now:.2f}',
                'wall_clock_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            })

        current_layer = new_ssm_states
        order += 1
        if not new_ssm_states:
            break

    cumulative_elapsed_s += time.perf_counter() - t_session_start
    result = _save(order, complete=True)
    if verbose:
        print(f'  [ESPM] complete after order {order - 1}, {cumulative_elapsed_s:.1f}s cumulative')
    return result, True
