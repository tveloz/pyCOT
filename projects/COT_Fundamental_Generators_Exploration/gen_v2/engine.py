"""
engine.py -- gen_v2's search engine, REBUILT (2026-10-04) around the
corrected architecture worked out with Tomas, replacing the first
attempt entirely (that attempt's two bugs, both traced to the same root
cause, are kept only in git history / memory notes, not in this file).

The corrected architecture, in one sentence: the frontier (which moves --
synergy, complementarity, lift -- are available next) is a function of
the CLOSED SET alone, never of which particular generator reached it.

Why this fixes what the first attempt got wrong
-------------------------------------------------
The same closed set can have more than one valid generator (the paper
itself says generators are not unique). The first attempt asked "what
does the GENERATOR I happen to have on record offer next" -- so whichever
generator the search found first for a given closed set determined, by
accident of search order, what could be discovered from there. On one
real network this silently lost 2 genuine semi-organizations: two
different generators built the identical intermediate closed set, only
one of them had a hierarchy climb available from it, and the search kept
the wrong one.

This version never asks the generator anything. For ANY closed set it
reaches, it recomputes "which ERCs are contained in this closed set" (via
present_ercs, below) directly from the species bitmask -- so two
generators reaching the same destination produce, by construction, the
exact same frontier, and the ambiguity cannot arise.

Lift and the fundamental relations read the closed set differently:
  - Lift only makes sense from the ERCs currently "on top" (maximal by
    containment among everything present) -- a smaller ERC whose bigger
    relative is already present has nothing new to offer by climbing,
    since its own direct parent is already covered.
  - Fundamental relations are typically witnessed by LOWER nodes
    (fundamentality is defined at the minimal witness, Def 26-27 of the
    paper), so that check scans every ERC present, not just the maximal
    ones.

Move order, as specified: while a closed set still has an unmet
requirement, try synergy before complementarity. Once it is already a
semi-organization and being grown further, try complementarity before
synergy (lots of complementarity steps from an already-reached
semi-organization land directly on another one), with lift as its own
third category, available only at this stage.

A generator is no longer load-bearing for correctness -- kept only as a
reporting trail (one example construction path per closed set, via
parent_of), matching the paper's own Definition 29 (an ordered partition
of extension steps), never consulted to decide what is reachable.

Memoization: the search space of closed sets is a DAG (species only ever
get added, never removed -- no cycles). Each DISTINCT closed set is
resolved -- its presence set computed, its frontier tried -- EXACTLY
ONCE, no matter how many different generators/paths would otherwise
rediscover it. This is the safe version of the pruning Tomas described:
nothing is ever skipped ahead of time as "probably redundant" -- every
candidate from every closed set is tried, but landing on a closed set
that's already fully resolved elsewhere is an O(1) lookup, not a
re-derivation. Implemented as an explicit-stack iterative DFS (not
Python recursion), to avoid recursion-depth limits on networks with deep
exploration chains.

What is reused, unmodified, from the production engine
---------------------------------------------------------
FundamentalGraph itself (species_mask, req_mask, prod_mask, hierarchy
ancestors/descendants/parents, syn_from, comp_by_species,
comp_consumers_by_species) and erc_syn_close. erc_syn_close in particular
needed NO changes at all: it was never the source of either bug found
earlier today -- it already only ever checks species coverage, never
generator membership, which is exactly the property now applied
consistently everywhere else too. It is used here exactly as production
uses it: given a closed set (a plain species bitmask) and one explicit
ERC to add, it returns the properly-closed resulting species set,
chasing any further fundamental-synergy chain reactions that follow.

What is gone, and why
-----------------------
- ClosureCompleteGraph (the ancestor-inheritance patch from the first
  attempt): no longer needed at all. "Which ERCs are present" is now
  computed directly from the closed set (present_ercs, below) rather than
  approximated by pre-augmenting a lookup table -- which is exactly what
  both of today's bugs traced back to.
- The min_ext/min_seed canonical-ordering and Pareto-dominance tie
  machinery: superseded by closed-set memoization, which de-duplicates
  more simply and directly than an artificial constraint on the order
  ERCs may be added in -- and was itself the proximate mechanism behind
  today's bug (two generators tied on a closed set, one arbitrarily
  "winning").

Public API
----------
explore(g: FundamentalGraph) -> ExploreResult
  all_so_masks, elementary_masks, so_by_order, stats, parent_of, complete

Checkpoint / resume (added 2026-10-04)
---------------------------------------
Genome-scale networks (iAB_RBC_283 and similar, 600+ reactions) can
legitimately need more wall-clock time than any single run should be
given. Rather than restart such a run from nothing every time, explore()
can be given a checkpoint_path: if a checkpoint already exists there, the
search resumes from it (memo/discovered_sos/parent_of/stats restored
exactly as a prior run left them); it is also saved periodically while
running (checkpoint_every_s) and, if time_budget_s is given, once more
right before stopping early.

This is safe because memo (see explore(), below) only ever holds CLOSED
SETS THAT ARE FULLY RESOLVED -- resolve_from's explicit stack is a
post-order traversal, so a closed set only gets written to memo after
every closed set reachable from it (through this seed) already has been.
A checkpoint taken at any moment -- between seeds, mid-seed on a timer,
or mid-seed because an external kill cut the run -- is therefore always
a safe, self-consistent subset of the final answer: nothing in it is
half-finished, and nothing genuine gets lost, only possibly redone a
second time (cheaply, since its own already-resolved descendants are
already sitting in memo).
"""
from __future__ import annotations

import os
import pickle
import time
from collections import deque
from dataclasses import dataclass, field

from pyCOT.analysis.organizations.fundamental_graph import FundamentalGraph


# ---------------------------------------------------------------------------
# TraceLogger -- staged progress + rolling window of the last N resolved
# closed sets, for watching a search that is still running. Requested
# explicitly: a marker at each stage, an immediate uncapped print the
# moment any new semi-organization is found, and a bounded-memory trace of
# recent activity refreshed periodically, so "where is it right now" is
# answerable from the log at any point during a long run.
# ---------------------------------------------------------------------------

class TraceLogger:
    def __init__(self, n_log: int = 30, print_every: int = 500, network_label: str = ""):
        self.n_log = n_log
        self.print_every = print_every
        self.network_label = network_label
        self.buf: deque = deque(maxlen=n_log)
        self.step_count = 0
        self.so_count = 0
        self.t0 = time.perf_counter()

    def _prefix(self) -> str:
        return f"[{self.network_label}] " if self.network_label else ""

    def stage(self, msg: str) -> None:
        print(f"{self._prefix()}[STAGE t={time.perf_counter() - self.t0:7.1f}s] {msg}", flush=True)

    def resolved(self, S: int, is_so: bool, n_candidates: int, n_present: int,
                 depth: int, req_bits: int) -> None:
        self.step_count += 1
        req_info = "is SO" if is_so else f"req: {req_bits} species still missing"
        self.buf.append(
            f"{self._prefix()}step#{self.step_count:>8}  {'SO ' if is_so else '   '}  "
            f"{n_present:>4} ERCs present  {bin(S).count('1'):>4} species  "
            f"depth={depth:>3} (additions from seed)  {n_candidates:>3} candidates tried  {req_info}"
        )
        if self.step_count % self.print_every == 0:
            self.flush()

    def flush(self) -> None:
        elapsed = time.perf_counter() - self.t0
        print(f"{self._prefix()}  -- progress: {self.step_count} closed sets resolved, "
              f"{self.so_count} SOs found, t={elapsed:.1f}s -- last {len(self.buf)}:", flush=True)
        for line in self.buf:
            print("    " + line, flush=True)

    def landmark_so(self, S: int, n_present: int, depth: int) -> None:
        self.so_count += 1
        elapsed = time.perf_counter() - self.t0
        print(f"{self._prefix()}  *** NEW SO #{self.so_count:<5} {n_present} ERCs, "
              f"{bin(S).count('1')} species, reached via {depth} additions from its seed "
              f"(t={elapsed:.1f}s, step#{self.step_count})", flush=True)


def _bits(mask: int):
    m = mask
    while m:
        lsb = m & (-m)
        yield lsb.bit_length() - 1
        m &= m - 1


# ---------------------------------------------------------------------------
# Closed-set-dependent presence and requirement/production
# ---------------------------------------------------------------------------

def present_ercs(g: FundamentalGraph, S: int) -> list[int]:
    """
    Which ERCs are contained in the closed set S (species_mask[i] subset
    of S), computed directly from S -- never from a generator.

    CORRECTNESS NOTE (found 2026-10-04, while retesting BIOMD0000000109):
    an earlier version of this also returned a "maximal_present" subset
    (ERCs not properly contained in any OTHER present ERC), on the
    assumption that only those could have a meaningful lift move -- "a
    smaller ERC on the same hierarchy chain as an already-present bigger
    one has nothing new to climb to." That assumption is FALSE in
    general: the ERC hierarchy is not one chain per node. A single ERC
    can have several direct parents that are mutually INCOMPARABLE to
    each other (confirmed directly on this network: ERC3's direct
    parents include ERC4 and ERC12, and ERC4/ERC12 do not contain one
    another). So being dominated by one present ancestor (ERC12) says
    nothing about whether a DIFFERENT branch (via ERC4) is still open --
    filtering to "maximal" silently discarded exactly that branch and
    lost 2 genuine semi-organizations. Lift candidates are now computed
    from every present ERC's own direct parents (see lift_candidates),
    each individually checked for whether it's already covered -- which
    was already correct at the per-candidate level; the bug was only in
    pre-filtering which ERCs got asked at all.
    """
    return [i for i in range(g.n) if (g.species_mask[i] & ~S) == 0]


def req_prod_of(g: FundamentalGraph, all_present: list[int]) -> tuple[int, int]:
    """
    req(S) and prod(S), computed as the union of req_mask/prod_mask over
    every ERC contained in S.

    This is exact, not an approximation: S is always closed in this
    engine (every state reached is produced by erc_syn_close, which
    always returns a properly-closed species set), so every reaction
    active within S has support reducing to some present ERC (r active in
    S means ERC(r) = closure(supp(r)) is itself closed and subset of S,
    since S is closed and contains supp(r) -- so ERC(r) is present), and
    conversely every reaction in R_i for a present ERC i is active in S
    (supp(r) subset of species_mask[i] subset of S). So the union over
    present ERCs' own req/prod masks (each already correct for that ERC's
    own closure, by construction in erc.py) is exactly req(S)/prod(S) --
    with no dependency at all on which generator built S.
    """
    prod = 0
    req_raw = 0
    for i in all_present:
        prod |= g.prod_mask[i]
        req_raw |= g.req_mask[i]
    return req_raw & ~prod, prod


# ---------------------------------------------------------------------------
# Frontiers -- both computed from the closed set (via present_ercs), never
# from a generator.
# ---------------------------------------------------------------------------

def fundamental_candidates(g: FundamentalGraph, S: int, all_present: list[int],
                            req: int, prod: int, is_so: bool) -> list[int]:
    """
    Synergy and complementarity candidates, in the specified move order:
    synergy-before-complementarity while still building (is_so=False);
    complementarity-before-synergy once already at a semi-organization
    (is_so=True).
    """
    syn_cand: dict[int, None] = {}
    for i in all_present:
        for (j, _k) in g.syn_from.get(i, []):
            if (g.species_mask[j] & ~S) == 0:
                continue  # j already fully present -- not a new move
            syn_cand.setdefault(j, None)

    comp_cand: dict[int, None] = {}
    if not is_so:
        # producer side (Mode-1-style): minimal producers of missing species
        for s_bit in _bits(req):
            for prod_idx in g.comp_by_species.get(s_bit, []):
                if (g.species_mask[prod_idx] & ~S) == 0:
                    continue
                comp_cand.setdefault(prod_idx, None)
    else:
        # consumer side (Mode-2-style): minimal consumers of species already produced
        for s_bit in _bits(prod):
            for cons_idx in g.comp_consumers_by_species.get(s_bit, []):
                if (g.species_mask[cons_idx] & ~S) == 0:
                    continue
                comp_cand.setdefault(cons_idx, None)

    if is_so:
        return list(comp_cand) + [j for j in syn_cand if j not in comp_cand]
    return list(syn_cand) + [j for j in comp_cand if j not in syn_cand]


def lift_candidates(g: FundamentalGraph, S: int, all_present: list[int]) -> list[int]:
    """
    Direct hierarchy parents of every present ERC -- only called once the
    closed set is already a semi-organization.

    Iterates ALL present ERCs, not a "maximal" subset (see present_ercs's
    docstring for why that filtering was wrong: a present ERC can have
    several mutually incomparable direct parents, so being dominated by
    one does not mean every branch is covered). The per-candidate check
    below ('already covered, skip') is what correctly avoids proposing a
    redundant climb when some direct parent IS already present -- that
    check was always correct; only the ERC-level pre-filtering was not.
    """
    cand: dict[int, None] = {}
    for i in all_present:
        for a in g.parents[i]:
            if (g.species_mask[a] & ~S) == 0:
                continue
            cand.setdefault(a, None)
    return list(cand)


# ---------------------------------------------------------------------------
# Checkpoint / resume -- see the module docstring for why this is safe.
# ---------------------------------------------------------------------------

@dataclass
class CheckpointState:
    memo: dict
    discovered_sos: set
    parent_of: dict
    stats: dict
    next_seed_index: int    # resume the seed loop from here
    complete: bool = False  # True only once every seed has been processed
    network_label: str = ""  # e.g. "iAB_RBC_283 (645 rxn, 161 ERCs)" -- so a cached-complete
                              # checkpoint can be reported on without rebuilding anything to
                              # re-derive it. Older checkpoints saved before this field existed
                              # simply won't have the attribute at all after unpickling -- always
                              # read it via getattr(state, "network_label", "") defensively.


def save_checkpoint(path: str, state: CheckpointState) -> None:
    """Atomic write (temp file + os.replace) so a kill mid-write can never
    leave a corrupt checkpoint behind -- the previous good one stays in
    place until the new one is fully written.

    On Windows, a synced folder (this project lives under Dropbox) can
    transiently hold its own read lock on the destination file at the
    exact moment of the rename, raising PermissionError / WinError 5 for
    no real reason -- retried with a short backoff rather than treated as
    a hard failure on the first try."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp_path = path + ".tmp"
    with open(tmp_path, "wb") as f:
        pickle.dump(state, f, protocol=pickle.HIGHEST_PROTOCOL)
    last_exc: OSError | None = None
    for attempt in range(8):
        try:
            os.replace(tmp_path, path)
            return
        except OSError as exc:
            last_exc = exc
            time.sleep(0.25 * (attempt + 1))
    raise last_exc


def load_checkpoint(path: str) -> CheckpointState | None:
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


class _TimeBudgetExceeded(Exception):
    pass


# ---------------------------------------------------------------------------
# Memoized exploration
# ---------------------------------------------------------------------------

@dataclass
class ExploreResult:
    all_so_masks: list
    elementary_masks: list
    so_by_order: dict
    stats: dict = field(default_factory=dict)
    parent_of: dict = field(default_factory=dict)   # sp -> (prev_sp, move_erc_idx)
    complete: bool = True   # False if a time budget cut the search short


def _assign_orders(all_sp_masks: list[int]) -> dict[int, int]:
    sorted_masks = sorted(all_sp_masks, key=lambda m: bin(m).count('1'))
    so_order: dict[int, int] = {}
    for cl in sorted_masks:
        max_sub = -1
        for sub_cl, sub_ord in so_order.items():
            if sub_cl != cl and (sub_cl & cl) == sub_cl and sub_ord > max_sub:
                max_sub = sub_ord
        so_order[cl] = max_sub + 1
    return so_order


def build_result(discovered_sos: set, stats: dict, parent_of: dict, complete: bool,
                  logger: TraceLogger | None = None) -> ExploreResult:
    """
    Public on purpose: builds an ExploreResult purely from already-known
    state (discovered_sos/stats/parent_of), independent of HOW that state
    was obtained. explore() calls this internally, but a caller that only
    wants to know "is this network already fully done, and if so what did
    it find" can call load_checkpoint() + build_result() directly, with
    complete=True, and never construct a FundamentalGraph or run any of
    the ERC/hierarchy/synergy/complementarity pipeline at all.
    """
    all_sp = sorted(discovered_sos)
    so_order = _assign_orders(all_sp)
    so_by_order: dict[int, list[int]] = {}
    for sp, k in so_order.items():
        so_by_order.setdefault(k, []).append(sp)
    for k in so_by_order:
        so_by_order[k].sort()
    elementary = sorted(sp for sp, k in so_order.items() if k == 0)

    if logger:
        logger.stage(f"{'Done' if complete else 'Stopped (incomplete)'} -- "
                     f"{len(all_sp)} total semi-organizations, {len(elementary)} elementary, "
                     f"{max(so_by_order.keys(), default=0)} higher orders")

    return ExploreResult(
        all_so_masks=all_sp,
        elementary_masks=elementary,
        so_by_order=so_by_order,
        stats=stats,
        parent_of=parent_of,
        complete=complete,
    )


def explore(g: FundamentalGraph, *, verbose: bool = False,
            n_log: int = 30, print_every: int = 500,
            checkpoint_path: str | None = None,
            time_budget_s: float | None = None,
            checkpoint_every_s: float = 30.0,
            network_label: str = "") -> ExploreResult:
    """
    Resolve every closed set reachable (via fundamental synergy,
    complementarity, and -- once at a semi-organization -- lift) starting
    from every single ERC's own closure, memoizing each distinct closed
    set exactly once.

    checkpoint_path, time_budget_s, checkpoint_every_s: see the module
    docstring's "Checkpoint / resume" section. All three are optional and
    off by default -- a plain explore(g) call behaves exactly as before.

    network_label: prefixed onto every verbose log line (e.g.
    "iAB_RBC_283 (645 rxn, 161 ERCs)") and saved into the checkpoint, so a
    long, scrolled-past log never loses track of which network a line
    belongs to, and a cached-complete checkpoint can report it back
    without recomputing anything.
    """
    logger = TraceLogger(n_log=n_log, print_every=print_every, network_label=network_label) if verbose else None

    memo: dict[int, bool] = {}       # closed set -> is_so (presence = fully resolved)
    pending: dict[int, tuple] = {}   # closed set -> (is_so, candidates, n_present, req_bits), between visits
    discovered_sos: set[int] = set()
    parent_of: dict[int, tuple] = {}
    depth: dict[int, int] = {}       # closed set -> number of fundamental/lift additions from its seed
    stats = {"states_resolved": 0, "presence_scans": 0}
    start_seed_index = 0

    if checkpoint_path is not None:
        loaded = load_checkpoint(checkpoint_path)
        if loaded is not None:
            memo, discovered_sos, parent_of, stats = (
                loaded.memo, loaded.discovered_sos, loaded.parent_of, loaded.stats)
            if loaded.complete:
                if logger:
                    logger.stage(f"Checkpoint already complete ({stats['states_resolved']} resolved, "
                                 f"{len(discovered_sos)} SOs) -- returning cached result, no work to do")
                return build_result(discovered_sos, stats, parent_of, complete=True, logger=logger)
            start_seed_index = loaded.next_seed_index
            if logger:
                logger.stage(f"Resuming from checkpoint: {stats['states_resolved']} closed sets "
                             f"already resolved, {len(discovered_sos)} SOs already found, "
                             f"continuing from seed {start_seed_index}/{g.n}")

    t_start = time.perf_counter()
    last_checkpoint_t = t_start

    def maybe_checkpoint(next_seed_index: int, force: bool = False) -> None:
        nonlocal last_checkpoint_t
        if checkpoint_path is None:
            return
        now = time.perf_counter()
        if not force and (now - last_checkpoint_t) < checkpoint_every_s:
            return
        # A save failure here (e.g. a synced folder's transient file lock,
        # even after save_checkpoint's own retries) must never crash a
        # long-running search -- the in-memory state is still fine, this
        # periodic save is just best-effort; the next one tries again.
        try:
            save_checkpoint(checkpoint_path, CheckpointState(
                memo=memo, discovered_sos=discovered_sos, parent_of=parent_of,
                stats=stats, next_seed_index=next_seed_index, complete=False,
                network_label=network_label,
            ))
            last_checkpoint_t = now
            if logger:
                logger.stage(f"Checkpoint saved -- {stats['states_resolved']} resolved, "
                             f"{len(discovered_sos)} SOs, resume point = seed {next_seed_index}/{g.n}")
        except OSError as exc:
            msg = f"Checkpoint save FAILED (continuing without it): {exc}"
            if logger:
                logger.stage(msg)
            else:
                print(msg, flush=True)
            last_checkpoint_t = now  # don't retry every single state until checkpoint_every_s passes again

    def resolve_from(seed_sp: int, current_seed_index: int) -> None:
        if seed_sp in memo:
            return
        depth.setdefault(seed_sp, 0)
        stack: list[tuple[int, bool]] = [(seed_sp, False)]
        in_progress: set[int] = set()

        while stack:
            S, expanded = stack.pop()
            if S in memo:
                continue

            if not expanded:
                if S in in_progress:
                    continue
                in_progress.add(S)

                all_present = present_ercs(g, S)
                stats["presence_scans"] += 1
                req, prod = req_prod_of(g, all_present)
                is_so = (req == 0)

                cand = fundamental_candidates(g, S, all_present, req, prod, is_so)
                if is_so:
                    cand = cand + lift_candidates(g, S, all_present)

                pending[S] = (is_so, cand, len(all_present), bin(req).count("1"))
                stack.append((S, True))

                for w in cand:
                    _implied, S_next = g.erc_syn_close(S, [w])
                    if S_next not in memo:
                        # parent_of is cross-run (loaded from a resumed
                        # checkpoint, so it may already have an entry here
                        # from a PRIOR run) -- depth is this-run-only (never
                        # persisted, see explore()'s docstring), so it needs
                        # its OWN guard. Reusing parent_of's guard for depth
                        # would skip setting it whenever resuming visits a
                        # state the previous run had already discovered but
                        # not yet resolved, leaving depth[S_next] unset and
                        # crashing the next level down with a KeyError.
                        if S_next not in parent_of:
                            parent_of[S_next] = (S, w)
                        if S_next not in depth:
                            depth[S_next] = depth.get(S, 0) + 1
                        if S_next not in in_progress:
                            stack.append((S_next, False))
            else:
                is_so, cand, n_present, req_bits = pending.pop(S)
                d = depth.get(S, 0)
                if is_so:
                    discovered_sos.add(S)
                    if logger is not None:
                        logger.landmark_so(S, n_present, d)
                if logger is not None:
                    logger.resolved(S, is_so, len(cand), n_present, d, req_bits)
                memo[S] = is_so
                stats["states_resolved"] += 1

                if stats["states_resolved"] % 256 == 0:
                    if time_budget_s is not None and (time.perf_counter() - t_start) > time_budget_s:
                        maybe_checkpoint(current_seed_index, force=True)
                        raise _TimeBudgetExceeded()
                    maybe_checkpoint(current_seed_index)

    complete = True
    if logger:
        logger.stage(f"Exploration: starting, {g.n} ERC seeds"
                     + (f" (resuming from seed {start_seed_index})" if start_seed_index else ""))
    try:
        for i in range(start_seed_index, g.n):
            if logger and g.species_mask[i] not in memo:
                logger.stage(f"Exploration: seed E{i} ({i + 1}/{g.n}) -- "
                             f"{stats['states_resolved']} resolved so far, "
                             f"{len(discovered_sos)} SOs found so far")
            resolve_from(g.species_mask[i], i)
            maybe_checkpoint(i + 1)
    except _TimeBudgetExceeded:
        complete = False
        if logger:
            logger.stage(f"Exploration: time budget ({time_budget_s}s) exceeded -- "
                         f"checkpoint saved, {stats['states_resolved']} resolved so far")
    if logger:
        logger.flush()

    if complete and checkpoint_path is not None:
        # Same non-fatal handling as maybe_checkpoint: the search itself
        # already fully succeeded (everything below depends only on
        # in-memory state) -- a failure to persist that MUST NOT cost the
        # caller the result they already have.
        try:
            save_checkpoint(checkpoint_path, CheckpointState(
                memo=memo, discovered_sos=discovered_sos, parent_of=parent_of,
                stats=stats, next_seed_index=g.n, complete=True,
                network_label=network_label,
            ))
        except OSError as exc:
            msg = f"Final checkpoint save FAILED (result below is still correct, just not cached): {exc}"
            if logger:
                logger.stage(msg)
            else:
                print(msg, flush=True)

    return build_result(discovered_sos, stats, parent_of, complete=complete, logger=logger)
