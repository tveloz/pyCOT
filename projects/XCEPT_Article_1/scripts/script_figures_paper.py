#!/usr/bin/env python3
"""
script_figures_paper.py
=======================
Generates four figures for the XCEPT Article 1 paper.
Uses the game-theoretic discrete-time kinetics (Reactive vs Reactive).

Figures saved to:  projects/XCEPT_Article_1/figures/
  loop_projections.png   – species time series + conflict/peace loop shadows
  timescale_loops.png    – loop shadows at per-step / 4-step / 12-step resolution
  two_paths.png          – LP: two solution paths in (A, I) space
  iterative_lp.png       – iterative LP: trajectory + budget allocation over 12 periods
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import linprog

# ── pyCOT path + STOICH ───────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))
from pyCOT.io.functions import read_txt

SPECIES   = ['G', 'A', 'I', 'T', 'Rc', 'Rp']
REACTIONS = ['u1', 'u2', 'u3', 'u4', 'c1', 'c2', 'c3', 'c4', 'f1', 'f2']

_RN_FILE  = os.path.join(_SCRIPT_DIR, 'data', 'Basic_example.txt')
_rn       = read_txt(_RN_FILE)
_sm       = _rn.stoichiometry_matrix()
_pycot_sp = list(_sm.species)
_pycot_rx = list(_sm.reactions)
_row_idx  = [_pycot_sp.index(s) for s in SPECIES]
_col_idx  = [_pycot_rx.index(r) for r in REACTIONS]
STOICH    = np.array(_sm, dtype=float)[np.ix_(_row_idx, _col_idx)]  # (6, 10)

# Reaction index shortcuts
U1, U2, U3, U4 = 0, 1, 2, 3
C1, C2, C3, C4 = 4, 5, 6, 7
F1, F2         = 8, 9
CONF_LOOP_IDX  = [U1, U2, U3]
PEACE_LOOP_IDX = [C1, C2, C3]
ACTIVE_IDX     = [U1, U2, U3, U4, C1, C2, C3, C4]

# Species index shortcuts
iG, iA, iI, iT, iRc, iRp = 0, 1, 2, 3, 4, 5

# Output directory
FIGURES_DIR = os.path.join(os.path.dirname(_SCRIPT_DIR), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── Colour palette ────────────────────────────────────────────────────────────
SP_COLORS = ['#d62728', '#ff7f0e', '#2ca02c', '#17becf', '#9467bd', '#1f77b4']
# G=red, A=orange, I=green, T=cyan, Rc=purple, Rp=blue

# ===========================================================================
# SIMULATION MACHINERY  (copied from script_basic_game.py)
# ===========================================================================

DEFAULT_K_MAX        = np.array([20.0, 15.0, 15.0, 15.0, 15.0, 20.0])
DEFAULT_K_AUTO       = dict(k_f1=0.001, k_f2=0.001)
DEFAULT_CONF_THR     = np.array([2.0, 1.0])
DEFAULT_CONF_PCT     = np.array([[0.3, 0.7], [0.5, 0.9]])
DEFAULT_PEACE_THR    = np.array([3.0, 2.0])
DEFAULT_PEACE_PCT    = np.array([[0.3, 0.7], [0.4, 0.8]])


def auto_decay(x, k_f1=0.001, k_f2=0.001):
    v = np.zeros(10)
    v[F1] = min(k_f1 * x[iG] ** 2, x[iG] / 2.0)
    v[F2] = min(k_f2 * x[iT] ** 2, x[iT] / 2.0)
    return v


def _investable(stock, threshold, pct_interval, rng):
    if stock <= threshold:
        return stock
    return rng.uniform(pct_interval[0], pct_interval[1]) * stock


def conflict_budgets(x, thr, pct, rng):
    return (_investable(x[iA],  thr[0], pct[0], rng),
            _investable(x[iRc], thr[1], pct[1], rng))


def peace_budgets(x, thr, pct, rng):
    return (_investable(x[iRp], thr[0], pct[0], rng),
            _investable(x[iI],  thr[1], pct[1], rng))


def cap_conflict(x, v_u1, v_u2, v_u3, v_u4):
    return (min(v_u1, x[iRp]),
            v_u2,
            min(v_u3, x[iG], x[iRc]),
            min(v_u4, x[iI]))


def cap_peace(x, v_c1, v_c2, v_c3, v_c4):
    return (min(v_c1, x[iT], x[iRp]),
            v_c2,
            min(v_c3, x[iRp]),
            min(v_c4, x[iA]))


def alloc_conflict_reactive(x, bA, bRc):
    eps = 1e-6
    w1 = 1.0 / (x[iRc] + eps)
    w2 = 1.0 / (x[iG]  + eps)
    w4 = x[iI] / (x[iA] + eps)
    tot = w1 + w2 + w4 + eps
    return cap_conflict(x, bA*w1/tot, bA*w2/tot, bRc, bA*w4/tot)


def alloc_peace_reactive(x, bRp, bI):
    eps = 1e-6
    w_c1 = 1.0 / (x[iRp] + eps)
    w_c3 = 1.0 / (x[iI]  + eps)
    tot_rp = w_c1 + w_c3 + eps
    w_c2 = 1.0 / (x[iT] + eps)
    w_c4 = x[iA] / (x[iI] + eps)
    tot_i = w_c2 + w_c4 + eps
    return cap_peace(x, bRp*w_c1/tot_rp, bI*w_c2/tot_i,
                        bRp*w_c3/tot_rp, bI*w_c4/tot_i)


def apply_ceilings(x, v, k_max):
    v = v.copy()
    for j in range(10):
        scale = 1.0
        for s in range(6):
            if STOICH[s, j] > 0:
                scale = min(scale, max(0.0, 1.0 - x[s] / k_max[s]))
        v[j] *= scale
    return v


def resolve_contention(x, v, max_iter=20):
    v = v.copy()
    for _ in range(max_iter):
        delta = STOICH @ v
        x_new = x + delta
        if np.all(x_new >= -1e-10):
            break
        for s in range(6):
            if x_new[s] < 0:
                cons = STOICH[s] < 0
                if cons.any():
                    drain = -(STOICH[s] @ v)
                    if drain > 1e-12:
                        v[cons] *= x[s] / drain
                        x_new = x + STOICH @ v
    return v


def step(x, conf_thr, conf_pct, peace_thr, peace_pct, k_auto, k_max, rng):
    v = auto_decay(x, **k_auto)
    bA, bRc = conflict_budgets(x, conf_thr, conf_pct, rng)
    v_u1, v_u2, v_u3, v_u4 = alloc_conflict_reactive(x, bA, bRc)
    v[U1], v[U2], v[U3], v[U4] = v_u1, v_u2, v_u3, v_u4
    bRp, bI = peace_budgets(x, peace_thr, peace_pct, rng)
    v_c1, v_c2, v_c3, v_c4 = alloc_peace_reactive(x, bRp, bI)
    v[C1], v[C2], v[C3], v[C4] = v_c1, v_c2, v_c3, v_c4
    v = apply_ceilings(x, v, k_max)
    v = resolve_contention(x, v)
    return np.maximum(x + STOICH @ v, 0.0), v


def run_simulation(n_steps, x0, seed=42,
                   conf_thr=None, conf_pct=None,
                   peace_thr=None, peace_pct=None,
                   k_auto=None, k_max=None):
    if conf_thr   is None: conf_thr   = DEFAULT_CONF_THR
    if conf_pct   is None: conf_pct   = DEFAULT_CONF_PCT
    if peace_thr  is None: peace_thr  = DEFAULT_PEACE_THR
    if peace_pct  is None: peace_pct  = DEFAULT_PEACE_PCT
    if k_auto     is None: k_auto     = DEFAULT_K_AUTO
    if k_max      is None: k_max      = DEFAULT_K_MAX
    rng = np.random.default_rng(seed)
    X = np.zeros((n_steps + 1, 6))
    V = np.zeros((n_steps, 10))
    X[0] = x0.copy()
    for t in range(n_steps):
        X[t+1], V[t] = step(X[t], conf_thr, conf_pct,
                             peace_thr, peace_pct, k_auto, k_max, rng)
    return X, V


# ===========================================================================
# LOOP SHADOW HELPERS
# ===========================================================================

def loop_shadow(V, loop_idx, d_coeff=None):
    """
    Per-step bottleneck projection π_d(v) = min_{j: d_j>0} v_j/d_j
    for the given loop reactions.
    d_coeff: array of mode coefficients for the loop reactions (default: all ones).
    Returns arrays: m (min rate), M (max rate), E (mean rate),
    and shadow = π_d(v) (bottleneck projection).
    """
    v_loop = V[:, loop_idx]                                # (T, k)
    if d_coeff is None:
        d_coeff = np.ones(len(loop_idx), dtype=float)
    else:
        d_coeff = np.asarray(d_coeff, dtype=float)
    shadow = np.min(v_loop / d_coeff, axis=1)              # (T,) — bottleneck
    m = v_loop.min(axis=1)
    M = v_loop.max(axis=1)
    E = v_loop.mean(axis=1)
    return m, M, E, shadow


def windowed_mean(arr, w):
    """Mean of arr over non-overlapping windows of size w. Returns (centers, values)."""
    T = len(arr)
    n = T // w
    vals    = np.array([arr[i*w:(i+1)*w].mean() for i in range(n)])
    centers = np.array([(i*w + w/2) for i in range(n)])
    return centers, vals


# ===========================================================================
# FIGURE 1 — loop_projections.png
# ===========================================================================

def fig_loop_projections(X, V, save_path):
    T   = len(V)
    t   = np.arange(T)

    # C2: d=[1,1,1]  →  π = min(v_u1, v_u2, v_u3)
    # P3: d=[2,2,1]  →  π = min(v_c1/2, v_c2/2, v_c3)
    mC, MC, EC, shC = loop_shadow(V, CONF_LOOP_IDX, d_coeff=[1, 1, 1])
    mP, MP, EP, shP = loop_shadow(V, PEACE_LOOP_IDX, d_coeff=[2, 2, 1])

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                             gridspec_kw={'height_ratios': [1.2, 1]})

    # ── Top panel: species trajectories ──────────────────────────────────────
    ax = axes[0]
    for i, (sp, col) in enumerate(zip(SPECIES, SP_COLORS)):
        ax.plot(t, X[:T, i], color=col, lw=1.5, label=sp)
    ax.set_ylabel('Population', fontsize=10)
    ax.legend(fontsize=8, ncol=6, loc='upper right')
    ax.grid(True, alpha=0.2)
    ax.set_title('Species trajectories', fontsize=9)

    # ── Bottom panel: loop min / mean / max + shading ────────────────────────
    ax = axes[1]

    # Conflict loop: red tones
    ax.fill_between(t, mC, MC, color='#FFCCCC', alpha=0.4)
    ax.plot(t, mC, color='#FF4444', lw=1.0, ls='--', alpha=0.7, label='m_C (min)')
    ax.plot(t, EC, color='#CC0000', lw=1.8,           label='E_C (mean)')
    ax.plot(t, MC, color='#FF9999', lw=1.0, ls=':',  alpha=0.7, label='M_C (max)')

    # Peace loop: blue tones
    ax.fill_between(t, mP, MP, color='#CCDEFF', alpha=0.4)
    ax.plot(t, mP, color='#4499FF', lw=1.0, ls='--', alpha=0.7, label='m_P (min)')
    ax.plot(t, EP, color='#0044CC', lw=1.8,           label='E_P (mean)')
    ax.plot(t, MP, color='#99CCFF', lw=1.0, ls=':',  alpha=0.7, label='M_P (max)')

    # Dominance shading
    conf_dom = mC > MP   # conflict minimum exceeds peace maximum
    peace_dom = mP > MC  # peace minimum exceeds conflict maximum
    ax.fill_between(t, 0, 1, where=conf_dom,
                    transform=ax.get_xaxis_transform(),
                    color='#FF0000', alpha=0.12, label='Conflict-dominated')
    ax.fill_between(t, 0, 1, where=peace_dom,
                    transform=ax.get_xaxis_transform(),
                    color='#0000FF', alpha=0.12, label='Peace-dominated')

    ax.set_ylabel('Rate / bottleneck projection', fontsize=10)
    ax.set_xlabel('Time step', fontsize=10)
    ax.legend(fontsize=7, ncol=4, loc='upper right')
    ax.grid(True, alpha=0.2)
    ax.set_title('C2 (u1–u2–u3, d=[1,1,1]) vs P3 (c1–c2–c3, d=[2,2,1])', fontsize=9)

    fig.suptitle('Loop bottleneck projections  π_d(v) = min v_j/d_j',
                 fontsize=10, fontweight='bold')
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {save_path}')


# ===========================================================================
# FIGURE 2 — timescale_loops.png
# ===========================================================================

def fig_timescale_loops(V, save_path):
    windows = [1, 4, 12]
    titles  = ['Per step', '4-step windows', '12-step windows']

    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=False)

    for ax, w, title in zip(axes, windows, titles):
        if w == 1:
            t_c = np.arange(len(V))
            ec  = V[:, CONF_LOOP_IDX].mean(axis=1)
            ep  = V[:, PEACE_LOOP_IDX].mean(axis=1)
        else:
            # Sum v over non-overlapping windows, then take mean across the 3 loop reactions
            T  = len(V)
            n  = T // w
            ec = np.array([V[i*w:(i+1)*w, CONF_LOOP_IDX].sum(axis=0).mean() for i in range(n)])
            ep = np.array([V[i*w:(i+1)*w, PEACE_LOOP_IDX].sum(axis=0).mean() for i in range(n)])
            t_c = np.array([(i*w + w/2) for i in range(n)])

        ax.plot(t_c, ec, color='#CC0000', lw=1.5, label='E_C (conflict)')
        ax.plot(t_c, ep, color='#0044CC', lw=1.5, label='E_P (peace)')
        ax.fill_between(t_c, ec, ep,
                        where=ec >= ep, color='#FFCCCC', alpha=0.35)
        ax.fill_between(t_c, ec, ep,
                        where=ep >  ec, color='#CCDEFF', alpha=0.35)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('Time step', fontsize=9)
        ax.set_ylabel('Mean loop rate', fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

    fig.suptitle('Conflict and peace loop activity at multiple time scales',
                 fontsize=10, fontweight='bold')
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {save_path}')


# ===========================================================================
# FIGURE 3 — two_paths.png
# ===========================================================================

def fig_two_paths(save_path):
    """
    Starting state x0 = (8,6,5,3,5,4): high G,A (conflict) but enough I,T,Rp
    for peace reactions to be kinetically active.
    Target: ΔG ≤ -2, ΔA ≤ -1.5, ΔRc ≤ -1 (achievable reductions in one period).
    Find two solution process vectors via LP with different objectives:
      Path 1 (confrontation): maximise v_c4 + v_f1
      Path 2 (peace-first):   maximise v_c1 + v_c2 + v_c3
    Split each into two phases, show trajectory in (A, I) space.
    """
    x0 = np.array([8.0, 6.0, 5.0, 3.0, 5.0, 4.0])  # G A I T Rc Rp

    conf_rows = [iG, iA, iRc]
    delta_tgt = np.array([-2.0, -1.5, -1.0])         # achievable at these x0 values

    k = dict(k_u1=0.05, k_u2=0.2, k_u3=0.05, k_u4=0.1,
             k_c1=0.1,  k_c2=0.25, k_c3=0.1, k_c4=0.1,
             k_f1=0.01, k_f2=0.01)
    G, A, I, T, Rc, Rp = x0
    v_ub = np.array([
        k['k_u1'] * A  * Rp,
        k['k_u2'] * A,
        k['k_u3'] * G  * Rc,
        k['k_u4'] * A  * I,
        k['k_c1'] * T  * Rp,
        k['k_c2'] * I,
        k['k_c3'] * I  * Rp,
        k['k_c4'] * I  * A,
        k['k_f1'] * G**2,
        k['k_f2'] * T**2,
    ])

    # LP: find v ≥ 0, S[conf_rows,:] @ v ≤ delta_tgt, v ≤ v_ub
    # Equality rows of interest in S:
    A_ub = STOICH[conf_rows, :]      # (3, 10)  S rows for G, A, Rc
    b_ub = delta_tgt                 # (3,)

    def solve_lp(objective_coeff):
        """Minimise objective_coeff @ v subject to constraints."""
        res = linprog(
            c      = objective_coeff,
            A_ub   = A_ub,
            b_ub   = b_ub,
            bounds = [(0, ub) for ub in v_ub],
            method = 'highs'
        )
        if not res.success:
            raise RuntimeError(f'LP failed: {res.message}')
        return res.x

    # Path 1: maximise v_c4 + v_f1  → minimise -(v_c4 + v_f1)
    c1_obj        = np.zeros(10)
    c1_obj[C4]    = -1.0
    c1_obj[F1]    = -1.0
    v1            = solve_lp(c1_obj)

    # Path 2: maximise v_c1 + v_c2 + v_c3  → minimise -(v_c1+v_c2+v_c3)
    c2_obj        = np.zeros(10)
    c2_obj[[C1, C2, C3]] = -1.0
    v2            = solve_lp(c2_obj)

    # Phase splits
    v1_ph1 = 0.6 * v1
    v1_ph2 = 0.4 * v1
    v2_ph1 = v2.copy(); v2_ph1[[U1,U2,U3,U4,F1,F2]] = 0.0   # peace-loop only
    v2_ph2 = v2.copy(); v2_ph2[[C1,C2,C3]]           = 0.0   # confrontation only

    def midpoint(x, v_phase):
        return np.maximum(x + STOICH @ v_phase, 0.0)

    x_mid1   = midpoint(x0, v1_ph1)
    x_mid2   = midpoint(x0, v2_ph1)
    x_target = midpoint(x0, v1)    # same endpoint for both (same total effect)

    k_c4 = k['k_c4']

    fig, ax = plt.subplots(figsize=(8, 7))

    # x0, x_target markers
    ax.scatter(*[x0[iA],     x0[iI]],     s=120, c='k',       zorder=5, label='x₀')
    ax.scatter(*[x_target[iA], x_target[iI]], s=180, c='green', marker='*',
               zorder=5, label='x_target')

    def draw_path(pts, color, ls, label):
        xs = [p[iA] for p in pts]
        ys = [p[iI] for p in pts]
        ax.plot(xs, ys, color=color, ls=ls, lw=2.0, label=label)
        for i in range(len(pts)-1):
            dx = xs[i+1]-xs[i]; dy = ys[i+1]-ys[i]
            ax.annotate('', xy=(xs[i+1], ys[i+1]), xytext=(xs[i], ys[i]),
                        arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

    draw_path([x0, x_mid1, x_target], '#CC0000', '--', 'Path 1 (confrontation)')
    draw_path([x0, x_mid2, x_target], '#0044CC', '-',  'Path 2 (peace-first)')

    # Circles proportional to kinetic capacity for c4 at each midpoint
    for xm, label_txt, color in [
        (x_mid1, 'Low I: c4 constrained',  '#CC0000'),
        (x_mid2, 'High I: c4 accessible',  '#0044CC'),
    ]:
        radius = k_c4 * xm[iI] * xm[iA] * 0.15   # scale for visibility
        circle = plt.Circle((xm[iA], xm[iI]), radius,
                             color=color, fill=False, lw=1.5, ls=':', zorder=4)
        ax.add_patch(circle)
        ax.annotate(label_txt, xy=(xm[iA], xm[iI]),
                    xytext=(xm[iA]+0.3, xm[iI]+0.2),
                    fontsize=8, color=color,
                    arrowprops=dict(arrowstyle='->', color=color, lw=0.8))

    ax.set_xlabel('A (armed actors)', fontsize=11)
    ax.set_ylabel('I (institutions)', fontsize=11)
    ax.set_title('Two paths to conflict reduction in (A, I) space\n'
                 'Circle radius ∝ kinetic capacity of c4 at midpoint', fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.2)
    ax.set_aspect('equal', adjustable='datalim')
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {save_path}')


# ===========================================================================
# FIGURE 4 — iterative_lp.png
# ===========================================================================

def fig_iterative_lp(save_path):
    """
    Iterative LP over 12 periods.
    Fixed adversarial v_C on (u1,u2,u3,u4) and decay v_N (f1,f2).
    LP maximises min_{i in peace_species} (S @ (v_C + v_P + v_N))_i
    subject to peace species not growing in conflict rows and budget cap.
    """
    x0 = np.array([8.0, 6.0, 1.0, 1.0, 5.0, 2.0])
    N_PERIODS = 12
    BUDGET    = 5.0
    COST_P    = np.array([1.0, 1.0, 2.0, 3.0])   # c1, c2, c3, c4

    k = dict(k_c1=0.1, k_c2=0.25, k_c3=0.1, k_c4=0.1,
             k_f1=0.01, k_f2=0.01)

    # Fixed adversarial process (scales with species availability)
    def v_conflict(x):
        v = np.zeros(10)
        v[U1] = 0.5 * min(x[iA], x[iRp])
        v[U2] = 0.3 * x[iA]
        v[U3] = 0.2 * min(x[iG], x[iRc])
        v[U4] = 0.1 * min(x[iA], x[iI])
        return v

    def v_decay(x):
        v = np.zeros(10)
        v[F1] = min(k['k_f1'] * x[iG]**2, x[iG]/2)
        v[F2] = min(k['k_f2'] * x[iT]**2, x[iT]/2)
        return v

    PEACE_SP = [iI, iT, iRp]    # rows we want to protect/grow
    CONF_SP  = [iG, iA, iRc]    # rows that should not grow from peace effort

    history_x  = [x0.copy()]
    history_vP = []

    x = x0.copy()
    for k_period in range(N_PERIODS):
        vC = v_conflict(x)
        vN = v_decay(x)
        vCN = vC + vN

        # Kinetic upper bounds for peace reactions at current state
        ub_c = np.array([
            k['k_c1'] * x[iT]  * x[iRp],
            k['k_c2'] * x[iI],
            k['k_c3'] * x[iI]  * x[iRp],
            k['k_c4'] * x[iI]  * x[iA],
        ])
        ub_c = np.maximum(ub_c, 0.0)

        # Peace reaction columns in STOICH: indices C1,C2,C3,C4
        S_peace = STOICH[:, [C1, C2, C3, C4]]   # (6, 4)

        # LP: maximise z (min peace species growth)
        # Variables: [v_c1, v_c2, v_c3, v_c4, z]  — 5 variables
        # Constraints:
        #   (S_peace @ v_P + S @ vCN)[i] >= z  for i in PEACE_SP
        #   ↔  -S_peace[i,:] @ v_P + z <= (S @ vCN)[i]  → negated
        #   Actually: S_peace[i,:] @ v_P - z >= -(S@vCN)[i]
        #   → -(S_peace[i,:]) @ v_P + z <= (S@vCN)[i]
        # Simplify: let base = S @ vCN
        base = STOICH @ vCN   # (6,)

        # For peace species: S_peace[i,:] @ v_P + base[i] >= z
        #   → -S_peace[i,:] @ v_P + z <= base[i]   for i in PEACE_SP
        n_ps = len(PEACE_SP)
        A_peace = np.zeros((n_ps, 5))
        b_peace = np.zeros(n_ps)
        for row_i, sp_i in enumerate(PEACE_SP):
            A_peace[row_i, :4] = -S_peace[sp_i, :]
            A_peace[row_i,  4] =  1.0
            b_peace[row_i]     = base[sp_i]

        # Budget: COST_P @ v_P <= BUDGET
        A_budget = np.zeros((1, 5))
        A_budget[0, :4] = COST_P
        b_budget = np.array([BUDGET])

        # Conflict species must not be increased by peace effort:
        #   S_peace[j,:] @ v_P <= -base[j]  for j in CONF_SP  (if base[j] < 0)
        #   only enforce when base[j] + ub contribution could be positive
        A_conf_rows = []
        b_conf_rows = []
        for sp_j in CONF_SP:
            row_c = np.zeros(5)
            row_c[:4] = S_peace[sp_j, :]
            A_conf_rows.append(row_c)
            b_conf_rows.append(-base[sp_j])

        A_ub = np.vstack([A_peace, A_budget] + A_conf_rows)
        b_ub = np.concatenate([b_peace, b_budget, b_conf_rows])

        # Bounds: v_ci in [0, ub_c[i]], z unbounded below
        bounds = [(0, float(ub_c[i])) for i in range(4)] + [(None, None)]

        # Objective: minimise -z  (maximise z)
        c_obj = np.array([0.0, 0.0, 0.0, 0.0, -1.0])

        res = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')

        if res.success:
            vP_star = res.x[:4]
        else:
            vP_star = np.zeros(4)

        # Assemble full reaction vector and update state
        v_full = vCN.copy()
        v_full[[C1, C2, C3, C4]] += vP_star

        # Contention resolution
        v_full = resolve_contention(x, v_full)
        x_new  = np.maximum(x + STOICH @ v_full, 0.0)

        history_vP.append(vP_star)
        x = x_new
        history_x.append(x.copy())

    history_x  = np.array(history_x)   # (N_PERIODS+1, 6)
    history_vP = np.array(history_vP)  # (N_PERIODS, 4)

    # ── Plot ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: trajectory in (A, I) space
    ax = axes[0]
    cmap   = plt.colormaps['RdYlBu_r']
    n_pts  = len(history_x)
    for k_p in range(n_pts - 1):
        col = cmap(k_p / max(n_pts - 2, 1))
        ax.annotate('', xy=(history_x[k_p+1, iA], history_x[k_p+1, iI]),
                    xytext=(history_x[k_p, iA],   history_x[k_p, iI]),
                    arrowprops=dict(arrowstyle='->', color=col, lw=1.8))
    sc = ax.scatter(history_x[:, iA], history_x[:, iI],
                    c=np.arange(n_pts), cmap='RdYlBu_r',
                    s=60, zorder=5, edgecolors='k', linewidths=0.5)
    plt.colorbar(sc, ax=ax, label='Period')
    ax.set_xlabel('A (armed actors)', fontsize=11)
    ax.set_ylabel('I (institutions)', fontsize=11)
    ax.set_title('LP trajectory in (A, I) space\n(early=red, late=blue)', fontsize=10)
    ax.grid(True, alpha=0.2)

    # Right: stacked bar of v_P components
    ax = axes[1]
    periods      = np.arange(N_PERIODS)
    bar_colors   = ['#2ca02c', '#bcbd22', '#1f77b4', '#ff7f0e']  # c1, c2, c3, c4
    bar_labels   = ['c1 (Rp regen)', 'c2 (trust)', 'c3 (inst.)', 'c4 (suppress)']
    bottom       = np.zeros(N_PERIODS)
    for j in range(4):
        ax.bar(periods, history_vP[:, j], bottom=bottom,
               color=bar_colors[j], label=bar_labels[j], alpha=0.85)
        bottom += history_vP[:, j]
    ax.set_xlabel('Period', fontsize=11)
    ax.set_ylabel('Optimal v_P (peace allocation)', fontsize=10)
    ax.set_title('Peace budget allocation over periods', fontsize=10)
    ax.legend(fontsize=9, loc='upper right')
    ax.set_xticks(periods)
    ax.grid(True, alpha=0.2, axis='y')

    fig.suptitle('Iterative LP — 12-period peace optimisation', fontsize=11, fontweight='bold')
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: {save_path}')


# ===========================================================================
# MAIN
# ===========================================================================

if __name__ == '__main__':
    # ── PARAMETERS ───────────────────────────────────────────────────────────
    N_STEPS = 200
    SEED    = 42
    X0      = np.array([2.0, 2.0, 4.0, 2.0, 2.0, 6.0])   # G A I T Rc Rp

    # Run single Reactive vs Reactive simulation
    print('Running simulation (Reactive vs Reactive) ...')
    X, V = run_simulation(N_STEPS, X0, seed=SEED)

    # Figure 1
    print('\nFigure 1: loop_projections.png')
    fig_loop_projections(X, V, os.path.join(FIGURES_DIR, 'loop_projections.png'))

    # Figure 2
    print('\nFigure 2: timescale_loops.png')
    fig_timescale_loops(V, os.path.join(FIGURES_DIR, 'timescale_loops.png'))

    # Figure 3
    print('\nFigure 3: two_paths.png')
    fig_two_paths(os.path.join(FIGURES_DIR, 'two_paths.png'))

    # Figure 4
    print('\nFigure 4: iterative_lp.png')
    fig_iterative_lp(os.path.join(FIGURES_DIR, 'iterative_lp.png'))

    print(f'\nDone. All figures saved to: {FIGURES_DIR}')
