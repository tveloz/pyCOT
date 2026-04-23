"""
script_basic_game.py — Discrete-time game-theoretic simulation
================================================================
Same 4+4+2 reaction network as script_basic_example.py, but with a
new kinetics: time advances in unit steps and two players (Conflict, Peace)
simultaneously allocate fractional portions of their budget stocks each step.

Reactions (indices 0–9):
    u1=0  A + Rp  => A + Rc   conflict: corruption  (A catalytic; Rp→Rc)
    u2=1  A       => A + G    conflict: propaganda   (A catalytic; G produced)
    u3=2  G + Rc  => A        conflict: recruitment  (G, Rc consumed; A gained)
    u4=3  A + I   => A        conflict: assault      (A catalytic; I destroyed)
    c1=4  T + Rp  => 2Rp      peace: Rp regeneration (T consumed; Rp net +1)
    c2=5  I       => I + T    peace: trust building  (I catalytic; T produced)
    c3=6  I + Rp  => 2I       peace: institutional   (I, Rp consumed; I net +1)
    c4=7  I + A   => I        peace: suppression     (I catalytic; A destroyed)
    f1=8  2G      =>          auto: grievance decay  (quadratic)
    f2=9  2T      =>          auto: trust decay      (quadratic)

Species (indices 0–5):  G=0  A=1  I=2  T=3  Rc=4  Rp=5

Budget pools:
    Conflict: A-pool → u1, u2, u4;  Rc-pool → u3
    Peace:    Rp-pool → c1, c3;     I-pool  → c2, c4

Secondary feasibility constraints (substrates that can go negative if uncapped):
    u1 consumes Rp;  u3 consumes G and Rc;  u4 consumes I
    c1 consumes T and Rp;  c4 consumes A

Strategies per player:
    0  Random    — Dirichlet allocation across active reactions
    1  Reactive  — weight reactions inversely by the target species level
                   (prioritise reactions that relieve the current bottleneck)
    2  Pure Cycle— fix budget ratios to balance the autocatalytic cycle;
                   enemy-destruction reaction receives residual

Investment rule (2-vector per player):
    min_threshold[i] — if stock i < threshold, invest 100 % (desperate mode)
    pct_invest[i]    — fraction to invest when stock i ≥ threshold

Simultaneous moves: both players compute their allocations from x(t) and the
combined changes are applied together.  Where both players consume the same
species, contention is resolved by proportional scaling so that x(t+1) ≥ 0.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')           # never open interactive windows; always write files
import matplotlib.pyplot as plt

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))

from pyCOT.io.functions import read_txt

# ── Stoichiometry — loaded from the reaction network file via pyCOT ───────────
# Canonical ordering (defines row/col indices for all downstream code)
SPECIES  = ['G', 'A', 'I', 'T', 'Rc', 'Rp']
RXNS     = ['u1', 'u2', 'u3', 'u4', 'c1', 'c2', 'c3', 'c4', 'f1', 'f2']

_RN_FILE  = os.path.join(_SCRIPT_DIR, 'data', 'Basic_example.txt')
_rn       = read_txt(_RN_FILE)
_sm       = _rn.stoichiometry_matrix()
_pycot_sp = list(_sm.species)
_pycot_rx = list(_sm.reactions)
_S_raw    = np.array(_sm, dtype=float)
_row_idx  = [_pycot_sp.index(s) for s in SPECIES]
_col_idx  = [_pycot_rx.index(r) for r in RXNS]
STOICH    = _S_raw[np.ix_(_row_idx, _col_idx)]   # shape (6, 10)

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_X0 = np.array([2.0, 3.0, 5.0, 3.0, 2.0, 8.0])   # G A I T Rc Rp

DEFAULT_K_AUTO = dict(k_f1=0.001, k_f2=0.001)

# Carrying-capacity ceilings — once a species reaches K_max, reactions that
# would produce it are fully suppressed.  The suppression is linear:
#   scale = max(0, 1 − x[s] / K_max[s])
# so production tapers smoothly to zero as x[s] → K_max[s].
# Order: G  A   I   T   Rc  Rp
DEFAULT_K_MAX = np.array([20.0, 15.0, 15.0, 15.0, 15.0, 20.0])

# Investment vectors: [A-pool, Rc-pool] for conflict; [Rp-pool, I-pool] for peace
DEFAULT_CONF_THRESHOLD = np.array([2.0, 1.0])
# Each row = one pool; columns = [pct_min, pct_max] — fraction sampled each step
DEFAULT_CONF_PCT       = np.array([[0.3, 0.7],   # A-pool
                                   [0.5, 0.9]])  # Rc-pool

DEFAULT_PEACE_THRESHOLD = np.array([3.0, 2.0])
DEFAULT_PEACE_PCT       = np.array([[0.3, 0.7],   # Rp-pool
                                    [0.4, 0.8]])  # I-pool

# Colour palette
_CONF_COL  = '#C0392B'
_PEACE_COL = '#2980B9'
_GREY      = '#7F8C8D'

# ── Automatic reactions ───────────────────────────────────────────────────────

def auto_decay(x, k_f1=0.001, k_f2=0.001):
    """
    Compute f1 and f2 firing amounts (fractional, not integer).
    f1 fires min(k_f1 * G^2, G/2) times → net change −2*f1_amt on G.
    f2 fires min(k_f2 * T^2, T/2) times → net change −2*f2_amt on T.
    Returns length-10 reaction vector with nonzero entries only at [8] and [9].
    """
    f1 = min(k_f1 * x[0] ** 2, x[0] / 2.0)
    f2 = min(k_f2 * x[3] ** 2, x[3] / 2.0)
    v = np.zeros(10)
    v[8] = f1
    v[9] = f2
    return v


# ── Investment budget ─────────────────────────────────────────────────────────

def _investable(stock, threshold, pct_interval, rng):
    """
    Sample a random fraction uniformly from pct_interval = [pct_min, pct_max]
    and invest that portion of the stock.
    If stock <= threshold (desperate mode) the entire stock is invested.
    """
    if stock <= threshold:
        return stock
    pct = rng.uniform(pct_interval[0], pct_interval[1])
    return pct * stock


def conflict_budgets(x, threshold, pct, rng):
    """Return (bA, bRc): stochastic investable amounts from A and Rc."""
    return (_investable(x[1], threshold[0], pct[0], rng),
            _investable(x[4], threshold[1], pct[1], rng))


def peace_budgets(x, threshold, pct, rng):
    """Return (bRp, bI): stochastic investable amounts from Rp and I."""
    return (_investable(x[5], threshold[0], pct[0], rng),
            _investable(x[2], threshold[1], pct[1], rng))


# ── Feasibility caps (secondary substrates) ───────────────────────────────────

def cap_conflict(x, v_u1, v_u2, v_u3, v_u4):
    """
    Ensure conflict reactions do not over-consume substrates.
    u1 needs Rp; u3 needs G and Rc; u4 needs I.
    A is the budget source and is catalytic in u1, u2, u4.
    """
    v_u1 = min(v_u1, x[5])               # Rp
    v_u3 = min(v_u3, x[0], x[4])         # G and Rc
    v_u4 = min(v_u4, x[2])               # I
    return v_u1, v_u2, v_u3, v_u4


def cap_peace(x, v_c1, v_c2, v_c3, v_c4):
    """
    Ensure peace reactions do not over-consume substrates.
    c1 needs T and Rp; c4 needs A.
    I is catalytic in c2, c3, c4.
    """
    v_c1 = min(v_c1, x[3], x[5])         # T and Rp
    v_c3 = min(v_c3, x[5])               # Rp (I is catalytic)
    v_c4 = min(v_c4, x[1])               # A
    return v_c1, v_c2, v_c3, v_c4


# ── Strategies ────────────────────────────────────────────────────────────────

def _dirichlet(n, rng, concentration=1.0):
    """Symmetric Dirichlet sample of length n."""
    return rng.dirichlet([concentration] * n)


# -- Strategy 0: Random --------------------------------------------------------

def alloc_conflict_random(bA, bRc, rng):
    w = _dirichlet(3, rng)
    v_u1, v_u2, v_u4 = w * bA
    return v_u1, v_u2, bRc, v_u4


def alloc_peace_random(bRp, bI, rng):
    w_rp = _dirichlet(2, rng)
    w_i  = _dirichlet(2, rng)
    v_c1, v_c3 = w_rp * bRp
    v_c2, v_c4 = w_i  * bI
    return v_c1, v_c2, v_c3, v_c4


# -- Strategy 1: Reactive ------------------------------------------------------
# Weight each reaction inversely by the level of its target/product species.
# Reactions that produce scarce species get more budget; enemy-destruction
# reactions get more budget when the enemy is strong.

def alloc_conflict_reactive(x, bA, bRc):
    G, A, I, Rc, Rp = x[0], x[1], x[2], x[4], x[5]
    eps = 1e-6
    # A-pool weights
    w_u1 = 1.0 / (Rc + eps)   # corrupt Rp→Rc when Rc scarce
    w_u2 = 1.0 / (G  + eps)   # propaganda when G scarce
    w_u4 = I   / (A  + eps)   # assault when I high relative to A
    tot  = w_u1 + w_u2 + w_u4 + eps
    v_u1 = bA * w_u1 / tot
    v_u2 = bA * w_u2 / tot
    v_u4 = bA * w_u4 / tot
    # Rc-pool: all to recruitment
    return v_u1, v_u2, bRc, v_u4


def alloc_peace_reactive(x, bRp, bI):
    A, I, T, Rp = x[1], x[2], x[3], x[5]
    eps = 1e-6
    # Rp-pool weights: c1 regenerates Rp (urgent when Rp low), c3 grows I (urgent when I low)
    w_c1 = 1.0 / (Rp + eps)
    w_c3 = 1.0 / (I  + eps)
    tot_rp = w_c1 + w_c3 + eps
    v_c1 = bRp * w_c1 / tot_rp
    v_c3 = bRp * w_c3 / tot_rp
    # I-pool weights: c2 builds T (urgent when T low), c4 suppresses A (urgent when A high)
    w_c2 = 1.0 / (T + eps)
    w_c4 = A   / (I + eps)
    tot_i = w_c2 + w_c4 + eps
    v_c2 = bI * w_c2 / tot_i
    v_c4 = bI * w_c4 / tot_i
    return v_c1, v_c2, v_c3, v_c4


# -- Strategy 2: Pure Cycle ----------------------------------------------------
# Conflict: balance u1+u2+u3 cycle; u4 receives residual A.
# Peace:    balance c1+c2+c3 cycle; c4 receives residual I.

def alloc_conflict_cycle(bA, bRc):
    # A split: 40 % u1, 40 % u2, 20 % u4 — equalises u1≈u2; u3 limited by Rc
    v_u1 = 0.40 * bA
    v_u2 = 0.40 * bA
    v_u4 = 0.20 * bA
    return v_u1, v_u2, bRc, v_u4


def alloc_peace_cycle(bRp, bI):
    # Rp split: 60 % c1 (Rp regeneration), 40 % c3 (I growth)
    v_c1 = 0.60 * bRp
    v_c3 = 0.40 * bRp
    # I split: 70 % c2 (trust building), 30 % c4 (suppression)
    v_c2 = 0.70 * bI
    v_c4 = 0.30 * bI
    return v_c1, v_c2, v_c3, v_c4


# ── Allocation dispatchers ────────────────────────────────────────────────────

def get_conflict_alloc(x, bA, bRc, strategy, rng):
    if strategy == 0:
        raw = alloc_conflict_random(bA, bRc, rng)
    elif strategy == 1:
        raw = alloc_conflict_reactive(x, bA, bRc)
    else:
        raw = alloc_conflict_cycle(bA, bRc)
    return cap_conflict(x, *raw)


def get_peace_alloc(x, bRp, bI, strategy, rng):
    if strategy == 0:
        raw = alloc_peace_random(bRp, bI, rng)
    elif strategy == 1:
        raw = alloc_peace_reactive(x, bRp, bI)
    else:
        raw = alloc_peace_cycle(bRp, bI)
    return cap_peace(x, *raw)


# ── Ceiling suppression ───────────────────────────────────────────────────────

def apply_ceilings(x, v, k_max):
    """
    For each reaction j, compute a ceiling scale factor as:
        scale_j = min over all s where STOICH[s,j] > 0  of  max(0, 1 − x[s]/K_max[s])

    A reaction that would produce a species already at its ceiling fires at
    rate 0.  A reaction that produces multiple species is gated by the most-
    saturated product.  Reactions that produce nothing (f1, f2 and consuming-
    only reactions) are unaffected.
    """
    v = v.copy()
    for j in range(10):
        scale = 1.0
        for s in range(6):
            if STOICH[s, j] > 0:            # reaction j produces species s
                scale = min(scale, max(0.0, 1.0 - x[s] / k_max[s]))
        v[j] *= scale
    return v


# ── Contention resolution ─────────────────────────────────────────────────────

def resolve_contention(x, v, max_iter=20):
    """
    Both players have already respected their own substrate constraints, but
    shared species (Rp, A in c4, etc.) may still be over-committed when both
    players' moves are combined.  Iterate: for each species that would go
    negative, scale all consuming reactions proportionally until x + S@v >= 0.
    """
    v = v.copy()
    for _ in range(max_iter):
        delta  = STOICH @ v
        x_new  = x + delta
        if np.all(x_new >= -1e-10):
            break
        for s in range(len(x)):
            if x_new[s] < 0:
                cons = STOICH[s] < 0    # which reactions consume species s
                if cons.any():
                    total_drain = -(STOICH[s] @ v)   # > 0
                    if total_drain > 1e-12:
                        scale = x[s] / total_drain
                        v[cons] *= scale
                        delta = STOICH @ v
                        x_new = x + delta
    return v


# ── Single step ───────────────────────────────────────────────────────────────

def step(x,
         conflict_strategy, peace_strategy,
         conf_threshold, conf_pct,
         peace_threshold, peace_pct,
         k_auto, k_max, rng):
    """Advance one discrete time step.  Returns (x_new, v)."""
    v = auto_decay(x, **k_auto)

    bA, bRc  = conflict_budgets(x, conf_threshold, conf_pct, rng)
    v_u1, v_u2, v_u3, v_u4 = get_conflict_alloc(x, bA, bRc, conflict_strategy, rng)
    v[0], v[1], v[2], v[3] = v_u1, v_u2, v_u3, v_u4

    bRp, bI  = peace_budgets(x, peace_threshold, peace_pct, rng)
    v_c1, v_c2, v_c3, v_c4 = get_peace_alloc(x, bRp, bI, peace_strategy, rng)
    v[4], v[5], v[6], v[7] = v_c1, v_c2, v_c3, v_c4

    v     = apply_ceilings(x, v, k_max)
    v     = resolve_contention(x, v)
    x_new = np.maximum(x + STOICH @ v, 0.0)
    return x_new, v


# ── Simulation ────────────────────────────────────────────────────────────────

def run_simulation(n_steps=300,
                   x0=None,
                   conflict_strategy=1,
                   peace_strategy=1,
                   conf_threshold=None,
                   conf_pct=None,
                   peace_threshold=None,
                   peace_pct=None,
                   k_auto=None,
                   k_max=None,
                   seed=81):
    """
    Run the game for n_steps discrete time steps.

    Parameters
    ----------
    k_max : array-like (6,), optional
        Carrying-capacity ceiling for each species [G,A,I,T,Rc,Rp].
        Production of a species is linearly suppressed as it approaches
        k_max and fully halted at k_max.  Defaults to DEFAULT_K_MAX.

    Returns
    -------
    X : ndarray (n_steps+1, 6)  — species concentrations at each step
    V : ndarray (n_steps, 10)   — reaction amounts fired each step
    """
    if x0              is None: x0              = DEFAULT_X0.copy()
    if conf_threshold  is None: conf_threshold  = DEFAULT_CONF_THRESHOLD
    if conf_pct        is None: conf_pct        = DEFAULT_CONF_PCT
    if peace_threshold is None: peace_threshold = DEFAULT_PEACE_THRESHOLD
    if peace_pct       is None: peace_pct       = DEFAULT_PEACE_PCT
    if k_auto          is None: k_auto          = DEFAULT_K_AUTO
    if k_max           is None: k_max           = DEFAULT_K_MAX

    rng = np.random.default_rng(seed)
    X   = np.zeros((n_steps + 1, 6))
    V   = np.zeros((n_steps, 10))
    X[0] = x0.copy()

    for t in range(n_steps):
        X[t + 1], V[t] = step(
            X[t],
            conflict_strategy, peace_strategy,
            conf_threshold, conf_pct,
            peace_threshold, peace_pct,
            k_auto, k_max, rng
        )
    return X, V


# ── Productive strength metrics ───────────────────────────────────────────────

def compute_strengths(v):
    """
    Productive strengths from a single reaction vector v (length 10).
    Cycle bottleneck = min over cycle reactions (weakest link sets throughput).
    """
    C_res = float(v[0])                          # conflict resource capacity: v(u1)
    P_res = float(min(v[4], v[5]))               # peace resource capacity: min(v(c1),v(c2))
    C_grp = float(min(v[0], v[1], v[2]))         # conflict cycle: min(v(u1),v(u2),v(u3))
    P_grp = float(min(v[4], v[5], v[6]))         # peace cycle:    min(v(c1),v(c2),v(c3))
    C_dst = float(v[3])                          # conflict enemy destruction: v(u4)
    P_dst = float(v[7])                          # peace enemy destruction:    v(c4)
    return dict(C_res=C_res, P_res=P_res,
                C_grp=C_grp, P_grp=P_grp,
                C_dst=C_dst, P_dst=P_dst,
                A_res=C_res - P_res,
                A_grp=C_grp - P_grp,
                A_dst=C_dst - P_dst)


def strengths_array(V):
    """Compute strengths for every step.  Returns dict of arrays (length n_steps)."""
    rows = [compute_strengths(V[t]) for t in range(len(V))]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


# ── Figure 1: Stock dynamics (4 panels) ──────────────────────────────────────

def plot_stocks(X, title_suffix='', save_path=None):
    """
    4-panel comparison of stock dynamics:
      [0,0] Social memory  — G (grievances) and T (trust)
      [0,1] Actors         — A (armed groups) and I (institutions)
      [1,0] Resources      — Rc (conflict res.) and Rp (peace res.)
      [1,1] Totals         — conflict mass (G+A+Rc) vs peace mass (I+T+Rp)
    """
    n  = X.shape[0]
    t  = np.arange(n)
    G, A, I, T, Rc, Rp = X.T

    conf_total  = G + A + Rc
    peace_total = I + T + Rp

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('Stock Dynamics' + title_suffix, fontsize=13, fontweight='bold')

    # ── Panel 0,0: Social memory ─────────────────────────────────────────────
    ax = axes[0, 0]
    ax.plot(t, G, color='#E67E22', lw=1.8, label='G  Grievances')
    ax.plot(t, T, color='#1ABC9C', lw=1.8, label='T  Trust')
    ax.set_title('Social Memory', fontsize=10)
    ax.set_xlabel('Step')
    ax.set_ylabel('Level')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)

    # ── Panel 0,1: Actors ────────────────────────────────────────────────────
    ax = axes[0, 1]
    ax.plot(t, A, color=_CONF_COL,  lw=1.8, label='A  Armed groups')
    ax.plot(t, I, color=_PEACE_COL, lw=1.8, label='I  Institutions')
    ax.set_title('Actors', fontsize=10)
    ax.set_xlabel('Step')
    ax.set_ylabel('Level')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)

    # ── Panel 1,0: Resources ─────────────────────────────────────────────────
    ax = axes[1, 0]
    ax.plot(t, Rc, color=_CONF_COL,  lw=1.8, ls='--', label='Rc  Conflict res.')
    ax.plot(t, Rp, color=_PEACE_COL, lw=1.8, ls='--', label='Rp  Peace res.')
    ax.set_title('Resources', fontsize=10)
    ax.set_xlabel('Step')
    ax.set_ylabel('Level')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)

    # ── Panel 1,1: Totals ────────────────────────────────────────────────────
    ax = axes[1, 1]
    ax.plot(t, conf_total,  color=_CONF_COL,  lw=2.0, label='Conflict  G+A+Rc')
    ax.plot(t, peace_total, color=_PEACE_COL, lw=2.0, label='Peace    I+T+Rp')
    ax.fill_between(t, conf_total, peace_total,
                    where=conf_total  >= peace_total, color=_CONF_COL,  alpha=0.12)
    ax.fill_between(t, conf_total, peace_total,
                    where=peace_total >= conf_total,  color=_PEACE_COL, alpha=0.12)
    ax.set_title('Totals', fontsize=10)
    ax.set_xlabel('Step')
    ax.set_ylabel('Mass')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.25)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)
    return fig


# ── Figure 2: Productive strengths (3 panels) ────────────────────────────────

def plot_strengths(V, title_suffix='', save_path=None):
    """
    3-panel productive strength time series.
    Each panel: conflict (red) vs peace (blue); shaded advantage region.
    """
    S    = strengths_array(V)
    t    = np.arange(len(V))

    panels = [
        ('Resource capacity',
         S['C_res'], S['P_res'],
         r'$C_{res}=v(u_1)$', r'$P_{res}=\min(v(c_1),v(c_2))$'),
        ('Group strength',
         S['C_grp'], S['P_grp'],
         r'$C_{grp}=\min(v(u_1),v(u_2),v(u_3))$',
         r'$P_{grp}=\min(v(c_1),v(c_2),v(c_3))$'),
        ('Enemy destruction',
         S['C_dst'], S['P_dst'],
         r'$C_{dst}=v(u_4)$', r'$P_{dst}=v(c_4)$'),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    fig.suptitle('Productive Strengths' + title_suffix, fontsize=13, fontweight='bold')

    for ax, (title, C, P, lc, lp) in zip(axes, panels):
        ax.plot(t, C, color=_CONF_COL,  lw=1.6, label=lc)
        ax.plot(t, P, color=_PEACE_COL, lw=1.6, label=lp)
        ax.fill_between(t, C, P, where=C >= P, color=_CONF_COL,  alpha=0.12)
        ax.fill_between(t, C, P, where=P >= C, color=_PEACE_COL, alpha=0.12)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel('Strength')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.25)
    axes[-1].set_xlabel('Step')

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)
    return fig


# ── Figure 3: Phase portrait in (A, I) space ─────────────────────────────────

def plot_phase(X, V, title_suffix='', save_path=None):
    """
    Phase portrait in (A, I) space.
    Colour encodes time (viridis: dark purple = early, yellow = late).
    Arrows show direction of motion; start and end are marked.
    """
    A_traj = X[:, 1]
    I_traj = X[:, 2]
    n      = len(V)

    cmap = plt.colormaps['viridis']
    norm = plt.Normalize(vmin=0, vmax=n - 1)

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle('Phase Portrait: Armed Groups vs Institutions' + title_suffix,
                 fontsize=12, fontweight='bold')

    for i in range(n):
        ax.plot(A_traj[i:i+2], I_traj[i:i+2],
                color=cmap(norm(i)), lw=1.3, alpha=0.75)

    arrow_gap = max(1, n // 20)
    for i in range(0, n, arrow_gap):
        dA = A_traj[i+1] - A_traj[i]
        dI = I_traj[i+1] - I_traj[i]
        if abs(dA) + abs(dI) > 1e-6:
            ax.annotate('',
                        xy=(A_traj[i] + dA, I_traj[i] + dI),
                        xytext=(A_traj[i], I_traj[i]),
                        arrowprops=dict(arrowstyle='->', color=_GREY, lw=0.9))

    ax.scatter(A_traj[0],  I_traj[0],  color='black', s=80,  zorder=5, label='Start')
    ax.scatter(A_traj[-1], I_traj[-1], color='gold',  s=100, zorder=5,
               marker='*', label='End')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02)
    cb.set_label('Time  (early → late)', fontsize=9)

    ax.set_xlabel('A  (Armed groups)',   fontsize=10)
    ax.set_ylabel('I  (Institutions)',   fontsize=10)
    ax.legend(fontsize=9, loc='upper right', bbox_to_anchor=(1.0, 1.0))
    ax.grid(True, alpha=0.25)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)
    return fig


# ── Figure 4: Stock comparison (no policy vs policy) ────────────────────────

def plot_stocks_comparison(X_A, X_B, label_A='No policy', label_B='Peace policy',
                           title_suffix='', save_path=None):
    """
    4-panel stock dynamics comparison between two runs.
    Dashed lines = Sim A (no policy); solid lines = Sim B (peace policy).
    Colour identifies the species; linestyle identifies the scenario.
    Mirrors the layout of plot_comparison() in script_basic_example.py.
    """
    n  = X_A.shape[0]
    t  = np.arange(n)
    GA, AA, IA, TA, RcA, RpA = X_A.T
    GB, AB, IB, TB, RcB, RpB = X_B.T

    pairs = [
        ('G', 'Grievances',         '#C0392B',
         'T', 'Trust',              '#2980B9',
         'G vs T  (social memory)'),
        ('A', 'Armed groups',       '#E74C3C',
         'I', 'Institutions',       '#27AE60',
         'A vs I  (actors)'),
        ('Rc', 'Conflict res.',     '#F39C12',
         'Rp', 'Peace res.',        '#1ABC9C',
         'Rc vs Rp  (resources)'),
    ]
    series_A = dict(G=GA, T=TA, A=AA, I=IA, Rc=RcA, Rp=RpA)
    series_B = dict(G=GB, T=TB, A=AB, I=IB, Rc=RcB, Rp=RpB)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_flat = axes.flatten()

    for ax, (csym, clabel, ccol, psym, plabel, pcol, title) in zip(ax_flat[:3], pairs):
        ax.plot(t, series_A[csym], color=ccol, ls='--', lw=1.5,
                label=f'{clabel} ({label_A})')
        ax.plot(t, series_B[csym], color=ccol, ls='-',  lw=2.0,
                label=f'{clabel} ({label_B})')
        ax.plot(t, series_A[psym], color=pcol, ls='--', lw=1.5,
                label=f'{plabel} ({label_A})')
        ax.plot(t, series_B[psym], color=pcol, ls='-',  lw=2.0,
                label=f'{plabel} ({label_B})')
        ax.set_title(title, fontsize=9)
        ax.set_xlabel('Step')
        ax.set_ylabel('Level')
        ax.legend(fontsize=7, ncol=2)
        ax.grid(True, alpha=0.25)

    # Panel 4 — aggregate totals
    ax4 = ax_flat[3]
    confA = GA + AA + RcA;  confB = GB + AB + RcB
    peacA = IA + TA + RpA;  peacB = IB + TB + RpB
    ax4.plot(t, confA, color='#C0392B', ls='--', lw=1.5, label=f'Σ conflict ({label_A})')
    ax4.plot(t, confB, color='#C0392B', ls='-',  lw=2.0, label=f'Σ conflict ({label_B})')
    ax4.plot(t, peacA, color='#27AE60', ls='--', lw=1.5, label=f'Σ peace ({label_A})')
    ax4.plot(t, peacB, color='#27AE60', ls='-',  lw=2.0, label=f'Σ peace ({label_B})')
    ax4.set_title('Aggregate  Σ(G+A+Rc) vs Σ(I+T+Rp)', fontsize=9)
    ax4.set_xlabel('Step')
    ax4.set_ylabel('Mass')
    ax4.legend(fontsize=7, ncol=2)
    ax4.grid(True, alpha=0.25)

    fig.suptitle(
        f'Stock Dynamics — policy comparison{title_suffix}\n'
        f'-- = {label_A}   |   — = {label_B}',
        fontsize=11, fontweight='bold'
    )
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)
    return fig


# ── Figure 5: Phase portrait comparison (no policy vs policy) ────────────────

def plot_phase_comparison(X_A, V_A, X_B, V_B,
                          label_A='No policy', label_B='Peace policy',
                          title_suffix='', save_path=None):
    """
    Phase portrait in (A, I) space showing both trajectories.
    Colour encodes time within each trajectory:
        Blues  (light→dark) = Sim A (no policy)   -- dashed
        Oranges (light→dark) = Sim B (peace policy) — solid
    This way identity (which run) and progression (when) are both visible.
    """
    n_A = len(V_A);  n_B = len(V_B)

    cmap_A = plt.colormaps['Blues']
    cmap_B = plt.colormaps['Oranges']
    # Avoid the very light end (near white) — start at 0.35
    def col_A(i): return cmap_A(0.35 + 0.65 * i / max(n_A - 1, 1))
    def col_B(i): return cmap_B(0.35 + 0.65 * i / max(n_B - 1, 1))

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle(
        f'Phase Portrait: Armed Groups vs Institutions{title_suffix}\n'
        f'Blues (--) = {label_A}   |   Oranges (—) = {label_B}',
        fontsize=11, fontweight='bold'
    )

    # Trajectory A — Blues, dashed
    Aa = X_A[:, 1];  Ia = X_A[:, 2]
    for i in range(n_A):
        ax.plot(Aa[i:i+2], Ia[i:i+2], color=col_A(i), lw=1.3, ls='--', alpha=0.75)

    # Trajectory B — Oranges, solid
    Ab = X_B[:, 1];  Ib = X_B[:, 2]
    for i in range(n_B):
        ax.plot(Ab[i:i+2], Ib[i:i+2], color=col_B(i), lw=1.8, ls='-',  alpha=0.85)

    # Arrows for A
    gap_A = max(1, n_A // 15)
    for i in range(0, n_A - 1, gap_A):
        dA, dI = Aa[i+1] - Aa[i], Ia[i+1] - Ia[i]
        if abs(dA) + abs(dI) > 1e-6:
            ax.annotate('', xy=(Aa[i]+dA, Ia[i]+dI), xytext=(Aa[i], Ia[i]),
                        arrowprops=dict(arrowstyle='->', color='steelblue', lw=0.8))

    # Arrows for B
    gap_B = max(1, n_B // 15)
    for i in range(0, n_B - 1, gap_B):
        dA, dI = Ab[i+1] - Ab[i], Ib[i+1] - Ib[i]
        if abs(dA) + abs(dI) > 1e-6:
            ax.annotate('', xy=(Ab[i]+dA, Ib[i]+dI), xytext=(Ab[i], Ib[i]),
                        arrowprops=dict(arrowstyle='->', color='darkorange', lw=0.8))

    # Start (shared) and end markers
    ax.scatter(Aa[0],  Ia[0],  color='black',      s=80,  zorder=6,
               marker='o', label='Start (shared)')
    ax.scatter(Aa[-1], Ia[-1], color='steelblue',  s=90,  zorder=6,
               marker='X', label=f'End — {label_A}')
    ax.scatter(Ab[-1], Ib[-1], color='darkorange', s=110, zorder=6,
               marker='*', label=f'End — {label_B}')

    # Two colourbars side by side
    sm_A = plt.cm.ScalarMappable(cmap=cmap_A, norm=plt.Normalize(vmin=0, vmax=1))
    sm_A.set_array([])
    sm_B = plt.cm.ScalarMappable(cmap=cmap_B, norm=plt.Normalize(vmin=0, vmax=1))
    sm_B.set_array([])
    cb_A = fig.colorbar(sm_A, ax=ax, pad=0.02, fraction=0.035)
    cb_A.set_label(f'{label_A}  (early→late)', fontsize=8)
    cb_B = fig.colorbar(sm_B, ax=ax, pad=0.12, fraction=0.035)
    cb_B.set_label(f'{label_B}  (early→late)', fontsize=8)

    ax.set_xlabel('A  (Armed groups)', fontsize=10)
    ax.set_ylabel('I  (Institutions)', fontsize=10)
    ax.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1.0, 1.0))
    ax.grid(True, alpha=0.25)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────

STRATEGY_LABELS = {0: 'Random', 1: 'Reactive', 2: 'Pure Cycle'}

SCENARIOS = [
    dict(conflict_strategy=1, peace_strategy=1, tag='reactive_vs_reactive',
         label='Conflict=Reactive  |  Peace=Reactive'),
    dict(conflict_strategy=2, peace_strategy=1, tag='cycle_vs_reactive',
         label='Conflict=Pure Cycle  |  Peace=Reactive'),
    dict(conflict_strategy=1, peace_strategy=2, tag='reactive_vs_cycle',
         label='Conflict=Reactive  |  Peace=Pure Cycle'),
    dict(conflict_strategy=2, peace_strategy=2, tag='cycle_vs_cycle',
         label='Conflict=Pure Cycle  |  Peace=Pure Cycle'),
]


if __name__ == '__main__':

    # ════════════════════════════════════════════════════════════════════════
    # PARAMETERS  — edit here; everything else runs automatically
    # ════════════════════════════════════════════════════════════════════════

    # ── Simulation length and random seed ───────────────────────────────────
    N_STEPS = 300
    SEED    = 81

    # ── Initial stocks  [G, A, I, T, Rc, Rp] ───────────────────────────────
    # Canonical order: [G, A, I, T, Rc, Rp]
    # Same x0 as script_basic_example.py for cross-kinetics comparison in paper.
    X0 = np.array([2.0, 2.0, 4.0, 2.0, 2.0, 6.0])

    # ── Carrying-capacity ceilings  [G, A, I, T, Rc, Rp] ───────────────────
    # Production of a species tapers to zero as it approaches K_MAX and stops
    # completely at K_MAX.  Set any entry to np.inf to remove that ceiling.
    K_MAX = np.array([20.0, 15.0, 15.0, 20.0, 15.0, 15.0])

    # ── Automatic decay rate constants ──────────────────────────────────────
    # f1: 2G => ;   rate = k_f1 * G^2
    # f2: 2T => ;   rate = k_f2 * T^2
    K_AUTO = dict(k_f1=0.001, k_f2=0.001)

    # ── Conflict investment rule  [A-pool, Rc-pool] ─────────────────────────
    # If stock <= threshold → invest 100 % (desperate mode)
    # If stock >  threshold → draw a random fraction from [pct_min, pct_max]
    CONF_THRESHOLD = np.array([0.5, 0.5])      # [A threshold, Rc threshold]
    CONF_PCT = np.array([[0.05, 0.15],           # A-pool:  invest 30%–70% each step
                         [0.15, 0.25]])          # Rc-pool: invest 50%–90% each step

    # ── Peace investment rule  [Rp-pool, I-pool] ────────────────────────────
    PEACE_THRESHOLD = np.array([0.5, 0.5])     # [Rp threshold, I threshold]
    PEACE_PCT = np.array([[0.3, 0.5],          # Rp-pool: invest 30%–70% each step
                          [0.3, 0.5]])         # I-pool:  invest 40%–80% each step

    # ── Scenarios to run ────────────────────────────────────────────────────
    # Strategies: 0=Random  1=Reactive  2=Pure Cycle
    # Add or remove dicts to control which combinations are executed.
    SCENARIOS = [
        dict(conflict_strategy=0, peace_strategy=0, tag='random_vs_random',
             label='Conflict=Random  |  Peace=Random'),
        dict(conflict_strategy=1, peace_strategy=1, tag='reactive_vs_reactive',
             label='Conflict=Reactive  |  Peace=Reactive'),
        dict(conflict_strategy=2, peace_strategy=2, tag='cycle_vs_reactive',
             label='Conflict=Pure Cycle  |  Peace=Reactive'),
        # dict(conflict_strategy=1, peace_strategy=2, tag='reactive_vs_cycle',
        #      label='Conflict=Reactive  |  Peace=Pure Cycle'),
        # dict(conflict_strategy=2, peace_strategy=2, tag='cycle_vs_cycle',
        #      label='Conflict=Pure Cycle  |  Peace=Pure Cycle'),
    ]

    # ── Policy comparison (mirrors no-policy vs peace-policy in det. script) ──
    # Policy lever = peace investment intervals.
    # Low fractions ≈ passive/no peacebuilding; high fractions ≈ active policy.
    # Conflict settings are identical in both runs.
    COMPARISON_STRATEGY = 1          # strategy used by both players in comparison
    PEACE_PCT_NOPOLICY  = np.array([[0.1, 0.2],   # Rp-pool: invest very little
                                    [0.01, 0.02]])  # I-pool:  invest very little
    PEACE_PCT_POLICY    = np.array([[0.1, 0.2],   # Rp-pool: active investment
                                    [0.15, 0.3]])  # I-pool:  active investment

    # ════════════════════════════════════════════════════════════════════════
    # END OF PARAMETERS
    # ════════════════════════════════════════════════════════════════════════

    OUT_DIR = os.path.normpath(
        os.path.join(_SCRIPT_DIR, '..', 'outputs', 'game_basic')
    )
    os.makedirs(OUT_DIR, exist_ok=True)

    for sc in SCENARIOS:
        tag   = sc['tag']
        label = sc['label']
        suf   = f'\n({label})'
        print(f'\nRunning: {label}')

        X, V = run_simulation(
            n_steps=N_STEPS,
            x0=X0.copy(),
            conflict_strategy=sc['conflict_strategy'],
            peace_strategy=sc['peace_strategy'],
            conf_threshold=CONF_THRESHOLD,
            conf_pct=CONF_PCT,
            peace_threshold=PEACE_THRESHOLD,
            peace_pct=PEACE_PCT,
            k_auto=K_AUTO,
            k_max=K_MAX,
            seed=SEED,
        )

        plot_stocks(X, title_suffix=suf,
                    save_path=os.path.join(OUT_DIR, f'{tag}_stocks.png'))
        plot_strengths(V, title_suffix=suf,
                       save_path=os.path.join(OUT_DIR, f'{tag}_strengths.png'))
        plot_phase(X, V, title_suffix=suf,
                   save_path=os.path.join(OUT_DIR, f'{tag}_phase.png'))

    # ── Policy comparison run ────────────────────────────────────────────────
    strat_label = STRATEGY_LABELS[COMPARISON_STRATEGY]
    print(f'\nRunning policy comparison  (strategy: {strat_label} vs {strat_label})')

    common_kwargs = dict(
        n_steps=N_STEPS,
        x0=X0.copy(),
        conflict_strategy=COMPARISON_STRATEGY,
        peace_strategy=COMPARISON_STRATEGY,
        conf_threshold=CONF_THRESHOLD,
        conf_pct=CONF_PCT,
        k_auto=K_AUTO,
        k_max=K_MAX,
        seed=SEED,
    )

    X_np, V_np = run_simulation(peace_threshold=PEACE_THRESHOLD,
                                peace_pct=PEACE_PCT_NOPOLICY, **common_kwargs)
    X_pp, V_pp = run_simulation(peace_threshold=PEACE_THRESHOLD,
                                peace_pct=PEACE_PCT_POLICY,   **common_kwargs)

    comp_suf = f'\n(Strategy: {strat_label} | -- = No policy, — = Peace policy)'
    plot_stocks_comparison(
        X_np, X_pp,
        label_A='No policy', label_B='Peace policy',
        title_suffix=comp_suf,
        save_path=os.path.join(OUT_DIR, 'comparison_stocks.png'),
    )
    plot_phase_comparison(
        X_np, V_np, X_pp, V_pp,
        label_A='No policy', label_B='Peace policy',
        title_suffix=comp_suf,
        save_path=os.path.join(OUT_DIR, 'comparison_phase.png'),
    )

    print('\nDone.')
