#!/usr/bin/env python3
"""
script_section_3_4.py
=====================
Section 3.4 — A Kinetic View.  Produces Figures 3 and 4 for the paper.

Figure 3: Deterministic Michaelis-Menten kinetics.
    Two scenarios from x0=(G=2,A=2,I=4,T=2,Rc=2,Rp=6):
      A — no policy  (alpha_inst=0, gamma_sec=0)
      B — active peacebuilding  (alpha_inst=0.4, gamma_sec=0.5)
    Output: 4-panel time-series + (A,I) phase portrait.

Figure 4: Stochastic Poisson tau-leaping kinetics.
    Propensity:  lambda_j = a_j * kappa_j * pi_in_j(x)
    Firings:     n_j ~ Poisson(lambda_j * dt), capped at floor(pi_in_j(x))
    Two scenarios, same x0 and seed:
      A — no policy   (peace allocation = NOPOLICY_SCALE * strategy allocation)
      B — peace strategy (full allocation, emphasising c3/c4)
    Conflict allocation a_C is identical in both scenarios.

Resource-burn diagnostic (Figure 4c):
    Higher a_j raises lambda_j and thus expected firings, but also depletes
    shared resources faster.  The enabling degree pi_in_j(x), which falls as
    reactants deplete, provides the natural self-limiting feedback.
    The Rp panel shows the burn-vs-replenishment tradeoff for c3 (burns Rp
    to grow I) vs c1 (replenishes Rp via trust once I is large enough).
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Path setup ────────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))

from pyCOT.io.functions     import read_txt
from pyCOT.simulations.ode  import simulation
from pyCOT.simulations.core import build_reaction_dict

# ── Output directory ──────────────────────────────────────────────────────────
_OUT_DIR = os.path.normpath(
    os.path.join(_SCRIPT_DIR, '..', 'outputs', 'section_3_4'))
os.makedirs(_OUT_DIR, exist_ok=True)

# ── Network + stoichiometry (canonical order) ─────────────────────────────────
_RN_FILE  = os.path.join(_SCRIPT_DIR, 'data', 'Basic_example.txt')
_rn       = read_txt(_RN_FILE)
_sm       = _rn.stoichiometry_matrix()

SPECIES   = ['G', 'A', 'I', 'T', 'Rc', 'Rp']   # row order
REACTIONS = ['u1','u2','u3','u4','c1','c2','c3','c4','f1','f2']  # col order

_pycot_sp  = list(_sm.species)
_pycot_rx  = list(_sm.reactions)
_row_idx   = [_pycot_sp.index(s) for s in SPECIES]
_col_idx   = [_pycot_rx.index(r) for r in REACTIONS]
STOICH     = np.array(_sm, dtype=float)[np.ix_(_row_idx, _col_idx)]  # (6,10)

# Species index shortcuts
iG, iA, iI, iT, iRc, iRp = 0, 1, 2, 3, 4, 5

# Initial state as stated in the paper
X0 = np.array([2.0, 2.0, 4.0, 2.0, 2.0, 6.0])   # [G, A, I, T, Rc, Rp]

# ── Input stoichiometry — reactant (consumption) coefficients ─────────────────
# Used exclusively for computing enabling degrees pi_in_j(x).
# Entry [s, j] = how many units of species s reaction j consumes per firing.
# Rows: G, A, I, T, Rc, Rp
# Cols: u1, u2, u3, u4, c1, c2, c3, c4, f1, f2
INPUT_STOICH = np.array([
#   u1  u2  u3  u4  c1  c2  c3  c4  f1  f2
    [ 0,  0,  1,  0,  0,  0,  0,  0,  2,  0],  # G
    [ 1,  1,  0,  1,  0,  0,  0,  1,  0,  0],  # A
    [ 0,  0,  0,  1,  0,  1,  1,  1,  0,  0],  # I
    [ 0,  0,  0,  0,  1,  0,  0,  0,  0,  2],  # T
    [ 0,  0,  1,  0,  0,  0,  0,  0,  0,  0],  # Rc
    [ 1,  0,  0,  0,  1,  0,  1,  0,  0,  0],  # Rp
], dtype=float)

# ── Colour palette ────────────────────────────────────────────────────────────
_CRED  = '#C0392B'   # conflict / armed groups
_CBLUE = '#2980B9'   # peace / trust
_CGRN  = '#27AE60'   # institutions
_CORNG = '#F39C12'   # Rc (conflict resources)
_CTEAL = '#1ABC9C'   # Rp / Trust


# ===========================================================================
# PART A — DETERMINISTIC MICHAELIS-MENTEN KINETICS  (Figure 3)
# ===========================================================================

# ── Kinetic rate functions (pyCOT calling convention) ─────────────────────────
def kin_u1(_, c, idx, sv):
    k, Km = sv
    return k * c[idx['A']] * c[idx['Rp']] / (c[idx['Rp']] + Km)

def kin_u2(_, c, idx, sv):
    k, Km = sv
    return k * c[idx['A']] / (c[idx['A']] + Km)

def kin_u3(_, c, idx, sv):
    k, Km = sv
    return k * c[idx['G']] * c[idx['Rc']] / (c[idx['G']] + Km)

def kin_u4(_, c, idx, sv):
    k, Km = sv
    return k * c[idx['A']] * c[idx['I']] / (c[idx['I']] + Km)

def kin_c1(_, c, idx, sv):
    k, Km = sv
    return k * c[idx['T']] * c[idx['Rp']] / (c[idx['Rp']] + Km)

def kin_c2(_, c, idx, sv):
    k, Km = sv
    return k * c[idx['I']] / (c[idx['I']] + Km)

def kin_c3(_, c, idx, sv):
    k, Km, alpha = sv
    return alpha * k * c[idx['I']] * c[idx['Rp']] / (c[idx['I']] + Km)

def kin_c4(_, c, idx, sv):
    k, Km, gamma = sv
    return gamma * k * c[idx['A']] * c[idx['I']] / (c[idx['A']] + Km)

def kin_f1(_, c, idx, sv):
    return sv[0] * c[idx['G']] ** 2

def kin_f2(_, c, idx, sv):
    return sv[0] * c[idx['T']] ** 2

CUSTOM_LAWS = {
    'kin_u1': kin_u1, 'kin_u2': kin_u2,
    'kin_u3': kin_u3, 'kin_u4': kin_u4,
    'kin_c1': kin_c1, 'kin_c2': kin_c2,
    'kin_c3': kin_c3, 'kin_c4': kin_c4,
    'kin_f1': kin_f1, 'kin_f2': kin_f2,
}
RATE_NAMES = [
    'kin_u1','kin_u2','kin_u3','kin_u4',
    'kin_c1','kin_c2','kin_c3','kin_c4',
    'kin_f1','kin_f2',
]

# Default kinetic parameters (from script_basic_example.py)
_DEFAULT_SV = [
    [0.05,  2.0],          # u1  [k, Km]
    [0.2,   2.0],          # u2  [k, Km]
    [0.05,  3.0],          # u3  [k, Km]
    [0.1,   2.0],          # u4  [k, Km]
    [0.1,   2.0],          # c1  [k, Km]
    [0.25,  2.0],          # c2  [k, Km]
    [0.1,   5.0, 0.4],     # c3  [k, Km, alpha_inst]
    [0.1,   2.0, 0.5],     # c4  [k, Km, gamma_sec]
    [0.01],                # f1  [k]
    [0.01],                # f2  [k]
]

def _make_sv(alpha_inst=0.4, gamma_sec=0.5):
    sv = [list(p) for p in _DEFAULT_SV]
    sv[6][-1] = alpha_inst
    sv[7][-1] = gamma_sec
    return sv


def run_ode(sv, t_end=300.0, n_steps=1500):
    """Run pyCOT ODE simulation; return DataFrame with Time + SPECIES columns."""
    pycot_order = list(_sm.species)
    x0_dict  = {sp: float(v) for sp, v in zip(SPECIES, X0)}
    x0_pycot = [x0_dict[sp] for sp in pycot_order]
    ts_df, _ = simulation(
        _rn,
        rate=RATE_NAMES, spec_vector=sv, x0=x0_pycot,
        t_span=(0.0, t_end), n_steps=n_steps,
        additional_laws=CUSTOM_LAWS, method='LSODA', verbose=False,
    )
    return ts_df


def plot_fig3_timeseries(ts_A, ts_B, save_path=None):
    """4-panel time-series comparison (no policy vs peace policy)."""
    time = ts_A['Time'].values
    pairs = [
        ('G', 'Grievances',   _CRED,  'T',  'Trust',        _CBLUE,
         'G vs T  (social memory)'),
        ('A', 'Armed groups', _CRED,  'I',  'Institutions',  _CGRN,
         'A vs I  (actors)'),
        ('Rc','Conflict res.',_CORNG, 'Rp', 'Peace res.',    _CTEAL,
         'Rc vs Rp  (resources)'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_flat   = axes.flatten()

    for ax, (cs, cl, cc, ps, pl, pc, title) in zip(ax_flat[:3], pairs):
        ax.plot(time, ts_A[cs], color=cc, ls='--', lw=1.5, label=f'{cl} (A)')
        ax.plot(time, ts_B[cs], color=cc, ls='-',  lw=2.0, label=f'{cl} (B)')
        ax.plot(time, ts_A[ps], color=pc, ls='--', lw=1.5, label=f'{pl} (A)')
        ax.plot(time, ts_B[ps], color=pc, ls='-',  lw=2.0, label=f'{pl} (B)')
        ax.set_title(title, fontsize=9)
        ax.set_xlabel('Time'); ax.set_ylabel('Concentration')
        ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)

    # Panel 4 — aggregate totals
    ax4 = ax_flat[3]
    for ts, ls, lw, lbl in [(ts_A, '--', 1.5, 'A'), (ts_B, '-', 2.0, 'B')]:
        conf = ts['G'] + ts['A'] + ts['Rc']
        peac = ts['I'] + ts['T'] + ts['Rp']
        ax4.plot(time, conf, color=_CRED, ls=ls, lw=lw,
                 label=f'Σ conflict ({lbl})')
        ax4.plot(time, peac, color=_CGRN, ls=ls, lw=lw,
                 label=f'Σ peace ({lbl})')
    ax4.set_title('Aggregate  Σ(G+A+Rc) vs Σ(I+T+Rp)', fontsize=9)
    ax4.set_xlabel('Time'); ax4.set_ylabel('Total')
    ax4.legend(fontsize=7, ncol=2); ax4.grid(True, alpha=0.3)

    fig.suptitle(
        'Conflict-Peace Dynamics — pairwise comparison\n'
        'A (--): α_inst=0, γ_sec=0  |  B (—): α_inst=0.4, γ_sec=0.5',
        fontsize=10, y=1.01)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)


def plot_fig3_phase(ts_A, ts_B, save_path=None):
    """Phase portrait (A, I) for deterministic MM — Blues vs Oranges."""
    A_a, I_a = ts_A['A'].values, ts_A['I'].values
    A_b, I_b = ts_B['A'].values, ts_B['I'].values
    n_a, n_b = len(A_a), len(A_b)

    cmap_A = plt.colormaps['Blues']
    cmap_B = plt.colormaps['Oranges']

    fig, ax = plt.subplots(figsize=(8, 6))
    for i in range(n_a - 1):
        c = cmap_A(0.35 + 0.65 * i / max(n_a - 1, 1))
        ax.plot(A_a[i:i+2], I_a[i:i+2], color=c, lw=1.2, ls='--', alpha=0.75)
    for i in range(n_b - 1):
        c = cmap_B(0.35 + 0.65 * i / max(n_b - 1, 1))
        ax.plot(A_b[i:i+2], I_b[i:i+2], color=c, lw=1.8, ls='-',  alpha=0.85)

    ax.scatter(A_a[0],  I_a[0],  color='black',     s=80,  zorder=6,
               marker='o', label='Start (shared)')
    ax.scatter(A_a[-1], I_a[-1], color='steelblue', s=90,  zorder=6,
               marker='X', label='End — no policy')
    ax.scatter(A_b[-1], I_b[-1], color='darkorange', s=110, zorder=6,
               marker='*', label='End — peace policy')

    sm_A = plt.cm.ScalarMappable(cmap=cmap_A, norm=plt.Normalize(0, 1))
    sm_B = plt.cm.ScalarMappable(cmap=cmap_B, norm=plt.Normalize(0, 1))
    sm_A.set_array([]); sm_B.set_array([])
    fig.colorbar(sm_A, ax=ax, pad=0.02,  fraction=0.035).set_label(
        'No policy (early→late)', fontsize=8)
    fig.colorbar(sm_B, ax=ax, pad=0.12, fraction=0.035).set_label(
        'Peace policy (early→late)', fontsize=8)

    ax.set_xlabel('A  (Armed groups)', fontsize=10)
    ax.set_ylabel('I  (Institutions)', fontsize=10)
    ax.set_title(
        'Phase Portrait: Armed Groups vs Institutions\n'
        'Blues (--): α_inst=0, γ_sec=0  (no policy)   |   '
        'Oranges (—): α_inst=0.4, γ_sec=0.5  (peace policy)',
        fontsize=9, fontweight='bold')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.25)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)


# ===========================================================================
# PART B — STOCHASTIC POISSON TAU-LEAPING KINETICS  (Figure 4)
# ===========================================================================

# ── Kinetic effectiveness coefficients kappa_j ────────────────────────────────
# Order: [u1,   u2,   u3,   u4,   c1,   c2,   c3,   c4,   f1,   f2]
KAPPA = np.array([0.10, 0.15, 0.08, 0.12,
                  0.15, 0.20, 0.12, 0.15,
                  0.05, 0.05])

# ── Allocation vectors ────────────────────────────────────────────────────────
# Conflict (fixed across both scenarios) — simplex over {u1,u2,u3,u4}
A_CONFLICT = np.array([0.30, 0.25, 0.25, 0.20])   # sum = 1.0

# Peace strategy — simplex over {c1,c2,c3,c4}; c3 and c4 weighted more
A_PEACE_STRATEGY  = np.array([0.30, 0.25, 0.30, 0.15])   # sum = 1.0

# No-policy: scaled-down version — all peace reactions fire but weakly
NOPOLICY_SCALE    = 0.30
A_PEACE_NOPOLICY  = NOPOLICY_SCALE * A_PEACE_STRATEGY     # sum = 0.30


def _full_allocation(a_peace):
    """
    Build full 10-element allocation vector:
      [u1, u2, u3, u4, c1, c2, c3, c4, f1, f2]
    Auto-decay reactions f1, f2 are uncontrolled: allocation = 1.0
    (they fire at kappa * pi_in, unattenuated by any player budget).
    """
    a = np.zeros(10)
    a[0:4] = A_CONFLICT      # conflict player
    a[4:8] = a_peace          # peace player
    a[8:10] = 1.0             # autonomous decay
    return a

A_NOPOLICY = _full_allocation(A_PEACE_NOPOLICY)
A_STRATEGY = _full_allocation(A_PEACE_STRATEGY)


def enabling_degrees(x):
    """
    Enabling degree vector pi_in_j(x) for all 10 reactions.
    pi_in_j(x) = min over s with INPUT_STOICH[s,j] > 0  of  x[s] / INPUT_STOICH[s,j]
    Returns shape (10,); reactions with no reactants get inf.
    """
    pi = np.full(10, np.inf)
    for j in range(10):
        col    = INPUT_STOICH[:, j]
        active = col > 0
        if active.any():
            pi[j] = np.min(x[active] / col[active])
    return np.clip(pi, 0.0, None)


def _resolve_contention(x, n, max_iter=20):
    """
    Proportionally scale down any reactions that would drive a species negative.
    Iterates until x + STOICH @ n >= 0 (within tolerance) or max_iter reached.
    """
    n = n.copy()
    for _ in range(max_iter):
        delta = STOICH @ n
        if np.all(x + delta >= -1e-10):
            break
        for s in range(6):
            if x[s] + delta[s] < -1e-10:
                consuming = STOICH[s] < 0
                if consuming.any():
                    total_drain = -(STOICH[s] @ n)
                    if total_drain > 1e-12:
                        scale = max(0.0, x[s] / total_drain)
                        n[consuming] *= scale
                        delta = STOICH @ n
    return n


def step_poisson(x, a, kappa, dt, rng):
    """
    One tau-leaping step.

    1. Compute enabling degrees pi_in_j(x).
    2. Propensity: lambda_j = a_j * kappa_j * pi_in_j(x).
    3. Sample: n_j ~ Poisson(lambda_j * dt), capped at floor(pi_in_j(x)).
    4. Resolve shared-resource contention.
    5. Update state: x_new = max(x + S @ n, 0).

    The cap at floor(pi_in) ensures no single reaction fires more times than
    the current reactant stock can support.  Contention resolution handles
    the residual case where two reactions sharing a reactant collectively
    over-commit after independent capping.

    Parameters
    ----------
    x     : state vector (6,)
    a     : full allocation vector (10,)
    kappa : effectiveness coefficients (10,)
    dt    : time increment
    rng   : numpy Generator

    Returns
    -------
    x_new : (6,)
    n     : actual firing counts (10,)
    pi    : enabling degrees used (10,)
    lam   : propensities used (10,)
    """
    pi  = enabling_degrees(x)                   # (10,)
    lam = a * kappa * pi                        # (10,)  propensities

    # Poisson draw, capped at floor(enabling degree)
    n_raw = rng.poisson(lam * dt).astype(float)
    n     = np.minimum(n_raw, np.floor(pi))

    n     = _resolve_contention(x, n)
    x_new = np.maximum(x + STOICH @ n, 0.0)
    return x_new, n, pi, lam


def run_poisson(a, n_steps=300, dt=1.0, seed=42):
    """
    Run the stochastic tau-leaping simulation for n_steps steps.

    Returns
    -------
    X   : (n_steps+1, 6)  species trajectories
    N   : (n_steps,   10) firing counts per step
    PI  : (n_steps,   10) enabling degrees per step
    LAM : (n_steps,   10) propensities per step
    """
    rng = np.random.default_rng(seed)
    X   = np.zeros((n_steps + 1, 6))
    N   = np.zeros((n_steps,    10))
    PI  = np.zeros((n_steps,    10))
    LAM = np.zeros((n_steps,    10))
    X[0] = X0.copy()

    for t in range(n_steps):
        X[t+1], N[t], PI[t], LAM[t] = step_poisson(X[t], a, KAPPA, dt, rng)
    return X, N, PI, LAM



def plot_fig4_stocks(X_A, X_B,
                     label_A='No policy', label_B='Peace strategy',
                     save_path=None):
    """
    Figure 4 top: 4-panel stock dynamics comparing two scenarios.
    Dashed = no policy; solid = peace strategy.
    """
    n = X_A.shape[0]
    t = np.arange(n)

    GA, AA, IA, TA, RcA, RpA = X_A.T
    GB, AB, IB, TB, RcB, RpB = X_B.T

    pairs = [
        ('G', GA, 'Grievances', _CRED,
         'T', TA, 'Trust',       _CBLUE,
         GB, TB,  'G vs T  (social memory)'),
        ('A', AA, 'Armed groups', _CRED,
         'I', IA, 'Institutions', _CGRN,
         AB, IB,  'A vs I  (actors)'),
        ('Rc', RcA, 'Conflict res.', _CORNG,
         'Rp', RpA, 'Peace res.',    _CTEAL,
         RcB, RpB,  'Rc vs Rp  (resources)'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_flat   = axes.flatten()

    for ax, (cs, dA, cl, cc, ps, dP, pl, pc, dAb, dPb, title) in \
            zip(ax_flat[:3], pairs):
        ax.plot(t, dA,  color=cc, ls='--', lw=1.5, label=f'{cl} ({label_A})')
        ax.plot(t, dAb, color=cc, ls='-',  lw=2.0, label=f'{cl} ({label_B})')
        ax.plot(t, dP,  color=pc, ls='--', lw=1.5, label=f'{pl} ({label_A})')
        ax.plot(t, dPb, color=pc, ls='-',  lw=2.0, label=f'{pl} ({label_B})')
        ax.set_title(title, fontsize=9)
        ax.set_xlabel('Step'); ax.set_ylabel('Level')
        ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.25)

    # Panel 4 — aggregate totals
    ax4  = ax_flat[3]
    confA = GA + AA + RcA;  peacA = IA + TA + RpA
    confB = GB + AB + RcB;  peacB = IB + TB + RpB
    ax4.plot(t, confA, color=_CRED, ls='--', lw=1.5,
             label=f'Σ conflict ({label_A})')
    ax4.plot(t, confB, color=_CRED, ls='-',  lw=2.0,
             label=f'Σ conflict ({label_B})')
    ax4.plot(t, peacA, color=_CGRN, ls='--', lw=1.5,
             label=f'Σ peace ({label_A})')
    ax4.plot(t, peacB, color=_CGRN, ls='-',  lw=2.0,
             label=f'Σ peace ({label_B})')
    ax4.fill_between(t, confB, peacB,
                     where=(confB >= peacB), color=_CRED,  alpha=0.10)
    ax4.fill_between(t, confB, peacB,
                     where=(peacB >= confB), color=_CGRN, alpha=0.10)
    ax4.set_title('Aggregate  Σ(G+A+Rc) vs Σ(I+T+Rp)', fontsize=9)
    ax4.set_xlabel('Step'); ax4.set_ylabel('Mass')
    ax4.legend(fontsize=7, ncol=2); ax4.grid(True, alpha=0.25)

    fig.suptitle(
        f'Stock Dynamics — policy comparison\n'
        f'Stochastic game-theoretic kinetics\n'
        f'-- = {label_A}   |   — = {label_B}',
        fontsize=11, fontweight='bold')
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)


def plot_fig4_phase(X_A, X_B,
                    label_A='No policy', label_B='Peace strategy',
                    save_path=None):
    """
    Figure 4 bottom: phase portrait (A, I) — same visual style as Fig 3 phase.
    Blues gradient (dashed)  = no policy trajectory (early→late).
    Oranges gradient (solid) = peace strategy trajectory (early→late).
    """
    A_a, I_a = X_A[:, iA], X_A[:, iI]
    A_b, I_b = X_B[:, iA], X_B[:, iI]
    n_a, n_b = len(A_a), len(A_b)

    cmap_A = plt.colormaps['Blues']
    cmap_B = plt.colormaps['Oranges']

    fig, ax = plt.subplots(figsize=(8, 6))

    for i in range(n_a - 1):
        c = cmap_A(0.35 + 0.65 * i / max(n_a - 1, 1))
        ax.plot(A_a[i:i+2], I_a[i:i+2], color=c, lw=1.2, ls='--', alpha=0.75)
    for i in range(n_b - 1):
        c = cmap_B(0.35 + 0.65 * i / max(n_b - 1, 1))
        ax.plot(A_b[i:i+2], I_b[i:i+2], color=c, lw=1.8, ls='-',  alpha=0.85)

    ax.scatter(A_a[0],  I_a[0],  color='black',      s=80,  zorder=6,
               marker='o', label='Start (shared)')
    ax.scatter(A_a[-1], I_a[-1], color='steelblue',  s=90,  zorder=6,
               marker='X', label=f'End — {label_A}')
    ax.scatter(A_b[-1], I_b[-1], color='darkorange',  s=110, zorder=6,
               marker='*', label=f'End — {label_B}')

    sm_A = plt.cm.ScalarMappable(cmap=cmap_A, norm=plt.Normalize(0, 1))
    sm_B = plt.cm.ScalarMappable(cmap=cmap_B, norm=plt.Normalize(0, 1))
    sm_A.set_array([]); sm_B.set_array([])
    fig.colorbar(sm_A, ax=ax, pad=0.02,  fraction=0.035).set_label(
        f'{label_A}  (early→late)', fontsize=8)
    fig.colorbar(sm_B, ax=ax, pad=0.12, fraction=0.035).set_label(
        f'{label_B}  (early→late)', fontsize=8)

    ax.set_xlabel('A  (Armed groups)', fontsize=10)
    ax.set_ylabel('I  (Institutions)', fontsize=10)
    ax.set_title(
        'Phase Portrait: Armed Groups vs Institutions\n'
        f'Blues (--): {label_A}   |   Oranges (—): {label_B}',
        fontsize=9, fontweight='bold')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.25)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)


def plot_resource_burn(LAM_A, LAM_B, PI_A, PI_B,
                       label_A='No policy', label_B='Peace strategy',
                       save_path=None):
    """
    Diagnostic — Figure 4c:  resource burn vs self-limiting feedback.

    Left panel:  Enabling degrees of the three Rp-consuming reactions over time.
        pi_in for u1 (corruption), c1 (trust→Rp), c3 (institutional growth).
        As Rp depletes these drop, throttling all three reactions automatically.

    Right panel: Propensities of c3 (burns Rp to grow I) vs c1 (replenishes Rp
        via trust).  Under peace strategy: c3 > c1 early (net Rp burn) but
        as I and T grow, c1 catches up, creating the eventual steady state.
        Under no-policy: both c3 and c1 are much weaker — Rp burn is slower
        but so is institutional growth.

    This illustrates the core resource-burn tradeoff:
        Higher a_j (peace strategy) → larger lambda_j → faster Rp depletion
        via c3, but also faster I growth → faster T production via c2 → faster
        Rp replenishment via c1.  The equilibrium depends on whether the I→T→Rp
        loop (c2+c1) closes fast enough to offset c3 consumption.
    """
    t_A = np.arange(LAM_A.shape[0])
    t_B = np.arange(LAM_B.shape[0])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: enabling degrees of Rp-consuming reactions
    ax = axes[0]
    # pi_in for u1 = PI[:,0], c1 = PI[:,4], c3 = PI[:,6]
    # (all three require Rp, so their enabling degrees reflect Rp availability)
    for PI, t, ls, lbl in [(PI_A, t_A, '--', label_A), (PI_B, t_B, '-', label_B)]:
        ax.plot(t, PI[:, 0], color=_CRED,  ls=ls, lw=1.5, label=f'π(u1) corruption — {lbl}')
        ax.plot(t, PI[:, 4], color=_CTEAL, ls=ls, lw=1.5, label=f'π(c1) trust→Rp — {lbl}')
        ax.plot(t, PI[:, 6], color=_CGRN,  ls=ls, lw=1.5, label=f'π(c3) inst. growth — {lbl}')
    ax.set_xlabel('Step')
    ax.set_ylabel('Enabling degree  π_in_j(x)')
    ax.set_title('Enabling degrees of Rp-consuming reactions\n'
                 '(drop as Rp depletes → automatic self-throttling)', fontsize=9)
    ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.25)

    # Right: c3 (burn) vs c1 (replenish) propensities
    ax = axes[1]
    for LAM, t, ls, lbl in [(LAM_A, t_A, '--', label_A), (LAM_B, t_B, '-', label_B)]:
        ax.plot(t, LAM[:, 6], color=_CGRN,  ls=ls, lw=2.0,
                label=f'λ(c3) inst. growth — burns Rp — {lbl}')
        ax.plot(t, LAM[:, 4], color=_CTEAL, ls=ls, lw=2.0,
                label=f'λ(c1) trust→Rp — replenishes — {lbl}')
    ax.set_xlabel('Step')
    ax.set_ylabel('Propensity  λ_j')
    ax.set_title('Resource burn vs replenishment tradeoff\n'
                 'c3 burns Rp to grow I;  c1 replenishes Rp via trust\n'
                 '(c1 catches up once I and T are sufficiently large)', fontsize=9)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.25)

    fig.suptitle(
        'Poisson kinetics — resource burn and self-limiting feedback\n'
        'Higher allocation a_j raises λ_j and resource consumption;\n'
        'enabling degrees fall as reactants deplete, closing the feedback loop.',
        fontsize=10, fontweight='bold')
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'  Saved: {save_path}')
    plt.close(fig)


# ===========================================================================
# MAIN
# ===========================================================================

if __name__ == '__main__':
    print('=' * 65)
    print('Section 3.4 — A Kinetic View')
    print('=' * 65)
    print(f'x0 = (G={X0[iG]}, A={X0[iA]}, I={X0[iI]}, '
          f'T={X0[iT]}, Rc={X0[iRc]}, Rp={X0[iRp]})')

    # ── FIGURE 3: Deterministic MM ─────────────────────────────────────────
    print('\n--- PART A: Deterministic MM kinetics (Figure 3) ---')

    sv_nopolicy = _make_sv(alpha_inst=0.0, gamma_sec=0.0)
    sv_peace    = _make_sv(alpha_inst=0.4, gamma_sec=0.5)

    print('  Running no-policy ODE (T=300, 1500 steps) ...')
    ts_A = run_ode(sv_nopolicy, t_end=300.0, n_steps=1500)
    print('  Running peace-policy ODE (T=300, 1500 steps) ...')
    ts_B = run_ode(sv_peace,    t_end=300.0, n_steps=1500)

    plot_fig3_timeseries(
        ts_A, ts_B,
        save_path=os.path.join(_OUT_DIR, 'fig3_deterministic_timeseries.png'))
    plot_fig3_phase(
        ts_A, ts_B,
        save_path=os.path.join(_OUT_DIR, 'fig3_deterministic_phase.png'))
    print('  Figure 3 done.')

    # ── FIGURE 4: Stochastic Poisson ──────────────────────────────────────
    print('\n--- PART B: Stochastic Poisson tau-leaping kinetics (Figure 4) ---')
    print(f'  Conflict allocation (fixed):   {dict(zip(["u1","u2","u3","u4"], A_CONFLICT))}')
    print(f'  Peace strategy allocation:     {dict(zip(["c1","c2","c3","c4"], A_PEACE_STRATEGY))}')
    print(f'  No-policy allocation (×{NOPOLICY_SCALE}):  '
          f'{dict(zip(["c1","c2","c3","c4"], A_PEACE_NOPOLICY.round(3)))}')
    print(f'  kappa = {dict(zip(REACTIONS, KAPPA))}')

    print('  Running no-policy Poisson simulation (N=300 steps) ...')
    X_A, N_A, PI_A, LAM_A = run_poisson(A_NOPOLICY, n_steps=300, seed=42)
    print('  Running peace-strategy Poisson simulation (N=300 steps) ...')
    X_B, N_B, PI_B, LAM_B = run_poisson(A_STRATEGY,  n_steps=300, seed=42)

    plot_fig4_stocks(
        X_A, X_B,
        save_path=os.path.join(_OUT_DIR, 'fig4_stochastic_stocks.png'))
    plot_fig4_phase(
        X_A, X_B,
        save_path=os.path.join(_OUT_DIR, 'fig4_stochastic_phase.png'))
    plot_resource_burn(
        LAM_A, LAM_B, PI_A, PI_B,
        save_path=os.path.join(_OUT_DIR, 'fig4_resource_burn_diagnostic.png'))
    print('  Figure 4 done.')

    print(f'\nAll outputs saved to: {_OUT_DIR}')
