#!/usr/bin/env python3
"""
XCEPT Article 1 — Basic Example Script
=======================================
Conflict-Peace Minimal Reaction Network

This script demonstrates the full workflow for the Basic_example.txt model:
  Part 2  — Load the network from file and inspect stoichiometric structure
  Part 3  — Attach custom kinetic functions (pyCOT interface)
  Part 4  — Compute process vectors and state changes
  Part 5  — Run time-course simulations (baseline vs peacebuilding policy)
  Part 6  — Interactive process-vector testing utility

All simulation work goes through pyCOT's simulation() engine and kinetics
interface. Custom kinetic functions follow the pyCOT calling convention:
    f(_, concentrations, species_idx, spec_vector) -> float

Run directly:
    python script_basic_example.py
"""

import os
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Path setup — allow import of pyCOT from the repo src/ folder.
# Script lives at  .../pyCOT/projects/XCEPT_Article_1/scripts/
# Three dirname() calls reach the pyCOT root.
# ---------------------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PYCOT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))
sys.path.insert(0, os.path.join(_PYCOT_ROOT, 'src'))

# --- pyCOT imports ---
from pyCOT.io.functions import read_txt                        # network loader
from pyCOT.simulations.ode import simulation                   # ODE integrator
from pyCOT.simulations.core import build_reaction_dict         # reaction→species map


# ===========================================================================
# PART 2 — Load reaction network and inspect stoichiometric structure
# ===========================================================================

def load_network(filepath):
    """Load a reaction network from a .txt file via pyCOT's read_txt parser."""
    return read_txt(filepath)


def get_stoichiometry(rn, canonical_species, canonical_reactions):
    """
    Extract pyCOT's stoichiometric matrix and reindex rows/columns to the
    caller-specified canonical orderings.

    pyCOT sorts species and reactions by their graph-node insertion index,
    so the internal order depends on the order reactions appear in the file.
    This function maps any internal ordering to the canonical one.

    Returns
    -------
    S          : np.ndarray (n_species × n_reactions)
    sp_names   : list[str]  — row labels  = canonical_species
    rx_names   : list[str]  — col labels  = canonical_reactions
    """
    sm = rn.stoichiometry_matrix()          # StoichiometryMatrix object
    pycot_sp = list(sm.species)             # pyCOT row order
    pycot_rx = list(sm.reactions)           # pyCOT col order
    S_raw    = np.array(sm)                 # underlying numpy array

    # Sanity-check that every canonical name is present
    missing_sp = [s for s in canonical_species   if s not in pycot_sp]
    missing_rx = [r for r in canonical_reactions if r not in pycot_rx]
    if missing_sp:
        raise ValueError(f"Species not in network: {missing_sp}")
    if missing_rx:
        raise ValueError(f"Reactions not in network: {missing_rx}")

    # Build permutation indices
    row_idx = [pycot_sp.index(s) for s in canonical_species]
    col_idx = [pycot_rx.index(r) for r in canonical_reactions]

    S = S_raw[np.ix_(row_idx, col_idx)]
    return S, canonical_species, canonical_reactions


def print_stoichiometry(S, sp_names, rx_names):
    """Pretty-print the stoichiometric matrix."""
    col_w = 6
    header = f"{'':>4s}  " + "  ".join(f"{r:>{col_w}s}" for r in rx_names)
    sep    = "-" * len(header)
    print(header)
    print(sep)
    for i, sp in enumerate(sp_names):
        row = f"{sp:>4s}  " + "  ".join(f"{int(S[i,j]):>{col_w}d}" for j in range(S.shape[1]))
        print(row)
    print(f"\nMatrix dimensions: {S.shape[0]} species x {S.shape[1]} reactions")


# ===========================================================================
# PART 3 — Custom kinetic functions (pyCOT interface)
# ===========================================================================
#
# Every function follows the pyCOT calling convention:
#
#   f(_, concentrations, species_idx, spec_vector) -> float
#
#   reactants    : list[(species_name, stoich_coeff)] from build_reaction_dict
#   concentrations: np.ndarray — all species concentrations (pyCOT internal order)
#   species_idx  : dict  — {species_name: index_into_concentrations}
#   spec_vector  : list  — ordered parameters for this reaction
#
# Because all species lookups go through the species_idx dict by name, the
# functions are agnostic to the internal species ordering — they work the same
# inside simulation() and in our manual compute_process_vector().
#
# Policy variables (alpha_inst, gamma_sec) are the LAST entries
# in the relevant spec_vectors, so they can be changed before calling
# simulation() without touching the function definitions.

def kin_u1(_, concentrations, species_idx, spec_vector):
    """
    u1: A + Rp => A + Rc
    Corruption: armed groups convert peace resources into conflict resources.
    A is catalytic (not consumed); Rp consumed, Rc produced.
    Inhibited by Rc: high existing Rc stock reduces marginal extraction gain.
    No peace mirror — conflict drains Rp directly; peace cycle cannot do this.
    v_u1 = k_u1 * A * Rp / (Rp + Km_u1)
    spec_vector = [k_u1, Km_u1]
    """
    k, Km = spec_vector
    A  = concentrations[species_idx['A']]
    Rp = concentrations[species_idx['Rp']]
    Rc = concentrations[species_idx['Rc']]
    return k * A * Rp / (Rp + Km)


def kin_u2(_, concentrations, species_idx, spec_vector):
    """
    u2: A => A + G
    Propaganda: armed groups generate grievances; MM saturation in A
    (diminishing returns at high A). Mirrors c2: I=>I+T — same form.
    Symmetric constants: k_u2 = k_c2, Km_u2 = Km_I.
    v_u2 = k_u2 * A / (A + Km_u2)
    spec_vector = [k_u2, Km_u2]
    """
    k, Km = spec_vector
    A = concentrations[species_idx['A']]
    return k * A / (A + Km)


def kin_u3(_, concentrations, species_idx, spec_vector):
    """
    u3: G + Rc => A
    Recruitment: grievances + conflict resources produce armed group members.
    Both G and Rc consumed; MM saturation in G.
    No direct peace mirror (peace cycle grows I via c3 autocatalysis, not a signal-conversion step).
    v_u3 = k_u3 * G * Rc / (G + Km_u3)
    spec_vector = [k_u3, Km_u3]
    """
    k, Km = spec_vector
    G  = concentrations[species_idx['G']]
    Rc = concentrations[species_idx['Rc']]
    return k * G * Rc / (G + Km)


def kin_u4(_, concentrations, species_idx, spec_vector):
    """
    u4: A + I => A
    Institutional assault: armed groups destroy institutional capacity.
    A is catalytic (not consumed); I is destroyed, nothing is gained.
    Purely destructive — symmetric counterpart to c4 (I+A=>I):
      u4 net: I(-1)   conflict destroys peace infrastructure, no resource gain
      c4 net: A(-1)   peace suppresses conflict, no resource gain either
    v_u4 = k_u4 * A * I / (I + Km_u4)
    spec_vector = [k_u4, Km_u4]
    """
    k, Km = spec_vector
    A = concentrations[species_idx['A']]
    I = concentrations[species_idx['I']]
    return k * A * I / (I + Km)


def kin_c1(_, concentrations, species_idx, spec_vector):
    """
    c1: T + Rp => 2Rp
    Trust-driven Rp regeneration: T amplifies Rp; MM saturation in Rp.
    Both T and Rp consumed; net Rp(+1). No conflict mirror — the LV asymmetry.
    Sub-path c2+c1 grows Rp with I as catalyst (I→T via c2, T+Rp→2Rp via c1).
    Rp caps when Rp >> Km_c1 (rate → k_c1*T, bounded via f2 on T).
    LV coupling: u4 steals I → c2 weakens → T falls → c1 slows → Rp collapses.
    v_c1 = k_c1 * T * Rp / (Rp + Km_c1)
    spec_vector = [k_c1, Km_c1]
    """
    k, Km = spec_vector
    T  = concentrations[species_idx['T']]
    Rp = concentrations[species_idx['Rp']]
    return k * T * Rp / (Rp + Km)


def kin_c2(_, concentrations, species_idx, spec_vector):
    """
    c2: I => I + T
    Trust building: institutions generate trust; saturates with I.
    Km_I calibrated near expected I operating range.
    v_c2 = k_c2 * I / (I + Km_I)
    spec_vector = [k_c2, Km_I]
    """
    k, Km_I = spec_vector
    I = concentrations[species_idx['I']]
    return k * I / (I + Km_I)


def kin_c3(_, concentrations, species_idx, spec_vector):
    """
    c3: I + Rp => 2I
    Institutional autocatalysis: existing institutions + peace resources → more
    institutions. I autocatalytic (both substrate and product); Rp consumed.
    MM saturation in I (autocatalytic species in denominator).
    Peace cycle growth: self-sustaining; no net resource loss per I gained.
    v_c3 = alpha_inst * k_c3 * I * Rp / (I + Km_c3)
    spec_vector = [k_c3, Km_c3, alpha_inst]
    alpha_inst is the last entry — change at runtime.
    """
    k, Km, alpha_inst = spec_vector
    I  = concentrations[species_idx['I']]
    Rp = concentrations[species_idx['Rp']]
    return alpha_inst * k * I * Rp / (I + Km)


def kin_c4(_, concentrations, species_idx, spec_vector):
    """
    c4: I + A => I
    Conflict reduction: institutions suppress armed groups. I is catalytic;
    A is eliminated. MM saturation in I, policy lever gamma_sec.
    v_c4 = gamma_sec * k_c4 * A * I / (I + Km_c4_I)
    spec_vector = [k_c4, Km_c4_I, gamma_sec]
    gamma_sec is the last entry — change at runtime.
    """
    k, Km, gamma_sec = spec_vector
    A = concentrations[species_idx['A']]
    I = concentrations[species_idx['I']]
    return gamma_sec * k * A * I / (A + Km)


def kin_f1(_, concentrations, species_idx, spec_vector):
    """
    f1: 2G => ;
    Grievance forgetting: quadratic decay — two units of G required to sustain
    each other in social memory; lone beliefs fade faster.
    v_f1 = k_f1 * G^2
    spec_vector = [k_f1]
    """
    k = spec_vector[0]
    G = concentrations[species_idx['G']]
    return k * G**2


def kin_f2(_, concentrations, species_idx, spec_vector):
    """
    f2: 2T => ;
    Trust forgetting: trust fades quadratically without social reinforcement.
    v_f2 = k_f2 * T^2
    spec_vector = [k_f2]
    """
    k = spec_vector[0]
    T = concentrations[species_idx['T']]
    return k * T**2


# Register custom laws under string keys for use in simulation(additional_laws=...)
# Keys are intentionally prefixed 'kin_' to avoid clashing with pyCOT's built-in
# law names ('mak', 'mmk', 'hill', etc.)
CUSTOM_LAWS = {
    'kin_u1': kin_u1,
    'kin_u2': kin_u2,
    'kin_u3': kin_u3,
    'kin_u4': kin_u4,
    'kin_c1': kin_c1,
    'kin_c2': kin_c2,
    'kin_c3': kin_c3,
    'kin_c4': kin_c4,
    'kin_f1': kin_f1,
    'kin_f2': kin_f2,
}

# Canonical reaction order and which law each reaction uses.
# Must match the order in stoichiometry_matrix().reactions (verified = canonical).
REACTION_NAMES = ['u1', 'u2', 'u3', 'u4', 'c1', 'c2', 'c3', 'c4', 'f1', 'f2']
RATE_NAMES     = ['kin_u1', 'kin_u2', 'kin_u3', 'kin_u4',
                  'kin_c1', 'kin_c2', 'kin_c3', 'kin_c4',
                  'kin_f1', 'kin_f2']

# Default parameter values.
# spec_vector[i] holds the parameters for REACTION_NAMES[i], in the order
# expected by the matching kinetic function above.
#
# Policy variables are the LAST element of each relevant sub-list:
#   c3 → [..., alpha_inst]    index 6
#   c4 → [..., gamma_sec]     index 7
#
# Kinetic constants:
#   u2 ↔ c2 : k=0.4, Km=2.0   (signal generation: catalyst MM-saturated)
#   u4 ↔ c4 : k=0.2, Km=2.0   (I-encounter: A×I/(I+Km), symmetric destruction)
# u1 (no peace mirror): k=0.2, Km_Rc=2.0  — Rp→Rc corruption, inhibited by Rc
# u3 (no peace mirror): k=0.25, Km=3.0
# c1 (no conflict mirror, LV asymmetry): k=0.2, Km=2.0
# c3 (no conflict mirror): k=0.2, Km=2.0, alpha_inst policy lever

DEFAULT_SPEC_VECTOR = [
    [0.05,  2.0],       # u1: A + Rp => A + Rc;      [k_u1,  Km_Rc]              ─  no mirror (corruption)
    [0.2,  2.0],       # u2: A => A + G;             [k_u2,  Km_u2]              ─┐ symmetric with c2
    [0.05, 3.0],       # u3: G + Rc => A;            [k_u3,  Km_u3]              ─  no mirror
    [0.1,  2.0],       # u4: A + I => A;             [k_u4,  Km_u4]              ─┐ symmetric with c4
    [0.1,  2.0],       # c1: T + Rp => 2Rp;          [k_c1,  Km_c1]              ─  LV asymmetry (no conflict mirror)
    [0.25,  2.0],       # c2: I => I + T;              [k_c2,  Km_I]               ─┘ ↑u2
    [0.1,  5.0, 0.4],  # c3: I + Rp => 2I;            [k_c3,  Km_c3, alpha_inst]  ─  no conflict mirror
    [0.1,  2.0, 0.5],  # c4: I + A => I;              [k_c4,  Km_c4, gamma_sec]   ─┘ ↑u4
    [0.01],           # f1: 2G  => ;                 [k_f1]  — grievance forgetting
    [0.01],           # f2: 2T  => ;                 [k_f2]  — trust forgetting
]


def make_spec_vector(alpha_inst=0.4, gamma_sec=0.5):
    """
    Build a spec_vector from scratch with the given policy allocations,
    using all other kinetic constants from DEFAULT_SPEC_VECTOR.
    Returns a new list (does not mutate DEFAULT_SPEC_VECTOR).

    Policy levers:
      alpha_inst — scales c3 (I+Rp=>2I institutional autocatalysis)  index 6
      gamma_sec  — scales c4 (I+A=>I  conflict suppression)          index 7
    """
    sv = [list(p) for p in DEFAULT_SPEC_VECTOR]  # deep copy
    sv[6][-1] = alpha_inst   # c3 (index 6): last element is alpha_inst
    sv[7][-1] = gamma_sec    # c4 (index 7): last element is gamma_sec
    return sv


def check_policy_constraint(spec_vector):
    """
    Warn if alpha_inst or gamma_sec are outside [0, 1].
    Returns a clipped copy if violated, otherwise returns spec_vector unchanged.

    alpha_inst (sv[6][-1]): scales c3 institutional autocatalysis.
    gamma_sec  (sv[7][-1]): scales c4 conflict suppression.
    Both are independent multipliers — no sum constraint.
    """
    alpha_inst = spec_vector[6][-1]
    gamma_sec  = spec_vector[7][-1]
    violated = []
    if not (0.0 <= alpha_inst <= 1.0 + 1e-9):
        violated.append(f"alpha_inst={alpha_inst:.3f}")
    if not (0.0 <= gamma_sec  <= 1.0 + 1e-9):
        violated.append(f"gamma_sec={gamma_sec:.3f}")
    if violated:
        warnings.warn(
            f"Policy parameter(s) outside [0,1]: {', '.join(violated)}. "
            "Clipping to [0,1].",
            UserWarning, stacklevel=2,
        )
        sv = [list(p) for p in spec_vector]
        sv[6][-1] = float(np.clip(alpha_inst, 0.0, 1.0))
        sv[7][-1] = float(np.clip(gamma_sec,  0.0, 1.0))
        return sv
    return spec_vector


# ===========================================================================
# PART 4 — Computing the process vector and state change
# ===========================================================================

def compute_process_vector(x, spec_vector, rn_dict, species_idx):
    """
    Evaluate all seven reaction rates at state x using the custom kinetics.

    This calls each kinetic function exactly as simulation() would internally,
    using the pyCOT calling convention.

    Parameters
    ----------
    x           : np.ndarray (n_species,) — concentrations in canonical order
    spec_vector : list-of-lists — parameters per reaction (canonical order)
    rn_dict     : dict from build_reaction_dict(rn) — {rx: (reactants, products)}
    species_idx : dict — {species_name: index_into_x}

    Returns
    -------
    v : np.ndarray (10,) — [v_u1, v_u2, v_u3, v_u4,
                             v_c1, v_c2, v_c3, v_c4, v_f1, v_f2]
    """
    v = np.zeros(len(REACTION_NAMES))
    for i, (rx, law_name) in enumerate(zip(REACTION_NAMES, RATE_NAMES)):
        reactants  = rn_dict[rx][0]           # list[(species_name, coeff)]
        fn         = CUSTOM_LAWS[law_name]
        v[i]       = fn(reactants, x, species_idx, spec_vector[i])
    return v


def compute_state_change(x, spec_vector, S, sp_names, rn_dict, species_idx, dt=1.0):
    """
    Compute and display the net change in each species over one timestep.

    delta_x = S @ v * dt

    Returns v and delta_x.
    """
    v       = compute_process_vector(x, spec_vector, rn_dict, species_idx)
    delta_x = S @ v * dt

    print("\n--- Process vector v ---")
    for name, rate in zip(REACTION_NAMES, v):
        print(f"  {name:>4s}: {rate:+.6f}")

    print(f"\n--- State change delta_x (dt={dt}) ---")
    for name, dx in zip(sp_names, delta_x):
        print(f"  {name:>4s}: {dx:+.6f}")

    return v, delta_x


# ===========================================================================
# PART 5 — Simulation using pyCOT's ODE engine and comparative plot
# ===========================================================================

def run_simulation(rn, spec_vector, x0_pycot, t_end=50.0, n_steps=500):
    """
    Run an ODE simulation via pyCOT's simulation() function.

    Parameters
    ----------
    rn           : ReactionNetwork
    spec_vector  : list-of-lists — parameters per reaction (canonical/pyCOT order)
    x0_pycot     : list — initial concentrations in pyCOT's internal species order
    t_end        : float
    n_steps      : int

    Returns
    -------
    ts_df  : pd.DataFrame  columns = ['Time', species...]
    flx_df : pd.DataFrame  columns = ['Time', reactions...]
    """
    spec_vector = check_policy_constraint(spec_vector)
    ts_df, flx_df = simulation(
        rn,
        rate           = RATE_NAMES,        # one law name per reaction (canonical order)
        spec_vector    = spec_vector,        # params per reaction
        x0             = x0_pycot,          # initial state (pyCOT species order)
        t_span         = (0.0, t_end),
        n_steps        = n_steps,
        additional_laws= CUSTOM_LAWS,       # inject our custom kinetics
        method         = 'LSODA',           # adaptive stiff-capable solver
        verbose        = False,
    )
    return ts_df, flx_df


def plot_comparison(ts_A, ts_B, sv_A, sv_B, save_path=None):
    """
    Four-panel (2×2) time-series comparison of two simulations.

    Each panel compares one conflict species against its natural peace counterpart,
    keeping the scales directly comparable within the pair:
      Panel 1 (top-left)  — G  vs T   (grievances vs trust)
      Panel 2 (top-right) — V  vs I   (violence vs institutions)
      Panel 3 (bot-left)  — Rc vs Rp  (conflict resources vs peace resources)
      Panel 4 (bot-right) — Σ(G+A+Rc) vs Σ(I+T+Rp)  (aggregate balance)

    Within every panel: dashed = Sim A (no policy), solid = Sim B (peace policy).
    Conflict species are drawn in warm colours (reds/orange),
    peace species in cool colours (greens/blue).
    """
    # Pair definitions: (conflict_sym, conflict_label, conflict_color,
    #                    peace_sym,    peace_label,    peace_color,
    #                    panel_title)
    pairs = [
        ('G',  'Grievances',        '#C0392B',
         'T',  'Trust',             '#2980B9',
         'G vs T  (grievances vs trust)'),
        ('A',  'Armed groups',      '#E74C3C',
         'I',  'Institutions',      '#27AE60',
         'A vs I  (armed groups vs institutions)'),
        ('Rc', 'Conflict resources','#F39C12',
         'Rp', 'Peace resources',   '#1ABC9C',
         'Rc vs Rp  (conflict vs peace resources)'),
    ]

    time = ts_A['Time'].values

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_flat = axes.flatten()   # [top-left, top-right, bot-left, bot-right]

    for ax, (csym, clabel, ccol, psym, plabel, pcol, title) in zip(ax_flat[:3], pairs):
        ax.plot(time, ts_A[csym], color=ccol, linestyle='--', linewidth=1.5,
                label=f'{clabel} (A)')
        ax.plot(time, ts_B[csym], color=ccol, linestyle='-',  linewidth=2.0,
                label=f'{clabel} (B)')
        ax.plot(time, ts_A[psym], color=pcol, linestyle='--', linewidth=1.5,
                label=f'{plabel} (A)')
        ax.plot(time, ts_B[psym], color=pcol, linestyle='-',  linewidth=2.0,
                label=f'{plabel} (B)')
        ax.set_ylabel('Concentration', fontsize=9)
        ax.set_xlabel('Time', fontsize=9)
        ax.set_title(title, fontsize=9)
        ax.legend(fontsize=8, ncol=2)
        ax.grid(True, alpha=0.3)

    # Panel 4 — aggregate conflict vs peace totals
    ax4 = ax_flat[3]
    conf_A = ts_A['G'] + ts_A['A'] + ts_A['Rc']
    conf_B = ts_B['G'] + ts_B['A'] + ts_B['Rc']
    peac_A = ts_A['I'] + ts_A['T'] + ts_A['Rp']
    peac_B = ts_B['I'] + ts_B['T'] + ts_B['Rp']

    ax4.plot(time, conf_A, color='#C0392B', linestyle='--', linewidth=1.5,
             label='Σ conflict (A)')
    ax4.plot(time, conf_B, color='#C0392B', linestyle='-',  linewidth=2.0,
             label='Σ conflict (B)')
    ax4.plot(time, peac_A, color='#27AE60', linestyle='--', linewidth=1.5,
             label='Σ peace (A)')
    ax4.plot(time, peac_B, color='#27AE60', linestyle='-',  linewidth=2.0,
             label='Σ peace (B)')
    ax4.set_ylabel('Total concentration', fontsize=9)
    ax4.set_xlabel('Time', fontsize=9)
    ax4.set_title('Aggregate balance  Σ(G+A+Rc) vs Σ(I+T+Rp)', fontsize=9)
    ax4.legend(fontsize=8, ncol=2)
    ax4.grid(True, alpha=0.3)

    # Suptitle with policy parameters
    a_inst_A, g_sec_A = sv_A[6][-1], sv_A[7][-1]
    a_inst_B, g_sec_B = sv_B[6][-1], sv_B[7][-1]
    fig.suptitle(
        f"Conflict-Peace Dynamics — pairwise comparison\n"
        f"A (--): α_inst={a_inst_A}, γ_sec={g_sec_A}  |  "
        f"B (—):  α_inst={a_inst_B}, γ_sec={g_sec_B}",
        fontsize=9, y=1.01,
    )
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved: {save_path}")

    plt.show()


# ===========================================================================
# Phase portrait — deterministic policy comparison
# ===========================================================================

def plot_phase_det(ts_A, ts_B, sv_A, sv_B, save_path=None):
    """
    Phase portrait in (A, I) space for the two deterministic simulations.
    Blues  (light→dark, dashed) = Sim A (no policy)
    Oranges (light→dark, solid) = Sim B (peace policy)
    Same colour scheme as plot_phase_comparison() in script_basic_game.py.
    """
    a_inst_A, g_sec_A = sv_A[6][-1], sv_A[7][-1]
    a_inst_B, g_sec_B = sv_B[6][-1], sv_B[7][-1]

    A_a = ts_A['A'].values;  I_a = ts_A['I'].values
    A_b = ts_B['A'].values;  I_b = ts_B['I'].values
    n_a = len(A_a);          n_b = len(A_b)

    cmap_A = plt.colormaps['Blues']
    cmap_B = plt.colormaps['Oranges']
    def col_A(i): return cmap_A(0.35 + 0.65 * i / max(n_a - 1, 1))
    def col_B(i): return cmap_B(0.35 + 0.65 * i / max(n_b - 1, 1))

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.suptitle(
        f'Phase Portrait: Armed Groups vs Institutions\n'
        f'Blues (--): α_inst={a_inst_A}, γ_sec={g_sec_A}  (no policy)   |   '
        f'Oranges (—): α_inst={a_inst_B}, γ_sec={g_sec_B}  (peace policy)',
        fontsize=9, fontweight='bold'
    )

    # Trajectory A — Blues, dashed
    for i in range(n_a - 1):
        ax.plot(A_a[i:i+2], I_a[i:i+2], color=col_A(i), lw=1.3, ls='--', alpha=0.75)

    # Trajectory B — Oranges, solid
    for i in range(n_b - 1):
        ax.plot(A_b[i:i+2], I_b[i:i+2], color=col_B(i), lw=1.8, ls='-',  alpha=0.85)

    # Arrows — Sim A
    gap_a = max(1, n_a // 15)
    for i in range(0, n_a - 1, gap_a):
        dA, dI = A_a[i+1] - A_a[i], I_a[i+1] - I_a[i]
        if abs(dA) + abs(dI) > 1e-9:
            ax.annotate('', xy=(A_a[i]+dA, I_a[i]+dI), xytext=(A_a[i], I_a[i]),
                        arrowprops=dict(arrowstyle='->', color='steelblue', lw=0.8))

    # Arrows — Sim B
    gap_b = max(1, n_b // 15)
    for i in range(0, n_b - 1, gap_b):
        dA, dI = A_b[i+1] - A_b[i], I_b[i+1] - I_b[i]
        if abs(dA) + abs(dI) > 1e-9:
            ax.annotate('', xy=(A_b[i]+dA, I_b[i]+dI), xytext=(A_b[i], I_b[i]),
                        arrowprops=dict(arrowstyle='->', color='darkorange', lw=0.8))

    # Start (shared) and end markers
    ax.scatter(A_a[0],  I_a[0],  color='black',      s=80,  zorder=6,
               marker='o', label='Start (shared)')
    ax.scatter(A_a[-1], I_a[-1], color='steelblue',  s=90,  zorder=6,
               marker='X', label='End — no policy')
    ax.scatter(A_b[-1], I_b[-1], color='darkorange', s=110, zorder=6,
               marker='*', label='End — peace policy')

    # Two colourbars — same layout as stochastic script
    sm_A = plt.cm.ScalarMappable(cmap=cmap_A, norm=plt.Normalize(vmin=0, vmax=1))
    sm_A.set_array([])
    sm_B = plt.cm.ScalarMappable(cmap=cmap_B, norm=plt.Normalize(vmin=0, vmax=1))
    sm_B.set_array([])
    cb_A = fig.colorbar(sm_A, ax=ax, pad=0.02,  fraction=0.035)
    cb_A.set_label('No policy  (early→late)', fontsize=8)
    cb_B = fig.colorbar(sm_B, ax=ax, pad=0.12, fraction=0.035)
    cb_B.set_label('Peace policy  (early→late)', fontsize=8)

    ax.set_xlabel('A  (Armed groups)', fontsize=10)
    ax.set_ylabel('I  (Institutions)', fontsize=10)
    ax.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1.0, 1.0))
    ax.grid(True, alpha=0.25)

    plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved: {save_path}")
    plt.show()


# ===========================================================================
# PART 6 — Process vector testing utility
# ===========================================================================

def test_process_vector(v_test, x_current, S, sp_names, label='test'):
    """
    Probe an arbitrary process vector to build intuition about the
    stoichiometric structure.

    Parameters
    ----------
    v_test    : np.ndarray (10,) — reaction rates to test (specified directly
                OR obtained from compute_process_vector for a given state)
    x_current : np.ndarray (6,) — reference species concentrations
    S         : np.ndarray (6,10) — stoichiometric matrix
    sp_names  : list[str]
    label     : str — descriptive name for this test
    """
    delta_x = S @ v_test   # net change per species (dt = 1 implied)

    def _mode(dx):
        if   dx >  1e-12: return 'producing'
        elif dx < -1e-12: return 'consuming'
        else:             return 'neutral'

    print(f"\n{'='*60}")
    print(f"  Process vector test: '{label}'")
    print(f"{'='*60}")

    print("\n  Reaction rates (v):")
    for r, vr in zip(REACTION_NAMES, v_test):
        print(f"    {r:>4s}: {vr:+.4f}")

    print(f"\n  {'Species':<8} {'x_current':>10} {'delta_x':>10}  {'mode':<12}")
    print(f"  {'-'*44}")
    for sp, xc, dx in zip(sp_names, x_current, delta_x):
        print(f"  {sp:<8} {xc:>10.4f} {dx:>10.4f}  {_mode(dx):<12}")

    pos  = np.sum(delta_x >  1e-12)
    neg  = np.sum(delta_x < -1e-12)
    zero = np.sum(np.abs(delta_x) <= 1e-12)

    if zero == len(sp_names):
        overall = 'STEADY-STATE (all delta_x == 0)'
    elif neg == 0 and pos > 0:
        overall = 'OVERPRODUCTION (all delta_x >= 0, at least one > 0)'
    elif pos == 0 and neg > 0:
        overall = 'PROBLEM (all delta_x <= 0, at least one < 0)'
    else:
        overall = 'CHALLENGE (mixed signs)'

    print(f"\n  Overall mode: {overall}")

    overproduced = [sp_names[i] for i in range(len(sp_names)) if delta_x[i] >  1e-12]
    depleted     = [sp_names[i] for i in range(len(sp_names)) if delta_x[i] < -1e-12]
    if overproduced:
        print(f"  Overproduced : {', '.join(overproduced)}")
    if depleted:
        print(f"  Depleted     : {', '.join(depleted)}")
    print()


# ===========================================================================
# PART 7 — Productive strength analysis
# ===========================================================================

def compute_strengths_from_v(v):
    """
    Extract six productive capacity scores from a process vector v (length 10).

    Canonical reaction indices:
      u1=0  u2=1  u3=2  u4=3  c1=4  c2=5  c3=6  c4=7  f1=8  f2=9

    Productive capacities
    ---------------------
    Resource production (gain own resource without losing any group species):
      C_res = v[u1]             — corruption rate A+Rp→A+Rc; A catalytic;
                                  single reaction, no cycle needed.
                                  Inhibited by Rc so self-limiting when Rc is high.
      P_res = min(v[c1],v[c2]) — bottleneck of c2+c1 sub-path:
                                  c2: I→I+T (I catalytic), c1: T+Rp→2Rp
                                  Net: Rp(+1), I fully catalytic; no group cost.

    Group strength (autocatalytic cycle throughput — min over bottleneck step):
      C_grp = min(v[u1],v[u2],v[u3]) — conflict cycle, net A(+1) per turn
                                        costs Rp(-1) each turn (parasitic).
      P_grp = min(v[c1],v[c2],v[c3]) — peace cycle, net I(+1) per turn;
                                        T and Rp internally balanced (self-sustaining).

    Enemy destruction (degrade the opposing force):
      C_dst = v[u4]  — A+I→A: institutional assault (I destroyed, nothing gained)
      P_dst = v[c4]  — I+A→I: armed-group suppression (A destroyed, nothing gained)

    Advantage scores (positive = conflict leads that dimension):
      A_res = C_res − P_res
      A_grp = C_grp − P_grp
      A_dst = C_dst − P_dst
    """
    v = np.asarray(v, dtype=float)
    C_res = float(v[0])
    P_res = float(min(v[4], v[5]))
    C_grp = float(min(v[0], v[1], v[2]))
    P_grp = float(min(v[4], v[5], v[6]))
    C_dst = float(v[3])
    P_dst = float(v[7])
    return dict(
        C_res=C_res, P_res=P_res,
        C_grp=C_grp, P_grp=P_grp,
        C_dst=C_dst, P_dst=P_dst,
        A_res=C_res - P_res,
        A_grp=C_grp - P_grp,
        A_dst=C_dst - P_dst,
    )


def analyse_productive_strengths(flx_df, save_path=None):
    """
    Productive strength analysis — Part 7.

    Reads the flux time series from a simulation (flx_df from simulation())
    and produces three complementary figures.

    FIGURE 1 — Strength time series  (3 × 1 stacked panels)
        Each panel tracks one productive dimension over time.
        Conflict strength (warm red) vs peace strength (cool blue) as solid
        lines; the shaded band between them is red when conflict leads and
        blue when peace leads.
          Resources     : C_res = v(u1)           vs  P_res = min(v(c1), v(c2))
          Group strength: C_grp = min(v(u1,u2,u3)) vs  P_grp = min(v(c1,c2,c3))
          Destruction   : C_dst = v(u4)           vs  P_dst = v(c4)

    FIGURE 2 — Advantage portrait  (hexbin + temporal trajectory)
        x-axis = A_res (resource advantage)
        y-axis = A_grp (group-strength advantage)
        Hexagonal bins coloured by median A_dst within each bin (RdBu_r).
        Subsampled trajectory overlaid as a grey→black path (early→late).
        The four quadrant zones are annotated with their strategic meaning.
        ◆ = start of simulation;  ★ = end.

    FIGURE 3 — One-step momentum field + strategic configuration distribution
        LEFT: Arrows in (A_res, A_grp) advantage space. Each arrow originates
              at the current position and points in the direction of the
              one-step change (δA_res, δA_grp). Arrow colour = δA_dst
              (how the destruction advantage changes in that step).
              Arrows are unit-normalised (uniform length) and subsampled
              to ~400 for readability.
        RIGHT: Horizontal bar chart showing the fraction of simulation time
               spent in each of the 8 binary strategic configurations
               defined by sign(A_res, A_grp, A_dst).
               Format: Res|Grp|Dst — C = conflict leads, P = peace leads.
    """
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize

    time = flx_df['Time'].values
    n    = len(time)

    # ------------------------------------------------------------------
    # Compute strength arrays from flux data
    # ------------------------------------------------------------------
    V_traj = np.column_stack([flx_df[rx].values for rx in REACTION_NAMES])  # (n,10)
    S_list = [compute_strengths_from_v(V_traj[i]) for i in range(n)]

    C_res = np.array([s['C_res'] for s in S_list])
    P_res = np.array([s['P_res'] for s in S_list])
    C_grp = np.array([s['C_grp'] for s in S_list])
    P_grp = np.array([s['P_grp'] for s in S_list])
    C_dst = np.array([s['C_dst'] for s in S_list])
    P_dst = np.array([s['P_dst'] for s in S_list])
    A_res = C_res - P_res
    A_grp = C_grp - P_grp
    A_dst = C_dst - P_dst

    # ==================================================================
    # FIGURE 1 — Strength time series
    # ==================================================================
    fig1, axes1 = plt.subplots(3, 1, figsize=(13, 11), sharex=True)
    dims = [
        ('Resources',
         C_res, P_res,
         'C_res  v(u1): A+Rp→A+Rc  (corruption)',
         'P_res  min(v(c1),v(c2)): c2→c1 Rp-bootstrap bottleneck'),
        ('Group strength',
         C_grp, P_grp,
         'C_grp  min(v(u1,u2,u3)): conflict-cycle throughput  [costs Rp each turn]',
         'P_grp  min(v(c1,c2,c3)): peace-cycle throughput  [self-sustaining]'),
        ('Enemy destruction',
         C_dst, P_dst,
         'C_dst  v(u4): A+I→A  (institutional assault)',
         'P_dst  v(c4): I+A→I  (armed-group suppression)'),
    ]
    for ax, (title, Cc, Pc, clabel, plabel) in zip(axes1, dims):
        ax.fill_between(time, Cc, Pc, where=(Cc >= Pc),
                        alpha=0.20, color='#C0392B', label='conflict leads')
        ax.fill_between(time, Cc, Pc, where=(Cc <  Pc),
                        alpha=0.20, color='#2980B9', label='peace leads')
        ax.plot(time, Cc, color='#C0392B', linewidth=1.8, label=clabel)
        ax.plot(time, Pc, color='#2980B9', linewidth=1.8, label=plabel)
        ax.axhline(0, color='black', linewidth=0.6)
        ax.set_ylabel('Strength  (flux units)', fontsize=9)
        ax.set_title(title, fontsize=9, fontweight='bold')
        ax.legend(fontsize=7, ncol=2, loc='upper right')
        ax.grid(True, alpha=0.2)
    axes1[-1].set_xlabel('Time', fontsize=9)
    fig1.suptitle(
        'Productive strength over time\n'
        'Red = conflict capacity  |  Blue = peace capacity  |  '
        'Shaded area = who leads each dimension',
        fontsize=10,
    )
    fig1.tight_layout()

    # ==================================================================
    # FIGURE 2 — Advantage portrait
    # ==================================================================
    fig2, ax2 = plt.subplots(figsize=(10, 8))

    # Hexbin: density coloured by median A_dst in each bin
    hb = ax2.hexbin(A_res, A_grp, C=A_dst,
                    reduce_C_function=np.median,
                    gridsize=30, cmap='RdBu_r', mincnt=1,
                    alpha=0.60, linewidths=0.2)
    dst_max = max(abs(A_dst.min()), abs(A_dst.max())) + 1e-12
    hb.set_clim(-dst_max, dst_max)
    cbar2 = fig2.colorbar(hb, ax=ax2, pad=0.01, shrink=0.70)
    cbar2.set_label('Median A_dst in bin\n(red = conflict destroys faster, '
                    'blue = peace suppresses faster)', fontsize=8)

    # Temporal trajectory (subsampled)
    step_t = max(1, n // 800)
    idx_t  = np.arange(0, n, step_t)
    t_norm = (time[idx_t] - time[0]) / (time[-1] - time[0] + 1e-12)
    pts    = np.c_[A_res[idx_t], A_grp[idx_t]]
    segs   = np.stack([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap='Greys', norm=Normalize(0.0, 1.4),
                        linewidth=1.2, alpha=0.70, zorder=4)
    lc.set_array(t_norm[:-1])
    ax2.add_collection(lc)

    # Start / end markers
    ax2.scatter(A_res[0],  A_grp[0],  s=160, c='lime',  marker='D', zorder=8,
                edgecolors='black', linewidths=0.8, label=f't = {time[0]:.0f}  ◆ start')
    ax2.scatter(A_res[-1], A_grp[-1], s=160, c='black', marker='*', zorder=8,
                label=f't = {time[-1]:.0f}  ★ end')

    # Quadrant dividers
    ax2.axhline(0, color='white', linewidth=1.2, linestyle='--', alpha=0.75)
    ax2.axvline(0, color='white', linewidth=1.2, linestyle='--', alpha=0.75)

    # Quadrant labels (placed in data-relative fractions)
    ax2.autoscale_view()
    xl, xr = ax2.get_xlim(); yb, yt = ax2.get_ylim()
    _qkw = dict(ha='center', va='center', fontsize=8, style='italic', alpha=0.80,
                bbox=dict(boxstyle='round,pad=0.18', facecolor='white', alpha=0.45))
    ax2.text(xl + 0.80*(xr-xl), yb + 0.88*(yt-yb),
             'C wins\nRes + Grp', color='#7B241C', **_qkw)
    ax2.text(xl + 0.20*(xr-xl), yb + 0.88*(yt-yb),
             'P wins Res\nC wins Grp', color='#6C3483', **_qkw)
    ax2.text(xl + 0.80*(xr-xl), yb + 0.12*(yt-yb),
             'C wins Res\nP wins Grp', color='#1A5276', **_qkw)
    ax2.text(xl + 0.20*(xr-xl), yb + 0.12*(yt-yb),
             'P wins\nRes + Grp', color='#1A5276', **_qkw)

    ax2.set_xlabel('A_res  =  C_res − P_res   (resource advantage)', fontsize=10)
    ax2.set_ylabel('A_grp  =  C_grp − P_grp   (group-strength advantage)', fontsize=10)
    ax2.set_title(
        'Advantage portrait — trajectory in 3-D productive-strength space\n'
        'Hexbin density (median destruction advantage)  |  '
        'Grey path = temporal trajectory (light→dark = early→late)  |  ◆ start  ★ end',
        fontsize=9,
    )
    ax2.legend(fontsize=8, loc='lower right')
    ax2.grid(False)
    fig2.tight_layout()

    # ==================================================================
    # FIGURE 3 — Momentum field + configuration distribution
    # ==================================================================
    fig3, (ax3L, ax3R) = plt.subplots(
        1, 2, figsize=(15, 8), gridspec_kw={'width_ratios': [2.5, 1]}
    )

    # One-step deltas in advantage space
    dA_res = np.diff(A_res)
    dA_grp = np.diff(A_grp)
    dA_dst = np.diff(A_dst)
    mag    = np.sqrt(dA_res**2 + dA_grp**2) + 1e-14

    # Subsample to ~400 arrows
    step_q = max(1, (n - 1) // 400)
    iq     = np.arange(0, n - 1, step_q)
    u_arr  = dA_res[iq] / mag[iq]          # unit direction x
    v_arr  = dA_grp[iq] / mag[iq]          # unit direction y
    dst_arr = dA_dst[iq]
    dst_lim = max(abs(dst_arr.min()), abs(dst_arr.max())) + 1e-12

    qv = ax3L.quiver(
        A_res[iq], A_grp[iq],
        u_arr, v_arr,
        dst_arr,
        cmap='RdBu_r', norm=Normalize(-dst_lim, dst_lim),
        scale=18, scale_units='width',
        width=0.004, headwidth=4, headlength=5,
        alpha=0.75, zorder=3,
    )
    cbar3 = fig3.colorbar(qv, ax=ax3L, pad=0.01, shrink=0.80)
    cbar3.set_label('δA_dst  (one-step change in destruction advantage)\n'
                    'red = conflict gains destruction lead, '
                    'blue = peace gains suppression lead', fontsize=8)

    ax3L.axhline(0, color='grey', linewidth=0.9, linestyle='--')
    ax3L.axvline(0, color='grey', linewidth=0.9, linestyle='--')
    ax3L.scatter(A_res[0],  A_grp[0],  s=140, c='lime',  marker='D',
                 zorder=8, edgecolors='black', linewidths=0.8)
    ax3L.scatter(A_res[-1], A_grp[-1], s=140, c='black', marker='*', zorder=8)
    ax3L.set_xlabel('A_res  (resource advantage)', fontsize=10)
    ax3L.set_ylabel('A_grp  (group-strength advantage)', fontsize=10)
    ax3L.set_title(
        'One-step momentum field in advantage space\n'
        'Each arrow = direction of advantage change per time step  '
        '|  Colour = δA_dst  |  ◆ start  ★ end',
        fontsize=9,
    )
    ax3L.grid(True, alpha=0.2)

    # Strategic configuration distribution (8 binary categories)
    def _sign(x): return 'C' if x >= 0 else 'P'
    keys = ['C|C|C','C|C|P','C|P|C','C|P|P','P|C|C','P|C|P','P|P|C','P|P|P']
    cat_labels = [
        'C|C|C  (conflict all three)',
        'C|C|P  (C res+grp,  P dst)',
        'C|P|C  (C res+dst,  P grp)',
        'C|P|P  (C res only)',
        'P|C|C  (P res,  C grp+dst)',
        'P|C|P  (P res+dst, C grp)',
        'P|P|C  (P res+grp, C dst)',
        'P|P|P  (peace all three)',
    ]
    bar_colors = [
        '#7B241C','#C0392B','#E74C3C','#F1948A',
        '#AED6F1','#2980B9','#1A5276','#0B3D91',
    ]
    cats  = [f"{_sign(A_res[i])}|{_sign(A_grp[i])}|{_sign(A_dst[i])}"
             for i in range(n)]
    fracs = np.array([cats.count(k) / n for k in keys])

    bars = ax3R.barh(cat_labels[::-1], fracs[::-1],
                     color=bar_colors[::-1], alpha=0.88, edgecolor='white')
    ax3R.set_xlabel('Fraction of simulation time', fontsize=9)
    ax3R.set_title('Strategic configuration\ndistribution\n(Res | Grp | Dst)',
                   fontsize=9)
    ax3R.set_xlim(0, max(fracs.max() * 1.20, 0.05))
    for bar, frac in zip(bars, fracs[::-1]):
        if frac > 0.005:
            ax3R.text(frac + 0.003, bar.get_y() + bar.get_height() / 2,
                      f'{frac:.1%}', va='center', fontsize=8)
    ax3R.grid(True, axis='x', alpha=0.25)

    fig3.suptitle(
        'Productive advantage dynamics\n'
        'Res|Grp|Dst — C = conflict leads, P = peace leads that dimension',
        fontsize=10,
    )
    fig3.tight_layout()

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        p1 = save_path.replace('.png', '_strengths.png')
        p2 = save_path.replace('.png', '_portrait.png')
        p3 = save_path.replace('.png', '_momentum.png')
        fig1.savefig(p1, dpi=150, bbox_inches='tight')
        fig2.savefig(p2, dpi=150, bbox_inches='tight')
        fig3.savefig(p3, dpi=150, bbox_inches='tight')
        print(f"  Strengths time series : {p1}")
        print(f"  Advantage portrait    : {p2}")
        print(f"  Momentum field        : {p3}")

    plt.show()

    # Console summary
    s_fin   = S_list[-1]
    dom_key = keys[int(np.argmax(fracs))]
    print(f"\n  --- Productive strengths at t = {time[-1]:.1f} ---")
    for dim, c_val, p_val, a_val in [
        ('Resources  ', s_fin['C_res'], s_fin['P_res'], s_fin['A_res']),
        ('Group str. ', s_fin['C_grp'], s_fin['P_grp'], s_fin['A_grp']),
        ('Destruction', s_fin['C_dst'], s_fin['P_dst'], s_fin['A_dst']),
    ]:
        leader = 'C leads' if a_val >= 0 else 'P leads'
        print(f"  {dim}: C={c_val:.5f}  P={p_val:.5f}  "
              f"Adv={a_val:+.5f}  ({leader})")
    print(f"\n  Dominant configuration : {dom_key}  ({fracs.max():.1%} of time)")
    print(f"  Configuration breakdown:")
    for k, f, lbl in zip(keys, fracs, cat_labels):
        if f > 0.0:
            print(f"    {lbl}  {f:.1%}")


# MAIN — run Parts 2, 4, 5, 6 sequentially
# ===========================================================================

if __name__ == '__main__':

    # -----------------------------------------------------------------------
    # PART 2: Load network, inspect structure
    # -----------------------------------------------------------------------
    DATA_FILE = os.path.join(_SCRIPT_DIR, 'data', 'Basic_example.txt')

    print("=" * 60)
    print("PART 2 — Network loading and stoichiometric structure")
    print("=" * 60)

    rn = load_network(DATA_FILE)

    # Canonical orderings as defined in the paper
    SPECIES   = ['G', 'A', 'I', 'T', 'Rc', 'Rp']
    REACTIONS = ['u1', 'u2', 'u3', 'u4', 'c1', 'c2', 'c3', 'c4', 'f1', 'f2']

    S, sp_names, rx_names = get_stoichiometry(rn, SPECIES, REACTIONS)

    # Print pyCOT's internal ordering (may differ from canonical)
    sm = rn.stoichiometry_matrix()
    print(f"\npyCOT internal species order   : {list(sm.species)}")
    print(f"pyCOT internal reaction order  : {list(sm.reactions)}")
    print(f"\nCanonical species order  (rows): {sp_names}")
    print(f"Canonical reaction order (cols): {rx_names}")
    print("\nStoichiometric matrix S:")
    print_stoichiometry(S, sp_names, rx_names)

    assert S.shape == (6, 10), f"Expected (6,10), got {S.shape}"
    print("\nDimension check: 6 species x 10 reactions  OK")

    # Build reaction dictionary (reactants/products per reaction, by name)
    rn_dict = build_reaction_dict(rn)
    print("\nReaction structure (from build_reaction_dict):")
    for rx in REACTIONS:
        reactants, products = rn_dict[rx]
        r_str = " + ".join(f"{c}{s}" if c != 1 else s for s, c in reactants) or "0"
        p_str = " + ".join(f"{c}{s}" if c != 1 else s for s, c in products)  or "0"
        print(f"  {rx}: {r_str} => {p_str}")

    # species_idx maps species name -> index in canonical x vector.
    # pyCOT's stoichiometry_matrix().species == rn.species() order (both sorted
    # by node insertion index), so the canonical order we verified above is
    # exactly what simulation() uses internally for x0 and concentrations.
    species_idx_canonical = {sp: i for i, sp in enumerate(sp_names)}

    # -----------------------------------------------------------------------
    # PART 4: Process vector and state change at the initial conflict state
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("PART 4 — Process vector and state change (default params)")
    print("=" * 60)

    # Initial state  [G, A, I, T, Rc, Rp]
    # Balanced start: equal amounts of all species so both conflict and peace
    # have something to work with from the outset.
    # Canonical order: [G, A, I, T, Rc, Rp]
    # Peace-favorable but balanced start: I > A, Rp > Rc, small social memory.
    # Same x0 used in script_basic_game.py for cross-kinetics comparison.
    x0_conflict = np.array([2.0, 2.0, 4.0, 2.0, 2.0, 6.0])

    v0, dx0 = compute_state_change(
        x0_conflict,
        DEFAULT_SPEC_VECTOR,
        S,
        sp_names,
        rn_dict,
        species_idx_canonical,
        dt=1.0,
    )

    # -----------------------------------------------------------------------
    # PART 5: Simulate two scenarios and compare
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("PART 5 — Simulation: baseline vs peacebuilding policy")
    print("=" * 60)

    # Simulation A — no policy (all allocation parameters = 0)
    sv_A = make_spec_vector(alpha_inst=0.0, gamma_sec=0.0)

    # Simulation B — active peacebuilding (default allocations)
    sv_B = make_spec_vector(alpha_inst=0.4, gamma_sec=0.5)

    # Build x0 in pyCOT's internal species order.
    # x0_conflict is in canonical order (SPECIES). pyCOT's internal order
    # follows first-appearance in the file (e.g. [G, Rc, V, I, Rp, T]),
    # which may differ from canonical. We reindex by name to be safe.
    pycot_sp_order = list(sm.species)   # sm already computed in Part 2
    x0_canon_dict  = {sp: float(v) for sp, v in zip(SPECIES, x0_conflict)}
    x0_list = [x0_canon_dict[sp] for sp in pycot_sp_order]
    print(f"\nx0 mapping to pyCOT internal order:")
    for sp, val in zip(pycot_sp_order, x0_list):
        print(f"  {sp}: {val}")

    T      = 300.0
    NSTEPS = 1500   # 2000 output points over T=200 → dt_output = 0.1

    print(f"\nRunning Simulation A (no policy)    T={T}, n_steps={NSTEPS} ...")
    ts_A, flx_A = run_simulation(rn, sv_A, x0_list, t_end=T, n_steps=NSTEPS)

    print(f"Running Simulation B (peace policy) T={T}, n_steps={NSTEPS} ...")
    ts_B, flx_B = run_simulation(rn, sv_B, x0_list, t_end=T, n_steps=NSTEPS)

    PLOT_PATH = os.path.normpath(
        os.path.join(_SCRIPT_DIR, '..', 'outputs', 'basic_example_simulation.png')
    )
    plot_comparison(ts_A, ts_B, sv_A, sv_B, save_path=PLOT_PATH)

    PHASE_PATH = os.path.normpath(
        os.path.join(_SCRIPT_DIR, '..', 'outputs', 'basic_example_phase.png')
    )
    plot_phase_det(ts_A, ts_B, sv_A, sv_B, save_path=PHASE_PATH)

    # -----------------------------------------------------------------------
    # PART 6: Process vector testing — three regime examples
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("PART 6 — Process vector testing utility")
    print("=" * 60)

    x_ref = x0_conflict.copy()  # reference state for all three tests

    # Conflict dominant: strong extraction/propaganda/recruitment/stealing, weak peace
    # [u1,  u2,  u3,  u4,  c1,  c2,  c3,  c4,  f1,   f2  ]
    v_conflict = np.array([1.5, 1.0, 0.5, 0.8, 0.1, 0.2, 0.2, 0.1, 0.05, 0.02])

    # Balanced contested: moderate rates across all reactions
    v_balanced = np.array([1.0, 0.8, 0.3, 0.4, 0.3, 0.3, 0.3, 0.2, 0.05, 0.02])

    # Peacebuilding trajectory: strong trust-driven Rp growth + development, weak conflict
    v_peace    = np.array([0.5, 0.5, 0.1, 0.2, 0.6, 0.5, 0.4, 0.4, 0.05, 0.02])

    test_process_vector(v_conflict, x_ref, S, sp_names, label='conflict dominant')
    test_process_vector(v_balanced, x_ref, S, sp_names, label='balanced contested')
    test_process_vector(v_peace,    x_ref, S, sp_names, label='peacebuilding trajectory')

    # -----------------------------------------------------------------------
    # PART 7: Productive strength analysis
    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("PART 7 — Productive strength analysis")
    print("=" * 60)

    STRENGTH_PATH = os.path.normpath(
        os.path.join(_SCRIPT_DIR, '..', 'outputs', 'productive_strengths.png')
    )
    analyse_productive_strengths(flx_B, save_path=STRENGTH_PATH)
