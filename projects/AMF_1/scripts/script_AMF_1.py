# ========================================
# AMF – TRADE BALANCE MODEL v6.1 (Corregido: Flujos de entrada)
# ======================================== 
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import sys
from matplotlib.colors import TwoSlopeNorm

# Rutas absolutas para importar pyCOT (Ajusta si es necesario)
sys.path.append('/Users/yvanomarbalderamoreno/Downloads/pyCOT/src')
sys.path.append('/Users/yvanomarbalderamoreno/Downloads/pyCOT')

from pyCOT.io.functions import read_txt
from pyCOT.simulations.ode import simulation

VIS_DIR = 'projects/AMF_1/outputs_1'
os.makedirs(VIS_DIR, exist_ok=True)

# ========================================
# PLOTTING FUNCTIONS
# ========================================
def plot_amf_dynamics_new(time_series, title="AMF Dynamics", save_path=None):
    time       = time_series['Time'].values
    all_species = [c for c in time_series.columns if c != 'Time']
    groups = {
        "Biomass (Plant & Fungus)": ['Plant_active', 'Plant_limited', 'Fungus', 'Mycelium'],
        "Carbon Pools":             ['C_plant', 'C_fungus'],
        "Nitrogen":                 ['N_plant', 'N_fungus', 'N_avail'],
        "Phosphorus":               ['P_plant', 'P_fungus', 'P_avail'],
    }
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    for ax, (group, group_vars) in zip(axes, groups.items()):
        for var in [v for v in group_vars if v in all_species]:
            lw = 2.5 if var in ['Plant_active', 'Fungus', 'C_plant'] else 1.8
            ax.plot(time, time_series[var], linewidth=lw, label=var)
        ax.set_title(group, fontweight='bold', fontsize=13)
        ax.set_xlabel("Time"); ax.set_ylabel("Concentration")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_flux_dynamics_new(flux_df, title="Flux Dynamics", save_path=None):
    time = flux_df['Time'].values
    all_reactions = [c for c in flux_df.columns if c != 'Time']
    groups = {
        "Plant Core":                    ['R1', 'R2', 'R3', 'R4'],
        "Plant Direct Uptake":           ['R5', 'R6'],
        "Mycelium Uptake":               ['R9', 'R10'],
        "Recycling & Mortality":         ['R7', 'R8', 'R12', 'R17'],
        "Symbiotic Exchange & Growth":   ['R11', 'R13', 'R14', 'R15', 'R16'],
        "Inputs":                        ['R18', 'R19'],
    }
    nrows = 2; ncols = 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 10))
    axes = axes.flatten()
    for ax, (group, rxns) in zip(axes, groups.items()):
        for rxn in [r for r in rxns if r in all_reactions]:
            ls = '--' if rxn in ['R11', 'R13', 'R14', 'R17'] else '-'
            ax.plot(time, flux_df[rxn], linewidth=1.8, linestyle=ls, label=rxn)
        ax.set_title(group, fontweight='bold', fontsize=12)
        ax.set_xlabel("Time"); ax.set_ylabel("Flux")
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

# ========================================
# BIOLOGICAL BASE PARAMETERS (ESCALADOS)
# ======================================== 
def create_base_params():
    return [
        [2.5, 10.0],   # R1(mmk):  Fotosíntesis. <--- AUMENTADO (antes 2.5) Da más "salario" a la planta.
        [1e-6],        # R2(mak):  Activación. 
        [0.02],        # R3(mak):  Senescencia. 
        [8e-6],        # R4(mak):  Crecimiento planta. 
        [5e-4],        # R5(mak):  Absorción directa N raíces. 
        [1e-2],        # R6(mak):  Absorción directa P raíces. (Mantenemos tu ajuste fuerte)
        [5e-8],        # R7(mak):  Reciclaje Plant_active. 
        [5e-8],        # R8(mak):  Reciclaje Plant_limited. 
        [1.5, 30.0],   # R9(mmk):  Captación N micelio. 
        [1.2, 10.0],   # R10(mmk): Captación P micelio. 
        [5e-7],        # R11(mak): Crecimiento micelio. 
        [0.03],        # R12(mak): Mortalidad micelio. 
        [0.505, 50.0], # R13(mmk): Captación C hongo. <--- REDUCIDO (antes 1.8) Evita bancarrota de C.
        [0.8, 20.0],   # R14(mmk): Transferencia N hongo→planta. 
        [0.4, 20.0],   # R15(mmk): Transferencia P hongo→planta. <--- REDUCIDO (antes 1.5) Permite a la raíz ganar en P abundante.
        [6e-7],        # R16(mak): Crecimiento hongo. 
        [0.15],        # R17(mak): Mortalidad hongo. 
        [0.0],         # R18(mak): Entrada N (Sobrescrito por escenario)
        [0.0],         # R19(mak): Entrada P (Sobrescrito por escenario)
    ]

def create_rate_list():
    return [
        'mmk', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak', 'mak',
        'mmk', 'mmk', 'mak', 'mak', 'mmk', 'mmk', 'mmk', 'mak', 
        'mak', 'mak', 'mak'
    ]

# ========================================
# SCENARIO CONFIGURATION (CORREGIDA Y BALANCEADA)
# ========================================
_SPECIES_ORDER = [
    'Plant_active', 'C_plant', 'Plant_limited', 'N_plant', 'P_plant',
    'N_avail', 'P_avail', 'Mycelium', 'N_fungus', 'P_fungus',
    'C_fungus', 'Fungus',
]

def get_scenario_config(scenario: str, species_order=None):
    if species_order is None:
        species_order = _SPECIES_ORDER

    params = create_base_params()
    
    IDX_R18 = 17 # N input
    IDX_R19 = 18 # P input

    # Condiciones iniciales base de biomasa
    x0_dict = {
        'Plant_active': 5.0,  'C_plant': 0.5,
        'Plant_limited': 0.1, 'N_plant': 25, 'P_plant': 25,
        'Mycelium': 0.05,     'N_fungus': 20, 'P_fungus': 20,
        'C_fungus': 0.02,     'Fungus': 0.01,
        'N_avail': 30,        # Valor por defecto (Medio)
        'P_avail': 25         # Valor por defecto (Alto)
    }

    # Trade Balance Model.
    # Limitación (Limitation) = "Bajo" en la tabla.
    #   N limitado: < 20 mg/Kg (Usamos 15).
    #   P limitado: < 10 mg/Kg (Usamos 5).
    # Exceso/Lujo (Luxury) = "Alto" en la tabla.
    #   N en exceso: > 40 mg/Kg (Usamos 45).
    #   P en exceso: > 20 mg/Kg (Usamos 30).
    
    if scenario == 'I':
        x0_dict['N_avail'] = 5 # N limitado
        x0_dict['P_avail'] = 2  # P limitado
        params[IDX_R18][0] = 0.005   # Input de N Bajo
        params[IDX_R19][0] = 0.002   # Input de P Bajo
        params_name   = "I – N↓ P↓"
        scenario_name = "I – Mutualismo C-limitado"

    elif scenario == 'II':
        x0_dict['N_avail'] = 45 # N en exceso
        x0_dict['P_avail'] = 5  # P limitado
        params[IDX_R18][0] = 0.50   # Input de N Alto
        params[IDX_R19][0] = 0.02   # Input de P Bajo
        params_name   = "II – N↑ P↓"
        scenario_name = "II – Mutualismo fuerte"

    elif scenario == 'III':
        x0_dict['N_avail'] = 15 # N limitado
        x0_dict['P_avail'] = 30 # P en exceso
        params[IDX_R18][0] = 0.05   # Input de N Bajo
        params[IDX_R19][0] = 0.30   # Input de P Alto
        params_name   = "III – N↓ P↑"
        scenario_name = "III – Comensalismo"

    elif scenario == 'IV':
        x0_dict['N_avail'] = 45 # N en exceso
        x0_dict['P_avail'] = 30 # P en exceso
        params[IDX_R18][0] = 0.50   # Input de N Alto
        params[IDX_R19][0] = 0.30   # Input de P Alto
        params_name   = "IV – N↑ P↑"
        scenario_name = "IV – Parasitismo"

    else:
        raise ValueError("Invalid scenario")

    x0 = [x0_dict[s] for s in species_order]
    return params, x0, params_name, scenario_name

# ========================================
# INDICATORS
# ========================================
def calcular_indicadores_amf(time_series, flux_vector, t_span, steady_frac=0.2):
    t_ss_init = t_span[1] - steady_frac * (t_span[1] - t_span[0])
    time_arr  = flux_vector['Time'].values
    n_crit    = np.searchsorted(time_arr, t_ss_init)

    def ss(rxn): return float(flux_vector[rxn].iloc[n_crit:].mean())

    v1  = ss('R1'); v5 = ss('R5'); v6 = ss('R6')
    v13 = ss('R13'); v14 = ss('R14'); v15 = ss('R15')

    P_benefit = v15 / (v6 + v15 + 1e-10)
    ind = {'C_supply': v1, 'C_demand': v13, 'P_benefit': P_benefit}
    return ind, n_crit


# ========================================
# PLOT TRADE BALANCE MODEL
# ========================================
def plot_trade_balance_model(rn, save_path=None):
    rate_list = create_rate_list()
    n=1
    t_span = (0, 360 * n) 
    n_steps = 500 * n
    species_order = [s.name for s in rn.species()]
    scenarios = ['III', 'IV', 'I', 'II']
    results = {}

    colors = {'I': '#4A90D9', 'II': '#27AE60', 'III': '#F39C12', 'IV': '#E74C3C'}
    labels = {
        'I':   'I – Mut. C-limitado (N↓P↓)',
        'II':  'II – Mutualismo fuerte (N↑P↓)',
        'III': 'III – Comensalismo (N↓P↑)',
        'IV':  'IV – Parasitismo (N↑P↑)',
    }

    print("=" * 65)
    print("  AMF TRADE BALANCE MODEL (Parámetros Biológicos)")
    print("=" * 65)

    ts_dict = {}; flux_dict = {}

    for scen in scenarios:
        params, x0, params_name, scenario_name = get_scenario_config(scen, species_order)
        print(f"\n[{scen}] {scenario_name}")

        # Se ajustan ligeramente las tolerancias para evitar warnings numéricos
        time_series, flux_vector = simulation(
            rn, x0=x0, rate=rate_list, spec_vector=params,
            t_span=t_span, n_steps=n_steps, method='LSODA', rtol=1e-8, atol=1e-10
        )
        ts_dict[scen] = time_series
        flux_dict[scen] = flux_vector
        
        plot_amf_dynamics_new(time_series, title=f"Scenario {scenario_name}", save_path=os.path.join(VIS_DIR, '1_time_series', f'{scen}_ts.png'))
        plot_flux_dynamics_new(flux_vector, title=f"Scenario {scenario_name}", save_path=os.path.join(VIS_DIR, '2_flux_vector', f'{scen}_fd.png'))
        
        ind, _ = calcular_indicadores_amf(time_series, flux_vector, t_span)
        results[scen] = ind

    # ── Figure 1: Indicator temporal trajectories ──────────
    ind_colors = {'C_supply': '#2ECC71', 'C_demand': '#E67E22', 'P_benefit': '#3498DB'}
    fig1, axes1 = plt.subplots(2, 2, figsize=(16, 12))
    axes1 = axes1.flatten()

    for ax, scen in zip(axes1, scenarios):
        _, _, params_name, scenario_name = get_scenario_config(scen, species_order)
        t   = flux_dict[scen]['Time'].values
        v1  = flux_dict[scen]['R1'].values
        v13 = flux_dict[scen]['R13'].values
        v15 = flux_dict[scen]['R15'].values
        v6  = flux_dict[scen]['R6'].values
        pb  = v15 / (v15 + v6 + 1e-10)

        ax.plot(t, v1,  color=ind_colors['C_supply'],  lw=2.2, label='v1 – C supply')
        ax.plot(t, v13, color=ind_colors['C_demand'],  lw=2.2, label='v13 – C demand')
        ax.plot(t, pb,  color=ind_colors['P_benefit'], lw=2.2, label='P benefit = v15/(v15+v6)')

        t_ss = t_span[1] - 0.2 * (t_span[1] - t_span[0])
        ax.axvspan(t_ss, t_span[1], alpha=0.07, color='gray', label='SS Zone')
        ax.axhline(0, color='black', lw=0.8, ls='--', alpha=0.5)

        ax.set_title(f"{scenario_name}", fontweight='bold', fontsize=12, color=colors[scen])
        ax.set_xlabel("Time"); ax.set_ylabel("Flux / Ratio")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.25)

    fig1.suptitle("AMF Trade Balance Model\nTrayectorias de Indicadores Biológicos", fontsize=15, fontweight='bold')
    plt.tight_layout()
    if save_path:
        p1 = save_path + '_indicadores_ts.png'
        os.makedirs(os.path.dirname(p1), exist_ok=True)
        fig1.savefig(p1, dpi=300, bbox_inches='tight')
        print(f"Guardado: {p1}")
    plt.close(fig1)

    # ── Figure 2: Phase map 2×2 ────────────────────────────────
    layout = {'III': (0, 0), 'IV': (0, 1), 'I': (1, 0), 'II': (1, 1)}
    fig2, axes2 = plt.subplots(2, 2, figsize=(14, 10))
    bar_keys    = ['C_supply', 'C_demand', 'P_benefit']
    bar_xlabels = ['C supply', 'C demand', 'P benefit']
    bar_cols    = [ind_colors[k] for k in bar_keys]

    for scen, (row, col) in layout.items():
        ax   = axes2[row][col]
        vals = [results[scen][k] for k in bar_keys]
        bars = ax.bar(bar_xlabels, vals, color=bar_cols, alpha=0.85, edgecolor='white', linewidth=1.2, width=0.55)
        for bar, val in zip(bars, vals):
            ypos = val + 0.003 if val >= 0 else val - 0.006
            ax.text(bar.get_x() + bar.get_width() / 2, ypos, f"{val:.4f}", ha='center', va='bottom' if val >= 0 else 'top', fontsize=9, fontweight='bold')
        ax.axhline(0, color='black', lw=1.0, ls='--', alpha=0.6)
        ax.set_ylabel("Estado Estacionario", fontsize=9)
        ax.set_title(labels[scen], fontweight='bold', fontsize=11, color=colors[scen])
        ax.grid(True, axis='y', alpha=0.25)
        ax.set_ylim(0, max(1.5, max(vals) * 1.05))

    fig2.text(0.5, 0.01, "Eje N → Limitación (izq) | Exceso (der)", ha='center', fontsize=11, style='italic', color='#555')
    fig2.text(0.01, 0.5, "Eje P → Limitación (abajo) | Exceso (arriba)", va='center', rotation='vertical', fontsize=11, style='italic', color='#555')
    fig2.suptitle("AMF Trade Balance Model – Mapa de Fase N × P", fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0.04, 0.04, 1, 1])
    if save_path:
        p2 = save_path + '_mapa_fase.png'
        fig2.savefig(p2, dpi=300, bbox_inches='tight')
        print(f"Guardado: {p2}")
    plt.close(fig2)

    return results

# ========================================
# PLOT MYCORRHIZAL FUNCTION
# ========================================
def plot_mycorrhizal_function(rn, save_path=None):
    rate_list     = create_rate_list()
    n=1
    t_span        = (0, 360 * n)
    n_steps       = 500 * n
    species_order = [s.name for s in rn.species()]
    scenarios     = ['I', 'II', 'III', 'IV']
    steady_frac   = 0.2

    colors = {'I': '#4A90D9', 'II': '#27AE60', 'III': '#F39C12', 'IV': '#E74C3C'}
    titles = {
        'I':   'I – Mut. C-lim\n(N↓P↓)',
        'II':  'II – Mutualismo\n(N↑P↓)',
        'III': 'III – Comens.\n(N↓P↑)',
        'IV':  'IV – Parasit.\n(N↑P↑)',
    }
    col_v1 = '#2ECC71'; col_v13 = '#E67E22'; col_pr = '#3498DB'

    print("\n" + "=" * 65)
    print("  AMF MYCORRHIZAL FUNCTION – MODULACIÓN POR N Y P")
    print("=" * 65)

    mf_results = {}

    for scen in scenarios:
        params, x0, _, _ = get_scenario_config(scen, species_order)
        _, flux_vector = simulation(
            rn, x0=x0, rate=rate_list, spec_vector=params,
            t_span=t_span, n_steps=n_steps, method='LSODA', rtol=1e-8, atol=1e-10
        )

        t_ss  = t_span[1] - steady_frac * (t_span[1] - t_span[0])
        n_c   = np.searchsorted(flux_vector['Time'].values, t_ss)
        ss    = lambda r: float(flux_vector[r].iloc[n_c:].mean())

        v1  = ss('R1'); v13 = ss('R13')
        v15 = ss('R15'); v6 = ss('R6')
        pb  = v15 / (v15 + v6 + 1e-10)

        mf_results[scen] = {'v1': v1, 'v13': v13, 'v15': v15, 'v6': v6, 'P_benefit': pb}

    # ── Heatmap ───────────────────────────────────
    layout = {'III': (1, 0), 'IV': (0, 0), 'I': (1, 1), 'II': (0, 1)}
    metrics    = ['v1', 'v13', 'P_benefit']
    met_labels = ['C supply\n(v1 = R1)', 'C demand\n(v13 = R13)', 'P exchange\n(v15/(v15+v6))']
    met_colors = [col_v1, col_v13, col_pr]

    matrices = {m: np.zeros((2, 2)) for m in metrics}
    for scen, (row, col) in layout.items():
        for m in metrics:
            matrices[m][row, col] = mf_results[scen][m]

    fig2, axes2 = plt.subplots(1, 3, figsize=(17, 6))
    scen_in_cell = {v: k for k, v in layout.items()}

    for ax, m, mlabel, mcol in zip(axes2, metrics, met_labels, met_colors):
        mat  = matrices[m]
        vmin, vmax = mat.min(), mat.max()
        vc   = 0.0 if vmin < 0 < vmax else (vmin + vmax) / 2
        try:
            norm = TwoSlopeNorm(vmin=vmin-1e-6, vcenter=vc, vmax=vmax+1e-6)
            cmap = 'RdYlGn' if m != 'v13' else 'RdYlGn_r'
            im = ax.imshow(mat, cmap=cmap, norm=norm, aspect='auto')
        except Exception:
            im = ax.imshow(mat, cmap='YlOrRd', aspect='auto')

        plt.colorbar(im, ax=ax, shrink=0.82)
        for r in range(2):
            for c in range(2):
                val = mat[r, c]
                ax.text(c, r, f"{val:.3f}", ha='center', va='center',
                        fontsize=13, fontweight='bold',
                        color='white' if abs(val) > 0.4*max(abs(vmin), abs(vmax)) else 'black')
                ax.text(c, r-0.38, f"Esc. {scen_in_cell.get((r,c),'?')}",
                        ha='center', va='center', fontsize=8, color='white', style='italic')

        ax.set_xticks([0, 1]); ax.set_xticklabels(['P lux', 'P lim'])
        ax.set_yticks([0, 1]); ax.set_yticklabels(['N lux', 'N lim'])
        ax.set_title(mlabel, fontweight='bold', fontsize=11, color=mcol, pad=10)

        for (row, col), scen in scen_in_cell.items():
            ax.add_patch(plt.Rectangle((col-0.5, row-0.5), 1, 1, linewidth=3, edgecolor=colors[scen], facecolor='none', zorder=5))

    legend_patches = [mpatches.Patch(color=colors[e], label=titles[e].replace('\n', ' ')) for e in scenarios]
    fig2.legend(handles=legend_patches, loc='lower center', ncol=4, fontsize=9, bbox_to_anchor=(0.5, -0.06), title='Escenarios', title_fontsize=9)
    fig2.suptitle("AMF Mycorrhizal Function – Mapa de Fase N × P\nC supply, C demand, P exchange ratio (estado estacionario)", fontsize=13, fontweight='bold')
    plt.tight_layout(rect=[0, 0.06, 1, 1])
    
    if save_path:
        p2 = save_path + '_mycorrhizal_function_heatmap.png'
        fig2.savefig(p2, dpi=300, bbox_inches='tight')
        print(f"Guardado: {p2}")
    plt.close(fig2)

    # ── Tabla resumen ──────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("  RESUMEN – FUNCIÓN MICORRÍCICA")
    print("=" * 72)
    print(f"  {'Escenario':<32} {'v1':>7} {'v13':>7} {'v15':>7} {'v6':>7} {'P_ben':>7}")
    print(f"  {'-'*62}")
    names = {
        'I':   'I   – Mut. C-limitado (N↓P↓)',
        'II':  'II  – Mutualismo fuerte (N↑P↓)',
        'III': 'III – Comensalismo    (N↓P↑)',
        'IV':  'IV  – Parasitismo     (N↑P↑)',
    }
    for scen in scenarios:
        r = mf_results[scen]
        print(f"  {names[scen]:<32} {r['v1']:>+7.4f} {r['v13']:>+7.4f} {r['v15']:>+7.4f} {r['v6']:>+7.4f} {r['P_benefit']:>+7.4f}")
    
    return mf_results

# ========================================
# MAIN
# ========================================
if __name__ == '__main__':
    FILE_PATH = 'data/Ecological_models/AMF_1.txt'
    rn = read_txt(FILE_PATH)
    
    print(f"Especies:    {[s.name for s in rn.species()]}")
    print(f"Reacciones:  {[r.name() for r in rn.reactions()]}")
    
    save_base = os.path.join(VIS_DIR, '3_trade_balance', 'trade_balance_v6')
    os.makedirs(os.path.dirname(save_base), exist_ok=True)
    
    results = plot_trade_balance_model(rn, save_path=save_base)
    mf_results = plot_mycorrhizal_function(rn, save_path=save_base)
    
    print("\nSimulación completada con éxito.")