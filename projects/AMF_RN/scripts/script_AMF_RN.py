# ========================================
# AMF – INDICADORES FUNCIONALES (TRADE BALANCE MODEL)
# ========================================

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from matplotlib.lines import Line2D 
import mplcursors

VIS_DIR = 'projects/AMF_RN/outputs'
os.makedirs(VIS_DIR, exist_ok=True)

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'src'))
from pyCOT.io.functions import read_txt
from pyCOT.simulations.ode import simulation

# ========================================
# 1. LOAD NETWORK
# ========================================
file_path = 'data/Ecological_models/AMF_RN.txt'
rn = read_txt(file_path)

species_amf = [specie.name for specie in rn.species()]
print(f"Especies en la red AMF: {species_amf}")

reaction_names = [reaction.name() for reaction in rn.reactions()]
print(f"Reacciones en la red AMF: {reaction_names}")

# ========================================
# 2. PLOTTING FUNCTIONS
# ========================================
def plot_amf_dynamics_new(time_series, title="AMF Dynamics", save_path=None):
    """Grafica todas las especies de la red AMF, agrupadas por tipo."""

    time = time_series['Time'].values
    # Todas las especies disponibles (excluye columna Time)
    all_species = [c for c in time_series.columns if c != 'Time']

    # Agrupación por tipo
    grupos = {
        "Biomass (Plant & Fungus)": ['Plant_active', 'Plant_limited', 'Fungus', 'Mycelium'],
        "Carbon Pools":             ['C_plant', 'C_fungus'],
        "Nitrogen":                 ['N_plant', 'N_fungus', 'N_soil'],
        "Phosphorus":               ['P_plant', 'P_fungus', 'P_avail', 'P_no_avail'],
    }

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for ax, (grupo, vars_grupo) in zip(axes, grupos.items()):
        presentes = [v for v in vars_grupo if v in all_species]
        # Si hay especies no asignadas a ningún grupo, las añadimos al primer panel
        for var in presentes:
            lw = 2.5 if var in ['Plant_active', 'Fungus', 'C_plant'] else 1.8
            ax.plot(time, time_series[var], linewidth=lw, label=var)
        ax.set_title(grupo, fontweight='bold', fontsize=13)
        ax.set_xlabel("Time", fontsize=11)
        ax.set_ylabel("Concentration", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # Especies no cubiertas por ningún grupo → avisamos en título
    cubiertas = {v for vs in grupos.values() for v in vs}
    no_cubiertas = [v for v in all_species if v not in cubiertas]
    if no_cubiertas:
        print(f"  ⚠ Especies no graficadas: {no_cubiertas}")

    fig.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    plt.close()
    return fig, axes


def plot_flux_dynamics_new(flux_df, title="Flux Dynamics", save_path=None):
    """Grafica todas las reacciones de la red AMF, agrupadas por proceso."""

    time = flux_df['Time'].values
    # Todas las reacciones disponibles (excluye columna Time)
    all_reactions = [c for c in flux_df.columns if c != 'Time']

    # Agrupación por proceso biológico
    grupos = {
        "Plant Core":           ['R1', 'R2', 'R3'], 
        "Mycelium Uptake":      ['R7', 'R8', 'R9'],
        "Plant & Mycelium & Fungus Growth":    ['R4', 'R10', 'R13'],
        "Plant & Mycelium & Fungus Recycling": ['R5', 'R6', 'R11', 'R14'],
        "Symbiotic Exchange":       ['R12', 'R15', 'R16'],
        "P Sequestration & Inputs": ['R17', 'R18', 'R19'],
    }

    n_grupos = len(grupos)
    ncols = 3
    nrows = (n_grupos + ncols - 1) // ncols  # techo de la división

    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 5 * nrows))
    axes = axes.flatten()

    for ax, (grupo, rxns) in zip(axes, grupos.items()):
        presentes = [r for r in rxns if r in all_reactions]
        for rxn in presentes:
            lw = 2.5 if rxn in ['R1'] else 1.8
            ls = '--' if rxn in ['R10', 'R13', 'R11', 'R14', 'R17'] else '-'
            ax.plot(time, flux_df[rxn], linewidth=lw, linestyle=ls, label=rxn)
        ax.set_title(grupo, fontweight='bold', fontsize=12)
        ax.set_xlabel("Time", fontsize=10)
        ax.set_ylabel("Flux", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # Ocultar ejes sobrantes
    for ax in axes[n_grupos:]:
        ax.set_visible(False)

    # Reacciones no cubiertas → aviso
    cubiertas = {r for rs in grupos.values() for r in rs}
    no_cubiertas = [r for r in all_reactions if r not in cubiertas]
    if no_cubiertas:
        print(f"  ⚠ Reacciones no graficadas: {no_cubiertas}")

    fig.suptitle(title, fontsize=16, fontweight='bold')
    plt.tight_layout()

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    plt.close()
    return fig, axes

# ========================================
# 3. PARAMETERS AND RATES
# ========================================
def create_base_params():
    """
    Parámetros base calibrados por balance de flujos en estado estacionario.
    Concentraciones objetivo en estado estacionario: ~0.3–0.5 para todas
    las especies principales. 
    """
    base_params = [
        [0.9],        # R1:  Fotosíntesis              (MAK)
        [0.5],        # R2:  Activación por nutrientes (MAK)
        [0.03],       # R3:  Limitación / Senescencia  (MAK)
        [0.3],        # R4:  Crecimiento planta        (MAK)    
        [0.02],       # R5:  Mortalidad Plant_active   (MAK)
        [0.02],       # R6:  Mortalidad Plant_limited  (MAK)
        [0.25, 0.1],  # R7:  Captación N micelio            (MMK) 
        [0.6,  0.1],  # R8:  Captación P micelio            (MMK)
        [0.3,  0.1],  # R9:  Movilización P_secuestrado     (MMK)
        [3.0],        # R10: Crecimiento micelio       (MAK)    
        [0.03],       # R11: Mortalidad micelio        (MAK)     
        [0.3,  0.15], # R12: Captación C del hongo          (MMK)
        [3.0],        # R13: Crecimiento hongo         (MAK)       
        [0.02],       # R14: Mortalidad hongo          (MAK)
        [0.5,  0.1],  # R15: Transferencia N hongo→planta   (MMK)
        [0.5,  0.1],  # R16: Transferencia P hongo→planta   (MMK)
        [0.05],       # R17: Secuestro de P            (MAK)
        [0.15],       # R18: Entrada N al suelo        (MAK)      
        [0.08],       # R19: Entrada P al suelo        (MAK)      
    ]
    return base_params

def get_param_version(version="base", base_params=None):
    """
    Devuelve una versión modificada de los parámetros base para simular diferentes escenarios de simbiosis AMF: 
    - "v1": Mutualismo limitado por C (N↓ P↓ C↓)
    - "v2": Mutualismo fuerte (N↑ P↓ C↑)
    - "v3": Comensalismo (N↓ P↑ C↓)
    - "v4": Parasitismo (N↑ P↑ C↑)
    Cada versión ajusta principalmente los Vmax de las reacciones clave (R7, R8, R12, R15, R16) para 
    reflejar las condiciones de N, P y C en cada escenario, 
    manteniendo las constantes de crecimiento (R4, R10, R13) iguales para 
    aislar el efecto de la disponibilidad de nutrientes y carbono.
    """
    if base_params is None:
        base_params = create_base_params()

    params = [p.copy() if isinstance(p, list) else p for p in base_params]

    if version == "base":
        return params

    elif version == "v1":
        # ── Escenario I: N↓ P↓ C↓ — Mutualismo limitado por C ──────────
        # Condición: suelo pobre en N y P; fotosíntesis baja → C escaso.
        # Hongo devuelve N y P a la planta. 
        # Beneficio mutuo pero limitado por la escasez de C que la planta puede ceder. 
        params[0][0]  = 0.6    # R1:  Fotosíntesis moderada-baja → C↓
        params[6][0]  = 0.15   # R7:  Vmax N_uptake bajo → N_soil↓
        params[7][0]  = 0.15   # R8:  Vmax P_uptake bajo → P_avail↓
        params[8][0]  = 0.2    # R9:  Vmax P_seq bajo
        params[11][0] = 0.3    # R12: Vmax C_uptake moderado
        params[14][0] = 2.0    # R15: Vmax N_transfer ALTO → hongo devuelve N
        params[15][0] = 2.0    # R16: Vmax P_transfer ALTO → hongo devuelve P
        params[17][0] = 0.05   # R18: Entrada N BAJA → N_soil↓
        params[18][0] = 0.02   # R19: Entrada P BAJA → P_avail↓
        return params

    elif version == "v2":
        # ── Escenario II: N↑ P↓ C↑ — Mutualismo fuerte ─────────────────
        # Condición: N abundante en suelo; P escaso; fotosíntesis alta → C↑.
        # Hongo aporta P (limitante) a cambio de C (abundante). 
        # Máximo beneficio mutuo: planta crece bien, hongo también. 
        params[0][0]  = 1.5    # R1:  Fotosíntesis alta → C↑
        params[1][0]  = 1.2    # R2:  Activación alta
        params[6][0]  = 1.5    # R7:  Vmax N_uptake alto (alta demanda)
        params[7][0]  = 0.2    # R8:  Vmax P_uptake bajo → P↓ en hongo
        params[8][0]  = 0.4    # R9:  Vmax P_seq moderado
        params[11][0] = 0.75   # R12: Vmax C_uptake alto
        params[14][0] = 2.0    # R15: Vmax N_transfer moderado
        params[15][0] = 1.0    # R16: Vmax P_transfer ALTO → hongo devuelve P
        params[17][0] = 0.35   # R18: Entrada N MUY ALTA → N_soil↑
        params[18][0] = 0.02   # R19: Entrada P BAJA → P_avail↓
        return params

    elif version == "v3":
        # ── Escenario III: N↓ P↑ C↓ — Comensalismo ─────────────────────
        # Condición: N escaso; P abundante; fotosíntesis moderada → C↓.
        # Hongo toma C de la planta pero NO devuelve N ni P en cantidad significativa. 
        # Planta no gana nada; hongo se beneficia. 
        params[0][0]  = 0.7    # R1:  Fotosíntesis moderada → C↓
        params[6][0]  = 0.15   # R7:  Vmax N_uptake bajo → N_soil↓
        params[7][0]  = 1.5    # R8:  Vmax P_uptake alto → P↑ en hongo
        params[8][0]  = 0.2    # R9:  Vmax P_seq bajo
        params[11][0] = 0.8    # R12: Vmax C_uptake moderado (hongo toma C)
        params[14][0] = 1.0    # R15: Vmax N_transfer MUY BAJO → planta no gana N
        params[15][0] = 2.0    # R16: Vmax P_transfer MUY BAJO → planta no gana P
        params[17][0] = 0.5    # R18: Entrada N BAJA → N_soil↓
        params[18][0] = 0.2    # R19: Entrada P MUY ALTA → P_avail↑
        return params

    elif version == "v4":
        # ── Escenario IV: N↑ P↑ C↑ — Parasitismo ───────────────────────
        # Condición: N y P abundantes; fotosíntesis alta → C↑. 
        # Hongo toma C de la planta pero NO devuelve N ni P.
        # Planta no gana nada; hongo se beneficia mucho.
        params[0][0]  = 1.5    # R1:  Fotosíntesis alta → C↑
        params[1][0]  = 1.0    # R2:  Activación moderada
        params[6][0]  = 1.5    # R7:  Vmax N_uptake alto → N↑
        params[7][0]  = 1.5    # R8:  Vmax P_uptake alto → P↑
        params[8][0]  = 1.0    # R9:  Vmax P_seq alto
        params[11][0] = 2.0    # R12: Vmax C_uptake MUY ALTO (2.0) (parasítico → La planta pierde C sin recibir nutrientes y el hongo extrae máximo C posible)
        params[14][0] = 9.0    # R15: Vmax N_transfer → planta no gana N
        params[15][0] = 8.0    # R16: Vmax P_transfer → planta no gana P
        params[17][0] = 0.15   # R18: Entrada N MUY ALTA → N_soil↑
        params[18][0] = 0.3    # R19: Entrada P MUY ALTA → P_avail↑
        return params

    else:
        raise ValueError(f"Versión '{version}' no válida.")

def create_rate_list():
    return [
        'mak', 'mak', 'mak', 'mak', 'mak', 'mak',
        'mmk', 'mmk', 'mmk',
        'mak', 'mak',
        'mmk',
        'mak', 'mak',
        'mmk', 'mmk',
        'mak', 'mak', 'mak',
    ]

def create_initial_conditions(escenario=1):
    if escenario == 1:
        x0_dict = {
            'Plant_active':  1.5, 'C_plant': 0.8,
            'Plant_limited': 0.1, 'N_plant': 0.2, 'P_plant': 0.2,
            'N_soil': 0.15,   # bajo → N↓
            'P_avail': 0.10,  # bajo → P↓
            'Mycelium': 0.4,  'N_fungus': 0.5, 'P_fungus': 0.4,
            'P_no_avail': 0.3,'C_fungus': 0.5, 'Fungus': 0.4,
        }
    elif escenario == 2:
        x0_dict = {
            'Plant_active':  1.5,
            'C_plant': 0.8,
            'Plant_limited': 0.1, 'N_plant': 0.4, 'P_plant': 0.2,
            'N_soil': 0.50,   # alto → N↑
            'P_avail': 0.10,  # bajo → P↓
            'Mycelium': 0.4,  'N_fungus': 0.5, 'P_fungus': 0.4,
            'P_no_avail': 0.3,'C_fungus': 0.5, 'Fungus': 0.4,
        }
    elif escenario == 3:
        x0_dict = {
            'Plant_active':  0.5, 'C_plant': 0.2,
            'Plant_limited': 0.3, 'N_plant': 0.4, 'P_plant': 0.8,
            'N_soil': 0.10,   # bajo → N↓
            'P_avail': 0.50,  # alto → P↑
            'Mycelium': 0.3,  'N_fungus': 0.2, 'P_fungus': 0.5,
            'P_no_avail': 0.4,'C_fungus': 0.3, 'Fungus': 0.3,
        }
    elif escenario == 4:
        x0_dict = {
            'Plant_active':  0.6, 'C_plant': 0.8,
            'Plant_limited': 0.6, 'N_plant': 0.2, 'P_plant': 0.8,
            'N_soil': 0.90,   # alto → N↑
            'P_avail': 0.50,  # alto → P↑
            'Mycelium': 0.6,  'N_fungus': 0.7, 'P_fungus': 0.7,
            'P_no_avail': 0.5,'C_fungus': 0.6, 'Fungus': 0.4,
        }
    else:
        raise ValueError(f"Escenario {escenario} no válido.")

    return [x0_dict[name] for name in species_amf]

def get_scenario_config(escenario=1):
    mapa_params = {1: "v1", 2: "v2", 3: "v3", 4: "v4"}
    nombres_escenario = {
        1: "Mutualismo limitado por C (N↓ P↓ C↓)",
        2: "Mutualismo fuerte (N↑ P↓ C↑)",
        3: "Comensalismo (N↓ P↑ C↓)",
        4: "Parasitismo (N↑ P↑ C↑)",
    }
    version_params = mapa_params.get(escenario, "base")
    nombre_params  = f"v{escenario}" if escenario in mapa_params else "base"
    nombre_escenario = nombres_escenario.get(escenario, f"Escenario {escenario}")

    base_params = create_base_params()
    params = get_param_version(version_params, base_params)
    x0 = create_initial_conditions(escenario)
    return params, x0, nombre_params, nombre_escenario

# ========================================
# 3. FUNCIÓN AUXILIAR
# ========================================
def _extract_time_and_data(arr, expected_cols):
    if arr.shape[1] == expected_cols + 1:
        return arr[:, 0], arr[:, 1:]
    elif arr.shape[1] == expected_cols:
        return None, arr
    else:
        raise ValueError(
            f"Array tiene {arr.shape[1]} columnas; "
            f"se esperaban {expected_cols} o {expected_cols + 1}."
        )


# ========================================
# 4. CALCULAR INDICADORES 
# ========================================
def calcular_indicadores_amf(time_series, flux_vector, species_names,
                              n_reactions=19,
                              t_span=(0, 360*10),
                              n_steps=500*3):   
    """
    Calcula indicadores funcionales de la simbiosis AMF.
    """

    ts_raw = time_series.values if isinstance(time_series, pd.DataFrame) else np.asarray(time_series)
    fv_raw = flux_vector.values  if isinstance(flux_vector,  pd.DataFrame) else np.asarray(flux_vector)

    n_species = len(species_names)
    tiempo_col, ts = _extract_time_and_data(ts_raw, n_species)
    _,          fv = _extract_time_and_data(fv_raw, n_reactions)

    tiempo_real = tiempo_col if tiempo_col is not None else \
                  np.linspace(t_span[0], t_span[1], n_steps)

    ts = np.clip(ts, 0.0, None)   # concentraciones ≥ 0, sin techo
    fv = np.clip(fv, 0.0, None)   # flujos ≥ 0

    sd = {name: idx for idx, name in enumerate(species_names)}

    # ------------------------------------------------------------------
    # 1. Extracción de concentraciones
    # ------------------------------------------------------------------
    plant_active  = ts[:, sd['Plant_active']]
    plant_limited = ts[:, sd['Plant_limited']]
    fungus        = ts[:, sd['Fungus']]
    mycelium      = ts[:, sd['Mycelium']]
    c_plant       = ts[:, sd['C_plant']]
    n_plant       = ts[:, sd['N_plant']]
    p_plant       = ts[:, sd['P_plant']]
    n_fungus      = ts[:, sd['N_fungus']]
    p_fungus      = ts[:, sd['P_fungus']]
    n_soil        = ts[:, sd['N_soil']]
    p_avail       = ts[:, sd['P_avail']]
    p_no_avail    = ts[:, sd['P_no_avail']]

    biomasa_planta = plant_active + plant_limited
    biomasa_hongo  = fungus + mycelium

    # ------------------------------------------------------------------
    # 2. Tasas de crecimiento
    # ------------------------------------------------------------------
    crecimiento_planta = np.gradient(biomasa_planta, tiempo_real)
    crecimiento_hongo  = np.gradient(biomasa_hongo,  tiempo_real)

    # ------------------------------------------------------------------
    # 3. Flujos clave
    # ------------------------------------------------------------------
    prod_c_planta = fv[:, 0]   # R1:  fotosíntesis
    adq_c_hongo   = fv[:, 11]  # R12: adquisición C por el hongo

    intercambio_c = fv[:, 11]  # R12: C planta → hongo
    intercambio_n = fv[:, 14]  # R15: N hongo  → planta
    intercambio_p = fv[:, 15]  # R16: P hongo  → planta

    demanda_c_hongo = fv[:, 10] + fv[:, 12]    # R11 + R13
    oferta_c_planta = prod_c_planta - fv[:, 3]  # R1  - R4

    p_total_suelo = p_avail + p_no_avail

    # ------------------------------------------------------------------
    # 4. Umbrales fijos  
    N_UMBRAL_BAJO  = 0.12
    P_UMBRAL_BAJO  = 0.10
    N_UMBRAL_ALTO  = 0.40
    P_UMBRAL_ALTO  = 0.35
    N_PUNTO_CORTE  = 0.30
    P_PUNTO_CORTE  = 0.25

    media_n = float(n_soil[n_soil > 0].mean()) if np.any(n_soil > 0) else 0.0
    media_p = float(p_avail[p_avail > 0].mean()) if np.any(p_avail > 0) else 0.0

    n_critico = N_UMBRAL_BAJO # N_UMBRAL_ALTO if media_n > N_PUNTO_CORTE else N_UMBRAL_BAJO # Si el suelo es rico en N, el umbral para considerar N↓ es más alto (0.40), y si es pobre, el umbral es más bajo (0.20).
    p_critico = P_UMBRAL_BAJO #P_UMBRAL_ALTO if media_p > P_PUNTO_CORTE else P_UMBRAL_BAJO # Si el suelo es rico en P, el umbral para considerar P↓ es más alto (0.35), y si es pobre, el umbral es más bajo (0.15).

    print(f"\nUmbrales fijos calibrados v6 (selección automática):")
    print(f"  media N_soil  = {media_n:.4f}  → n_critico = {n_critico:.4f}"
          f"  (min={n_soil.min():.4f}, max={n_soil.max():.4f})")
    print(f"  media P_avail = {media_p:.4f}  → p_critico = {p_critico:.4f}"
          f"  (min={p_avail.min():.4f}, max={p_avail.max():.4f})")

    # Variables binarias de limitación 
    n_lim = (n_soil  <= n_critico).astype(float)  # 1 → N limitado (N↓)
    p_lim = (p_avail <= p_critico).astype(float)  # 1 → P limitado (P↓)
    c_lim = (oferta_c_planta <= demanda_c_hongo).astype(float)  # 1 → C limitado (C↓)

    # Cuatro cuadrantes del Trade Balance Model
    #   Escenario I:   N↓ P↓ C↓  →  n_lim=1, p_lim=1, c_lim=1
    #   Escenario II:  N↑ P↓ C↑  →  n_lim=0, p_lim=1, c_lim=0
    #   Escenario III: N↓ P↑ C↓  →  n_lim=1, p_lim=0, c_lim=1
    #   Escenario IV:  N↑ P↑ C↑  →  n_lim=0, p_lim=0, c_lim=0
    #
    # Los 4 casos restantes de (n_lim, p_lim, c_lim) se asignan al
    # cuadrante más cercano en distancia L1 para cobertura total.
    condicion_I   = (n_lim == 1) & (p_lim == 1) & (c_lim == 1)
    condicion_II  = (n_lim == 0) & (p_lim == 1) & (c_lim == 0)
    condicion_III = (n_lim == 1) & (p_lim == 0) & (c_lim == 1)
    condicion_IV  = (n_lim == 0) & (p_lim == 0) & (c_lim == 0)

    # Casos mixtos (distancia L1 mínima a cada cuadrante canónico)
    # Cuadrante I   canónico: (1,1,1)
    # Cuadrante II  canónico: (0,1,0)
    # Cuadrante III canónico: (1,0,1)
    # Cuadrante IV  canónico: (0,0,0)
    estado = np.column_stack([n_lim, p_lim, c_lim])  # shape (n, 3)
    canonicos = np.array([[1,1,1],[0,1,0],[1,0,1],[0,0,0]])  # 4 cuadrantes

    dist = np.abs(estado[:, None, :] - canonicos[None, :, :]).sum(axis=2)  # (n, 4)
    idx_cercano = np.argmin(dist, axis=1)  # 0=I, 1=II, 2=III, 3=IV

    # Para instantes que no caen en ningún cuadrante canónico, asignar al más cercano
    sin_cuadrante = ~(condicion_I | condicion_II | condicion_III | condicion_IV)
    condicion_I   = condicion_I   | (sin_cuadrante & (idx_cercano == 0))
    condicion_II  = condicion_II  | (sin_cuadrante & (idx_cercano == 1))
    condicion_III = condicion_III | (sin_cuadrante & (idx_cercano == 2))
    condicion_IV  = condicion_IV  | (sin_cuadrante & (idx_cercano == 3))

    # Verificación: cobertura total
    cobertura = condicion_I | condicion_II | condicion_III | condicion_IV
    assert cobertura.all(), "ERROR: hay instantes sin clasificar"

    # Intensidades mecanísticas 
    eps = 1e-10
    max_c = intercambio_c.max() + eps
    max_p = intercambio_p.max() + eps
    max_n = intercambio_n.max() + eps

    km_n_fotos = np.median(n_plant) + eps
    factor_n_fotos = n_plant / (n_plant + km_n_fotos)

    factor_c_disponible = np.clip(
        oferta_c_planta / (demanda_c_hongo + eps), 0.0, 1.0
    )

    km_n_hongo = np.median(n_fungus) + eps
    factor_n_hongo_limitado = 1.0 - n_fungus / (n_fungus + km_n_hongo)
    factor_n_hongo_libre    =       n_fungus / (n_fungus + km_n_hongo)

    factor_sin_ganancia_p = 1.0 - intercambio_p / max_p

    # Escenario I: planta cede C limitado; hongo transfiere N y P
    intensidad_I = condicion_I * factor_n_fotos * factor_c_disponible * \
                   (intercambio_c / max_c)

    # Escenario II: alta fotosíntesis; hongo transfiere P
    intensidad_II = condicion_II * factor_n_fotos * \
                    (intercambio_c / max_c) * (intercambio_p / max_p)

    # Escenario III: hongo toma C pero no devuelve nada (N-limitado)
    intensidad_III = condicion_III * (intercambio_c / max_c) * \
                     factor_n_hongo_limitado

    # Escenario IV: hongo sin restricción; planta sin retorno de P
    intensidad_IV = condicion_IV * (intercambio_c / max_c) * \
                    factor_n_hongo_libre * factor_sin_ganancia_p

    # Clasificación final
    etiquetas = [
        'I: Mutualismo C-limitado',
        'II: Mutualismo fuerte',
        'III: Comensalismo',
        'IV: Parasitismo',
    ]

    # Construir matriz de intensidades con piso mínimo según cuadrante asignado 
    int_matrix = np.column_stack([
        intensidad_I,
        intensidad_II,
        intensidad_III,
        intensidad_IV,
    ])

    # Añadir piso pequeño para el cuadrante de pertenencia
    pertenencia = np.column_stack([
        condicion_I, condicion_II, condicion_III, condicion_IV
    ]).astype(float)
    int_matrix += pertenencia * 1e-12  # piso ínfimo pero no cero

    idx_dominante = np.argmax(int_matrix, axis=1)
    escenario_dominante = np.array([etiquetas[i] for i in idx_dominante])

    # 6. Eficiencias e índices ROI 
    roi_planta = np.clip((intercambio_n + intercambio_p) / (intercambio_c + eps), 0.0, 50.0)
    roi_hongo = np.clip(intercambio_c / ((intercambio_n + intercambio_p) + eps), 0.0, 50.0)
    roi_N_C = np.clip(intercambio_n / (intercambio_c + eps), 0.0, 50.0)
    roi_P_C = np.clip(intercambio_p / (intercambio_c + eps), 0.0, 50.0)

    # ------------------------------------------------------------------
    # 7. DataFrame de salida
    # ------------------------------------------------------------------
    indicadores = pd.DataFrame({
        'Tiempo': tiempo_real,

        # Biomasas
        'Biomasa_Planta':  biomasa_planta,
        'Biomasa_Hongo':   biomasa_hongo,
        'Plant_Active':    plant_active,
        'Plant_Limited':   plant_limited,
        'Fungus':          fungus,
        'Mycelium':        mycelium,

        # Crecimiento
        'Crecimiento_Planta': crecimiento_planta,
        'Crecimiento_Hongo':  crecimiento_hongo,

        # Flujos de C
        'Prod_C_Planta': prod_c_planta,
        'Adq_C_Hongo':   adq_c_hongo,

        # Intercambios simbióticos
        'Intercambio_C': intercambio_c,
        'Intercambio_N': intercambio_n,
        'Intercambio_P': intercambio_p,

        # ── Flujos por categoría para fila 1 de visualización ──────────
        # Crecimiento: v4 planta, v10 micelio, v13 hongo
        'Flujo_Crec_Planta':   fv[:, 3],   # v4:  crecimiento Plant_active
        'Flujo_Crec_Micelio':  fv[:, 9],   # v10: crecimiento Mycelium
        'Flujo_Crec_Hongo':    fv[:, 12],  # v13: crecimiento Fungus

        # Reciclaje/mortalidad: v5, v6 (planta), v11 (micelio), v14 (hongo)
        'Flujo_Mort_PlantA':   fv[:, 4],   # v5:  mortalidad Plant_active
        'Flujo_Mort_PlantL':   fv[:, 5],   # v6:  mortalidad Plant_limited
        'Flujo_Mort_Micelio':  fv[:, 10],  # v11: mortalidad/reciclaje Mycelium
        'Flujo_Mort_Hongo':    fv[:, 13],  # v14: mortalidad Fungus

        # Transferencia simbiótica: v12, v15, v16
        'Flujo_Trans_C':  fv[:, 11],  # v12: C planta → hongo (= Intercambio_C)
        'Flujo_Trans_N':  fv[:, 14],  # v15: N hongo  → planta (= Intercambio_N)
        'Flujo_Trans_P':  fv[:, 15],  # v16: P hongo  → planta (= Intercambio_P)

        # Otros: fotosíntesis, captación nutrientes suelo, secuestro P, entradas
        'Flujo_Fotosintesis':  fv[:, 0],   # v1:  fotosíntesis
        'Flujo_N_Uptake':      fv[:, 6],   # v7:  captación N del suelo por micelio
        'Flujo_P_Uptake':      fv[:, 7],   # v8:  captación P del suelo por micelio
        'Flujo_P_Secuestro':   fv[:, 16],  # v17: secuestro P → P_no_avail
        'Flujo_Entrada_N':     fv[:, 17],  # v18: entrada N al suelo
        'Flujo_Entrada_P':     fv[:, 18],  # v19: entrada P al suelo
        # ────────────────────────────────────────────────────────────────

        # Variables del trade balance
        'N_Suelo':          n_soil,
        'P_Disponible':     p_avail,
        'P_Total_Suelo':    p_total_suelo,
        'N_Limitado':       n_lim,
        'P_Limitado':       p_lim,
        'C_Limitado':       c_lim,
        'N_Critico':        np.full(len(ts), n_critico),
        'P_Critico':        np.full(len(ts), p_critico),
        'Demanda_C_Hongo':  demanda_c_hongo,
        'Oferta_C_Planta':  oferta_c_planta,

        # Factores mecanísticos
        'Factor_N_Fotos':          factor_n_fotos,
        'Factor_C_Disponible':     factor_c_disponible,
        'Factor_N_Hongo_Limitado': factor_n_hongo_limitado,
        'Factor_N_Hongo_Libre':    factor_n_hongo_libre,

        # Intensidades
        'Escenario_I_Mutualismo_Clim':    intensidad_I,
        'Escenario_II_Mutualismo_Fuerte': intensidad_II,
        'Escenario_III_Comensalismo':     intensidad_III,
        'Escenario_IV_Parasitismo':       intensidad_IV,

        # Clasificación
        'Escenario_Dominante': escenario_dominante,

        # ROI
        'ROI_Planta': roi_planta,
        'ROI_Hongo':  roi_hongo,
        'ROI_N_C':    roi_N_C,
        'ROI_P_C':    roi_P_C,
    })    

    return indicadores, n_critico, p_critico


# ========================================
# 5. VISUALIZAR INDICADORES
# ========================================
def visualizar_indicadores(indicadores, n_critico, p_critico,
                            guardar=False,
                            nombre_archivo='indicadores_amf.png'):
    fig, axes = plt.subplots(3, 4, figsize=(18, 15))
    tiempo = indicadores['Tiempo']

    # ── fila 0 ────────────────────────────────────────────────────────
    axes[0, 0].plot(tiempo, indicadores['Biomasa_Planta'], 'g-', label='Planta')
    axes[0, 0].plot(tiempo, indicadores['Biomasa_Hongo'],  'b-', label='Hongo')
    axes[0, 0].set_title('Biomasas')
    axes[0, 0].set_xlabel('Tiempo')
    axes[0, 0].set_ylabel('Biomasa')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].plot(tiempo, indicadores['Oferta_C_Planta'], 'g-', label='Oferta planta')
    axes[0, 1].plot(tiempo, indicadores['Demanda_C_Hongo'], 'b-', label='Demanda hongo')
    axes[0, 1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axes[0, 1].set_title('Balance de Carbono')
    axes[0, 1].set_xlabel('Tiempo')
    axes[0, 1].set_ylabel('Tasa')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)

    ax3 = axes[0, 2]
    ax3.plot(tiempo, indicadores['N_Suelo'],      'g-', label='N_soil')
    ax3.plot(tiempo, indicadores['P_Disponible'], 'b-', label='P_avail')
    ax3.axhline(n_critico, color='g', linestyle='--', alpha=0.9,
                label=f'n_crit={n_critico:.3f}')
    ax3.axhline(p_critico, color='b', linestyle='--', alpha=0.9,
                label=f'p_crit={p_critico:.3f}')
    ax3.set_title('N y P en suelo')
    ax3.set_xlabel('Tiempo')
    ax3.set_ylabel('Concentración')
    ax3.legend(fontsize=7)
    ax3.grid(True, alpha=0.3)

    axes[0, 3].plot(tiempo, indicadores['Intercambio_C'], 'r-', label='C (Planta→Hongo)')
    axes[0, 3].plot(tiempo, indicadores['Intercambio_N'], 'g-', label='N (Hongo→Planta)')
    axes[0, 3].plot(tiempo, indicadores['Intercambio_P'], 'b-', label='P (Hongo→Planta)')
    axes[0, 3].set_title('Flujos de Intercambio')
    axes[0, 3].set_xlabel('Tiempo')
    axes[0, 3].set_ylabel('Tasa')
    axes[0, 3].legend()
    axes[0, 3].grid(True, alpha=0.3)


    # ── fila 1: ────────────────────────────────────────────────────────
    # [1,0] Crecimiento — v4 (planta), v10 (micelio), v13 (hongo)
    axes[1, 0].plot(tiempo, indicadores['Flujo_Crec_Planta'],
                    color='#2ecc71', lw=1.5, label='Planta (v4)')
    axes[1, 0].plot(tiempo, indicadores['Flujo_Crec_Micelio'],
                    color='#3498db', lw=1.5, label='Micelio (v10)')
    axes[1, 0].plot(tiempo, indicadores['Flujo_Crec_Hongo'],
                    color='#9b59b6', lw=1.5, label='Hongo (v13)')
    axes[1, 0].set_title('Flujos de Crecimiento')
    axes[1, 0].set_xlabel('Tiempo')
    axes[1, 0].set_ylabel('Tasa (v)')
    axes[1, 0].legend(fontsize=7)
    axes[1, 0].grid(True, alpha=0.3)

    # [1,1] Reciclaje / mortalidad — v5+v6 (planta), v11 (micelio), v14 (hongo)
    # Estos flujos liberan biomasa de vuelta al suelo como N y P
    mort_planta  = indicadores['Flujo_Mort_PlantA'] + indicadores['Flujo_Mort_PlantL']
    mort_micelio = indicadores['Flujo_Mort_Micelio']
    mort_hongo   = indicadores['Flujo_Mort_Hongo']
    axes[1, 1].plot(tiempo, mort_planta,
                    color='#27ae60', lw=1.5, label='Planta (v5+v6)')
    axes[1, 1].plot(tiempo, mort_micelio,
                    color='#2980b9', lw=1.5, label='Micelio (v11)')
    axes[1, 1].plot(tiempo, mort_hongo,
                    color='#8e44ad', lw=1.5, label='Hongo (v14)')
    axes[1, 1].fill_between(tiempo,
                             mort_planta + mort_micelio + mort_hongo,
                             alpha=0.12, color='gray', label='Total')
    axes[1, 1].set_title('Flujos de Reciclaje / Mortalidad')
    axes[1, 1].set_xlabel('Tiempo')
    axes[1, 1].set_ylabel('Tasa (v)')
    axes[1, 1].legend(fontsize=7)
    axes[1, 1].grid(True, alpha=0.3)

    # [1,2] Otros flujos: , entradas externas, secuestro P
    axes[1, 2].plot(tiempo, indicadores['Flujo_Fotosintesis'],
                    color='#f1c40f', lw=1.8, label='Fotosíntesis (v1)')
    axes[1, 2].plot(tiempo, indicadores['Flujo_N_Uptake'],
                    color='#16a085', lw=1.5, label='Captación N suelo (v7)')
    axes[1, 2].plot(tiempo, indicadores['Flujo_P_Uptake'],
                    color='#2471a3', lw=1.5, label='Captación P suelo (v8)')
    axes[1, 2].set_title('Flujos de fotosíntesis y captación N/P suelo')
    axes[1, 2].set_xlabel('Tiempo')
    axes[1, 2].set_ylabel('Tasa (v)')
    axes[1, 2].legend(fontsize=7)
    axes[1, 2].grid(True, alpha=0.3)

    # [1,3] Otros flujos: fotosíntesis, captación nutrientes suelo, secuestro P, entradas
    axes[1, 3].plot(tiempo, indicadores['Flujo_P_Secuestro'],
                    color='#a04000', lw=1.2, linestyle='--', label='Secuestro P (v17)')    
    axes[1, 3].plot(tiempo, indicadores['Flujo_Entrada_N'],
                    color='#117a65', lw=1.2, linestyle=':', label='Entrada N (v18)')
    axes[1, 3].plot(tiempo, indicadores['Flujo_Entrada_P'],
                    color='#1a5276', lw=1.2, linestyle=':', label='Entrada P (v19)')
    axes[1, 3].set_title('Otros Flujos')
    axes[1, 3].set_xlabel('Tiempo')
    axes[1, 3].set_ylabel('Tasa (v)')
    axes[1, 3].legend(fontsize=7)
    axes[1, 3].grid(True, alpha=0.3) 

    # ── fila 2 ────────────────────────────────────────────────────────
    axes[2, 0].plot(tiempo, indicadores['ROI_Planta'], 'g-', label='Planta (N+P)/C')
    axes[2, 0].plot(tiempo, indicadores['ROI_Hongo'],  'b-', label='Hongo C/(N+P)')
    axes[2, 0].axhline(y=1, color='k', linestyle='--', alpha=0.3)
    axes[2, 0].set_title('Retorno de Inversión (ROI)')
    axes[2, 0].set_xlabel('Tiempo')
    axes[2, 0].set_ylabel('ROI')
    axes[2, 0].legend()
    axes[2, 0].grid(True, alpha=0.3)
    axes[2, 0].set_ylim(0, 10)

    axes[2, 1].plot(tiempo, indicadores['ROI_N_C'], 'g-', label='N/C')
    axes[2, 1].plot(tiempo, indicadores['ROI_P_C'], 'b-', label='P/C')
    axes[2, 1].axhline(y=1, color='k', linestyle='--', alpha=0.3)
    axes[2, 1].set_title('Eficiencia N/C y P/C')
    axes[2, 1].set_xlabel('Tiempo')
    axes[2, 1].set_ylabel('Ratio')
    axes[2, 1].legend()
    axes[2, 1].grid(True, alpha=0.3)
    axes[2, 1].set_ylim(0, 5)

    c_transfer = indicadores['Intercambio_C'].iloc[-1]
    n_transfer = indicadores['Intercambio_N'].iloc[-1]
    p_transfer = indicadores['Intercambio_P'].iloc[-1]

    categorias = ['C\n(Planta→Hongo)', 'N\n(Hongo→Planta)', 'P\n(Hongo→Planta)']
    valores    = [c_transfer, n_transfer, p_transfer]
    colores    = ['#e74c3c', '#2ecc71', '#3498db']

    bars = axes[2, 2].bar(categorias, valores, color=colores,
                           edgecolor='black', linewidth=0.8, width=0.5)
    for bar, val in zip(bars, valores):
        axes[2, 2].text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.002,
                        f'{val:.4f}', ha='center', va='bottom', fontsize=9)
    axes[2, 2].set_title('Flujos C, N y P — valor final')
    axes[2, 2].set_ylabel('Tasa de transferencia final')
    axes[2, 2].axhline(y=0, color='black', linewidth=0.8, linestyle='--')
    axes[2, 2].grid(True, alpha=0.3, axis='y')

    # Figura de Johnson — scatter N vs P con radio = C (valor final)      
    n_final = indicadores['Intercambio_N'].iloc[-1]
    p_final = indicadores['Intercambio_P'].iloc[-1]
    c_final = indicadores['Intercambio_C'].iloc[-1]

    n_max = n_critico*2 #indicadores['Intercambio_N'].max()
    p_max = p_critico*2 #indicadores['Intercambio_P'].max()
    radio = c_final * 1000 # / (c_final + 1e-6)  # Radio proporcional a C, escalado para ser visible

    axes[2, 3].scatter(n_final, p_final,
                       s=radio,
                       color='purple', alpha=0.6,
                       edgecolors='black', linewidth=1.2)

    # Líneas de umbrales críticos
    axes[2, 3].axvline(x=n_critico, color='green', linewidth=1.0,
                       linestyle='--', alpha=0.9)
    axes[2, 3].axhline(y=p_critico, color='blue', linewidth=1.0,
                       linestyle='--', alpha=0.9)

    # Etiqueta con los valores
    axes[2, 3].annotate(
        f'N={n_final:.4f}\nP={p_final:.4f}\nC={c_final:.4f}',
        xy=(n_final, p_final),
        xytext=(10, 10), textcoords='offset points',
        fontsize=8, color='black',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7)
    )

    axes[2, 3].set_title('Espacio N–P de Johnson') #\n(radio = tasa final de C)')
    axes[2, 3].set_xlabel('Tasa final N (Hongo→Planta)')
    axes[2, 3].set_ylabel('Tasa final P (Hongo→Planta)')
    axes[2, 3].legend(fontsize=7)
    axes[2, 3].set_xlim(0, max(n_max,n_final+0.1))    
    axes[2, 3].set_ylim(0, max(p_max,p_final+0.1))    
    axes[2, 3].grid(True, alpha=0.3) 

    plt.tight_layout()

    if guardar:
        os.makedirs(os.path.dirname(nombre_archivo), exist_ok=True)
        plt.savefig(nombre_archivo, dpi=300, bbox_inches='tight')
        print(f"Gráfico guardado como '{nombre_archivo}'")    

    plt.show()

def visualizar_comparativa_escenarios(diccionario_escenarios, n_critico, p_critico,
                                     n_puntos_tiempo=None,
                                     guardar=False,
                                     nombre_archivo=os.path.join(VIS_DIR, '4_comparacion_escenarios',
                                                                 'comparativa_final_johnson.png')):
    """
    Recibe un diccionario de DataFrames para comparar múltiples escenarios.
    - Si n_puntos_tiempo == 1 → muestra SOLO el punto final (sin línea de trayectoria)
    - Si n_puntos_tiempo > 1 → muestra trayectoria con puntos coloreados por tiempo
    """
    # Marcadores por cuadrante del espacio N-P de Johnson
    # El cuadrante se determina según la posición relativa a n_critico y p_critico
    def get_marker(nv, pv):
        """
        Cuadrante I   (N<n_crit, P<p_crit) → triángulo  ▲
        Cuadrante II  (N>n_crit, P<p_crit) → cuadrado   ■
        Cuadrante III (N<n_crit, P>p_crit) → asterisco  *
        Cuadrante IV  (N>n_crit, P>p_crit) → círculo    ●
        """
        if nv <= n_critico and pv <= p_critico:
            return '^'   # I:   triángulo
        elif nv >  n_critico and pv <= p_critico:
            return 's'   # II:  cuadrado
        elif nv <= n_critico and pv >  p_critico:
            return '*'   # III: asterisco
        else:
            return 'o'   # IV:  círculo

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    k_radio = 1000

    if n_puntos_tiempo is None:
        n_puntos_tiempo = 10
    n_puntos_tiempo = max(1, int(n_puntos_tiempo))

    primer_df = next(iter(diccionario_escenarios.values()))

    n_total = len(primer_df)
    if n_puntos_tiempo == 1:
        idx_pts = [n_total - 1]
    else:
        idx_pts = np.linspace(0, n_total - 1, n_puntos_tiempo, dtype=int)

    if n_puntos_tiempo == 1:
        colors_t    = ['#d62728'] * 1
        tiempos_leg = [primer_df['Tiempo'].iloc[-1]]
    else:
        colors_t    = plt.cm.viridis(np.linspace(0, 1, n_puntos_tiempo))
        tiempos_leg = primer_df['Tiempo'].iloc[idx_pts].values

    for i, (nombre, df) in enumerate(diccionario_escenarios.items()):
        if i >= 4:
            break

        ax = axes[i]
        ns, ps = [], []

        for j, idx in enumerate(idx_pts):
            nv = df['Intercambio_N'].iloc[idx]
            pv = df['Intercambio_P'].iloc[idx]
            cv = df['Intercambio_C'].iloc[idx]
            marker = get_marker(nv, pv)

            # Asterisco necesita tamaño mayor para ser visible
            size = max(cv * k_radio, 10) * (2.0 if marker == '*' else 1.0)

            ax.scatter(nv, pv,
                       s=size,
                       marker=marker,
                       color=colors_t[j],
                       alpha=0.8,
                       edgecolors='black',
                       linewidth=0.8,
                       zorder=3)
            ns.append(nv)
            ps.append(pv)

            # ax.annotate(
            #     f'N={nv:.4f}\nP={pv:.4f}\nC={cv:.4f}',
            #     xy=(nv, pv),
            #     xytext=(10, 10), textcoords='offset points',
            #     fontsize=8, color='black',
            #     bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7)
            # )
            scatter_objects = []  # acumular los scatter de este eje

            for j, idx in enumerate(idx_pts):
                nv = df['Intercambio_N'].iloc[idx]
                pv = df['Intercambio_P'].iloc[idx]
                cv = df['Intercambio_C'].iloc[idx]
                marker = get_marker(nv, pv)

                size = max(cv * k_radio, 10) * (2.0 if marker == '*' else 1.0)

                sc = ax.scatter(nv, pv,
                                s=size,
                                marker=marker,
                                color=colors_t[j],
                                alpha=0.8,
                                edgecolors='black',
                                linewidth=0.8,
                                zorder=3,
                                label=f'N={nv:.4f}\nP={pv:.4f}\nC={cv:.4f}')
                ns.append(nv)
                ps.append(pv)
                scatter_objects.append((sc, nv, pv, cv))

            # Tooltip interactivo: aparece solo al pasar el mouse
            for sc, nv, pv, cv in scatter_objects:
                cursor = mplcursors.cursor(sc, hover=True)

                @cursor.connect("add")
                def on_add(sel, _nv=nv, _pv=pv, _cv=cv):
                    sel.annotation.set_text(
                        f'N={_nv:.4f}\nP={_pv:.4f}\nC={_cv:.4f}'
                    )
                    sel.annotation.get_bbox_patch().set(
                        facecolor='lightyellow', alpha=0.9,
                        boxstyle='round,pad=0.4'
                    )
                    sel.annotation.set_fontsize(8)

                @cursor.connect("remove")
                def on_remove(sel):
                    sel.annotation.set_visible(False)

        if n_puntos_tiempo > 1:
            ax.plot(ns, ps, color='gray', lw=1, linestyle='--', alpha=0.5, zorder=2)

        ax.axvline(x=n_critico, color='green', lw=1.2, linestyle=':', alpha=0.7, label='N Crítico')
        ax.axhline(y=p_critico, color='blue',  lw=1.2, linestyle=':', alpha=0.7, label='P Crítico')

        # Etiquetas de cuadrante en el fondo
        x_lo = 0
        x_hi = max(max(ns) * 1.3, n_critico * 1.5, 0.01) if ns else 1.0
        y_lo = 0
        y_hi = max(max(ps) * 1.3, p_critico * 1.5, 0.01) if ps else 1.0

        offset = 0.02
        ax.text(n_critico * 0.5,        p_critico * 0.5,        'I ▲',  ha='center', va='center', fontsize=9, color='gray', alpha=0.5)
        ax.text(n_critico * 1.5,        p_critico * 0.5,        'II ■', ha='center', va='center', fontsize=9, color='gray', alpha=0.5)
        ax.text(n_critico * 0.5,        p_critico * 1.5,        'III ✳', ha='center', va='center', fontsize=9, color='gray', alpha=0.5)
        ax.text(n_critico * 1.5,        p_critico * 1.5,        'IV ●', ha='center', va='center', fontsize=9, color='gray', alpha=0.5)

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_lo, y_hi)
        ax.set_title(f"Escenario: {nombre}", fontsize=12, fontweight='bold')
        ax.set_xlabel('N (Hongo→Planta)')
        ax.set_ylabel('P (Hongo→Planta)')
        ax.grid(True, alpha=0.25)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # Leyenda temporal (color)
    if n_puntos_tiempo == 1:
        title_legend = "Estado Final (tamaño = Flujo C)"
    else:
        title_legend = "Progreso Temporal (Color) y Flujo C (Tamaño)"

    legend_time = [
        Line2D([0], [0], marker='o', color='w',
               label=f't={int(round(t))}' if n_puntos_tiempo > 1 else f't_final={int(round(t))}',
               markerfacecolor=colors_t[k], markersize=9, markeredgecolor='k')
        for k, t in enumerate(tiempos_leg)
    ]

    # Leyenda de marcadores por cuadrante
    legend_markers = [
        Line2D([0], [0], marker='^', color='w', label='I: N↓ P↓',  markerfacecolor='gray', markersize=9, markeredgecolor='k'),
        Line2D([0], [0], marker='s', color='w', label='II: N↑ P↓', markerfacecolor='gray', markersize=9, markeredgecolor='k'),
        Line2D([0], [0], marker='*', color='w', label='III: N↓ P↑', markerfacecolor='gray', markersize=11, markeredgecolor='k'),
        Line2D([0], [0], marker='o', color='w', label='IV: N↑ P↑', markerfacecolor='gray', markersize=9, markeredgecolor='k'),
    ]

    leg1 = fig.legend(handles=legend_time,
                      loc='lower center',
                      ncol=min(5, len(legend_time)),
                      bbox_to_anchor=(0.5, 0.01),
                      title=title_legend,
                      fontsize=8)

    fig.legend(handles=legend_markers,
               loc='lower center',
               ncol=4,
               bbox_to_anchor=(0.5, 0.06),
               title='Cuadrante (forma)',
               fontsize=8)

    fig.add_artist(leg1)  # preservar ambas leyendas

    plt.tight_layout(rect=[0, 0.13, 1, 0.94])

    if guardar:
        if not nombre_archivo.endswith('.png'):
            nombre_archivo += '.png'
        archivo_completo = os.path.join(VIS_DIR, nombre_archivo)
        os.makedirs(os.path.dirname(archivo_completo), exist_ok=True)
        plt.savefig(archivo_completo, dpi=300, bbox_inches='tight')
        print(f"Comparativa guardada en: '{archivo_completo}'")

    plt.show()

# ========================================
# 6. EXPORTAR RESULTADOS
# ========================================
def exportar_indicadores(indicadores, n_critico, p_critico,
                          nombre_archivo='indicadores_amf.csv'):
    dir_out        = os.path.dirname(nombre_archivo)
    nombre_base    = os.path.basename(nombre_archivo)
    nombre_sin_ext = nombre_base.replace('.csv', '')

    os.makedirs(dir_out, exist_ok=True)

    indicadores.to_csv(nombre_archivo, index=False)
    print(f"Indicadores exportados a '{nombre_archivo}'")

    resumen = indicadores.describe()
    archivo_resumen = os.path.join(dir_out, f'resumen_{nombre_base}')
    resumen.to_csv(archivo_resumen)
    print(f"Estadísticas resumidas en '{archivo_resumen}'")

    proporciones = indicadores['Escenario_Dominante'].value_counts(normalize=True) * 100
    print("\nProporción de tiempo en cada escenario:")
    for esc, pct in proporciones.items():
        print(f"  {esc}: {pct:.1f}%")

    # umbrales = pd.DataFrame({
    #     'Variable': ['n_critico', 'p_critico'],
    #     'Valor':    [n_critico, p_critico]
    # })
    # archivo_umbrales = os.path.join(dir_out, f'umbrales_{nombre_sin_ext}.csv')
    # umbrales.to_csv(archivo_umbrales, index=False)
    # print(f"Umbrales dinámicos guardados en '{archivo_umbrales}'")


# ========================================
# 7. GENERAR REPORTE .TXT
# ========================================
VIS_DIR2 = 'projects/AMF_RN'
REPORT_DIR = os.path.join(VIS_DIR2, 'reports')
os.makedirs(REPORT_DIR, exist_ok=True)

def generar_reporte_escenario(indicadores, n_critico, p_critico,
                               escenario, nombre_escenario, nombre_params,
                               t_span, n_steps):
    """Genera un reporte .txt con los resultados de un escenario."""
    archivo = os.path.join(REPORT_DIR, f'reporte_escenario{escenario}.txt')

    proporciones = indicadores['Escenario_Dominante'].value_counts(normalize=True) * 100
    dominante    = proporciones.index[0] if len(proporciones) > 0 else 'N/A'

    with open(archivo, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write(f"REPORTE ESCENARIO {escenario}: {nombre_escenario}\n")
        f.write("=" * 60 + "\n\n")

        # ── Configuración ──────────────────────────────────────────────
        f.write("CONFIGURACIÓN\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Versión de parámetros : {nombre_params}\n")
        f.write(f"  Tiempo de simulación  : {t_span[0]} – {t_span[1]} días\n")
        f.write(f"  Pasos de integración  : {n_steps}\n")
        f.write(f"  Umbral N crítico      : {n_critico:.4f}\n")
        f.write(f"  Umbral P crítico      : {p_critico:.4f}\n\n")

        # ── Clasificación dominante ────────────────────────────────────
        f.write("CLASIFICACIÓN FUNCIONAL (Trade Balance Model)\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Escenario dominante   : {dominante}\n\n")
        f.write("  Proporción de tiempo en cada escenario:\n")
        for esc, pct in proporciones.items():
            marca = " ◄ dominante" if esc == dominante else ""
            f.write(f"    {esc:<35} {pct:6.1f}%{marca}\n")
        f.write("\n")

        # ── Biomasas finales ───────────────────────────────────────────
        f.write("BIOMASAS (valor final)\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Biomasa Planta : {indicadores['Biomasa_Planta'].iloc[-1]:.6f}\n")
        f.write(f"  Biomasa Hongo  : {indicadores['Biomasa_Hongo'].iloc[-1]:.6f}\n")
        f.write(f"  Plant_Active   : {indicadores['Plant_Active'].iloc[-1]:.6f}\n")
        f.write(f"  Plant_Limited  : {indicadores['Plant_Limited'].iloc[-1]:.6f}\n")
        f.write(f"  Fungus         : {indicadores['Fungus'].iloc[-1]:.6f}\n")
        f.write(f"  Mycelium       : {indicadores['Mycelium'].iloc[-1]:.6f}\n\n")

        # ── Intercambios simbióticos ───────────────────────────────────
        f.write("INTERCAMBIOS SIMBIÓTICOS (valor final)\n")
        f.write("-" * 40 + "\n")
        f.write(f"  C (Planta → Hongo) : {indicadores['Intercambio_C'].iloc[-1]:.6f}\n")
        f.write(f"  N (Hongo → Planta) : {indicadores['Intercambio_N'].iloc[-1]:.6f}\n")
        f.write(f"  P (Hongo → Planta) : {indicadores['Intercambio_P'].iloc[-1]:.6f}\n\n")

        # ── ROI ────────────────────────────────────────────────────────
        f.write("RETORNO DE INVERSIÓN (media temporal)\n")
        f.write("-" * 40 + "\n")
        f.write(f"  ROI Planta (N+P)/C : {indicadores['ROI_Planta'].mean():.6f}\n")
        f.write(f"  ROI Hongo  C/(N+P) : {indicadores['ROI_Hongo'].mean():.6f}\n")
        f.write(f"  Eficiencia N/C     : {indicadores['ROI_N_C'].mean():.6f}\n")
        f.write(f"  Eficiencia P/C     : {indicadores['ROI_P_C'].mean():.6f}\n\n")

        # ── Nutrientes en suelo ────────────────────────────────────────
        f.write("NUTRIENTES EN SUELO (estadísticas)\n")
        f.write("-" * 40 + "\n")
        for var, label in [('N_Suelo', 'N_soil'), ('P_Disponible', 'P_avail')]:
            s = indicadores[var]
            f.write(f"  {label}:\n")
            f.write(f"    min={s.min():.4f}  max={s.max():.4f}  "
                    f"media={s.mean():.4f}  std={s.std():.4f}\n")
        f.write("\n")

        # ── Intensidades medias ────────────────────────────────────────
        f.write("INTENSIDADES MEDIAS POR CUADRANTE\n")
        f.write("-" * 40 + "\n")
        cols = {
            'I: Mutualismo C-limitado':  'Escenario_I_Mutualismo_Clim',
            'II: Mutualismo fuerte':     'Escenario_II_Mutualismo_Fuerte',
            'III: Comensalismo':         'Escenario_III_Comensalismo',
            'IV: Parasitismo':           'Escenario_IV_Parasitismo',
        }
        for label, col in cols.items():
            f.write(f"  {label:<35} {indicadores[col].mean():.6f}\n")
        f.write("\n")

        f.write("=" * 60 + "\n")
        f.write("FIN DEL REPORTE\n")
        f.write("=" * 60 + "\n")

    print(f"Reporte guardado en '{archivo}'")
    return archivo


def generar_reporte_comparativo(resultados):
    """Genera un reporte .txt comparando todos los escenarios."""
    archivo = os.path.join(REPORT_DIR, 'reporte_comparativo.txt')

    nombres = {
        1: "Mutualismo limitado por C (N↓ P↓ C↓)",
        2: "Mutualismo fuerte (N↑ P↓ C↑)",
        3: "Comensalismo (N↓ P↑ C↓)",
        4: "Parasitismo (N↑ P↑ C↑)",
    }

    with open(archivo, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("REPORTE COMPARATIVO — TODOS LOS ESCENARIOS AMF\n")
        f.write("Trade Balance Model\n")
        f.write("=" * 60 + "\n\n")

        for esc in [1, 2, 3, 4]:
            ind   = resultados[esc]['indicadores']
            n_crit = resultados[esc]['n_critico']
            p_crit = resultados[esc]['p_critico']
            props  = ind['Escenario_Dominante'].value_counts(normalize=True) * 100
            dominante = props.index[0] if len(props) > 0 else 'N/A'

            f.write(f"── ESCENARIO {esc}: {nombres.get(esc, '')} ──\n")
            f.write(f"  Escenario dominante  : {dominante}\n")
            f.write(f"  Biomasa Planta final : {ind['Biomasa_Planta'].iloc[-1]:.4f}\n")
            f.write(f"  Biomasa Hongo  final : {ind['Biomasa_Hongo'].iloc[-1]:.4f}\n")
            f.write(f"  ROI Planta (media)   : {ind['ROI_Planta'].mean():.4f}\n")
            f.write(f"  ROI Hongo  (media)   : {ind['ROI_Hongo'].mean():.4f}\n")
            f.write(f"  Intercambio C (media): {ind['Intercambio_C'].mean():.4f}\n")
            f.write(f"  Intercambio N (media): {ind['Intercambio_N'].mean():.4f}\n")
            f.write(f"  Intercambio P (media): {ind['Intercambio_P'].mean():.4f}\n")
            f.write(f"  n_critico            : {n_crit:.4f}\n")
            f.write(f"  p_critico            : {p_crit:.4f}\n")
            f.write("  Proporciones:\n")
            for label, pct in props.items():
                f.write(f"    {label:<35} {pct:6.1f}%\n")
            f.write("\n")

        f.write("=" * 60 + "\n")
        f.write("FIN DEL REPORTE COMPARATIVO\n")
        f.write("=" * 60 + "\n")

    print(f"Reporte comparativo guardado en '{archivo}'")
    return archivo

# ========================================
# 8. EJECUTAR UN ESCENARIO
# ========================================

def ejecutar_escenario(escenario=1, guardar=True):
    print(f"\n{'='*55}")
    print(f"EJECUTANDO ESCENARIO {escenario}")
    print(f"{'='*55}")

    params, x0, nombre_params, nombre_escenario = get_scenario_config(escenario)
    rate_list = create_rate_list()

    print(f"Escenario: {nombre_escenario}")
    print(f"Versión de parámetros: {nombre_params}")

    t_span  = (0, 360*10)
    n_steps = 500*3  # [FIX-3] paso temporal más fino

    print("\nEjecutando simulación...")
    time_series, flux_vector = simulation(
        rn, x0=x0, rate=rate_list,
        spec_vector=params,
        t_span=t_span,
        n_steps=n_steps,
        method='LSODA',   
        rtol=1e-8,        
        atol=1e-10,        
    )
    
    # ── Rutas usando VIS_DIR ───────────────────────────────────────────
    ts_path   = os.path.join(VIS_DIR, '1_time_series',  f'{escenario}_ts.png')
    flux_path = os.path.join(VIS_DIR, '2_flux_vector',  f'{escenario}_fd.png')
    ind_path  = os.path.join(VIS_DIR, '3_indicadores',  f'escenario{escenario}')

    plot_amf_dynamics_new(
        time_series,
        title=f"Escenario: {nombre_escenario} - Temporal Dynamics",
        save_path=ts_path
    )
    plot_flux_dynamics_new(
        flux_vector,
        title=f"Escenario: {nombre_escenario} - Flux Dynamics",
        save_path=flux_path
    )
    print("Simulación completada!")

    print("\nCalculando indicadores...")
    indicadores, n_critico, p_critico = calcular_indicadores_amf(
        time_series, flux_vector, species_amf,
        n_reactions=len(rate_list),
        t_span=t_span,
        n_steps=n_steps,
    )

    proporciones = indicadores['Escenario_Dominante'].value_counts(normalize=True) * 100
    print("\nProporción de tiempo en cada escenario:")
    for esc, pct in proporciones.items():
        print(f"  {esc}: {pct:.1f}%")

    esperados = {1: 'I', 2: 'II', 3: 'III', 4: 'IV'}
    dominante_label = proporciones.index[0] if len(proporciones) > 0 else 'N/A'
    dominante_num   = dominante_label.split(':')[0].strip()
    if dominante_num == esperados.get(escenario, '?'):
        print(f"  ✓ Escenario dominante coincide con el esperado ({dominante_label})")
    else:
        print(f"  ⚠ Escenario dominante '{dominante_label}' no coincide "
              f"con el esperado (Escenario {esperados.get(escenario)})")

    if guardar:
        print("\nGenerando visualizaciones...")
        os.makedirs(os.path.join(VIS_DIR, '3_indicadores'), exist_ok=True)

        visualizar_indicadores(
            indicadores, n_critico, p_critico,
            guardar=True,
            nombre_archivo=f'{ind_path}.png',
        )
        print("\nExportando resultados...")
        exportar_indicadores(
            indicadores, n_critico, p_critico,
            nombre_archivo=f'{ind_path}.csv',
        )    
        generar_reporte_escenario(
            indicadores, n_critico, p_critico,
            escenario=escenario,
            nombre_escenario=nombre_escenario,
            nombre_params=nombre_params,
            t_span=t_span,
            n_steps=n_steps,
            )        

    return indicadores, n_critico, p_critico


def ejecutar_todos_escenarios():
    resultados = {}

    for escenario in [1, 2, 3, 4]:
        print(f"\n{'#'*60}")
        print(f"# EJECUTANDO ESCENARIO {escenario}")
        print(f"{'#'*60}")
        ind, n_crit, p_crit = ejecutar_escenario(escenario, guardar=True)
        resultados[escenario] = {
            'indicadores': ind,
            'n_critico':   n_crit,
            'p_critico':   p_crit,
        }

    print("\n" + "=" * 60)
    print("RESUMEN COMPARATIVO DE ESCENARIOS")
    print("=" * 60)

    comparacion = []
    for esc in [1, 2, 3, 4]:
        ind  = resultados[esc]['indicadores']
        props = ind['Escenario_Dominante'].value_counts(normalize=True) * 100

        fila = {
            'Escenario':        esc,
            'Biomasa_Planta_f': ind['Biomasa_Planta'].iloc[-1],
            'Biomasa_Hongo_f':  ind['Biomasa_Hongo'].iloc[-1],
            'ROI_Planta':    ind['ROI_Planta'].mean(),
            'ROI_Hongo':     ind['ROI_Hongo'].mean(),
            'Intercambio_C': ind['Intercambio_C'].mean(),
            'Intercambio_P': ind['Intercambio_P'].mean(),
            'n_critico':     resultados[esc]['n_critico'],
            'p_critico':     resultados[esc]['p_critico'],
            '% I':   props.get('I: Mutualismo C-limitado', 0),
            '% II':  props.get('II: Mutualismo fuerte',    0),
            '% III': props.get('III: Comensalismo',        0),
            '% IV':  props.get('IV: Parasitismo',          0),
        }
        comparacion.append(fila)

    df_comparacion = pd.DataFrame(comparacion)
    print("\n", df_comparacion.to_string())

    dir_indicadores = os.path.join(VIS_DIR, '4_comparacion_escenarios')
    os.makedirs(dir_indicadores, exist_ok=True)
    archivo_comp = os.path.join(dir_indicadores, 'comparacion_escenarios.csv')
    df_comparacion.to_csv(archivo_comp, index=False)
    print(f"\nTabla comparativa guardada en '{archivo_comp}'")

    return resultados


# ========================================
# 9. EJECUCIÓN PRINCIPAL
# ========================================

if __name__ == "__main__":
    print("=" * 55)
    print("ANÁLISIS DE INDICADORES FUNCIONALES AMF")
    print("TRADE BALANCE MODEL — Versión 6 (CORREGIDA)")
    print("=" * 55)
    print(f"\nLos archivos se guardarán en: {VIS_DIR}/")

    ESCENARIO_A_EJECUTAR = 4
    EJECUTAR_TODOS = True

    if EJECUTAR_TODOS:
        print("\nEjecutando TODOS los escenarios (1, 2, 3, 4)...")
        resultados = ejecutar_todos_escenarios()
        # FIX: Pass a dict of the 'indicadores' DataFrames only
        dic_indicadores = {k: v['indicadores'] for k, v in resultados.items()}
        visualizar_comparativa_escenarios(
            diccionario_escenarios=dic_indicadores,  # Use this instead of resultados
            n_critico=0.12,
            p_critico=0.1,
            n_puntos_tiempo = 5,
            guardar=True,
            nombre_archivo='4_comparacion_escenarios/comparativa_final_johnson.png'
        )
    else:
        print(f"\nEjecutando ESCENARIO {ESCENARIO_A_EJECUTAR}...")
        indicadores, n_critico, p_critico = ejecutar_escenario(
            ESCENARIO_A_EJECUTAR, 
            guardar=True
        )

    print("\n" + "=" * 55)
    print("ANÁLISIS COMPLETADO")
    print(f"Todos los archivos se guardaron en: {VIS_DIR}/")
    print("=" * 55)