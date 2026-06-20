# ==========================================
# MONTE CARLO ROBUSTNESS ANALYSIS FOR AMF
# WITH COLOR-CODED TRAJECTORIES BY INDICATOR - CORREGIDO
# ==========================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from pyCOT.simulations.ode import simulation

# ========================================
# BIOLOGICAL BASE PARAMETERS
# ========================================

x0_dict = {
    'Pact': 5.5,    # [1.00-7.50] g biomass/kg soil
    'Plim': 1.5,    # [0.50-2.50] g biomass/kg soil
    'Cp':   45.0,   # [35.0-50.0] % leaf C
    'Np':   2.0,    # [1.00-3.00] % leaf N
    'Pp':   0.4,    # [0.05-0.88] % leaf P
    'F':    20,     # [6.00-90.0] % root length colonized
    'M':    0.5,    # [0.10-3.00] m mycelium/g soil
    'Cf':   40.0,   # [35.0-50.0] % fungal C
    'Nf':   5.0,    # [1.00-6.00] % fungal N
    'Pf':   0.5,    # [0.10-2.00] % fungal P
    'Ns':   30,     # [10-40] mg N/kg dry soil
    'Ps':   25      # [5-20] mg P/kg dry soil
}


def get_base_params_dict():
    return {
        'r1': [3.5, 10.0], 
        'r2': [1e-5], 'r3': [0.015], 'r4': [100e-6],
        'r5': [5e-4], 'r6': [0.005], 'r7': [5e-8], 'r8': [5e-3],
        'r9': [1.5, 30.0], 'r10': [1.2, 10.0], 
        'r11': [5e-7], 'r12': [0.01], 'r13': [0.0007], 
        'r14': [0.8, 20.0],
        'r15': [0.5, 20.0], 
        'r16': [6e-7], 'r17': [0.001],
        'r18': [0.0], 'r19': [0.0]
    }

def get_rate_list_dict():
    return {
        'r1': 'mmk', 'r2': 'mak', 'r3': 'mak', 'r4': 'mak',
        'r5': 'mak', 'r6': 'mak', 'r7': 'mak', 'r8': 'mak',
        'r9': 'mmk', 'r10': 'mmk', 'r11': 'mak', 'r12': 'mak',
        'r13': 'mak', 'r14': 'mmk', 'r15': 'mmk', 'r16': 'mak',
        'r17': 'mak', 'r18': 'mak', 'r19': 'mak'
    }

def get_scenario_config(rn, scenario: str):
    species_order = [s.name for s in rn.species()]
    reaction_order = [r.name() if callable(r.name) else r.name for r in rn.reactions()]

    params_dict = get_base_params_dict()
    local_x0_dict = x0_dict.copy()

    if scenario == 'I':
        local_x0_dict['Ns'] = 10
        local_x0_dict['Ps'] = 5
        params_dict['r18'][0] = 0.001
        params_dict['r19'][0] = 0.01
        params_name   = "I – N↓ P↓"
        scenario_name = "I – C-limited mutualism"

    elif scenario == 'II':
        local_x0_dict['Ns'] = 45
        local_x0_dict['Ps'] = 5
        params_dict['r18'][0] = 0.5
        params_dict['r19'][0] = 0.1
        params_name   = "II – N↑ P↓"
        scenario_name = "II – Strong mutualism"

    elif scenario == 'III':
        local_x0_dict['Ns'] = 10
        local_x0_dict['Ps'] = 45
        params_dict['r18'][0] = 0.001
        params_dict['r19'][0] = 2.0
        params_name   = "III – N↓ P↑"
        scenario_name = "III – Commensalism"

    elif scenario == 'IV':
        local_x0_dict['Ns'] = 90 #100 # 90
        local_x0_dict['Ps'] = 90 #95 # 100
        params_dict['r18'][0] = 0.01 # 0.01
        params_dict['r19'][0] = 4.5 # 4
        params_name   = "IV – N↑ P↑"
        scenario_name = "IV – Parasitism"
    else:
        raise ValueError("Invalid scenario")

    x0 = [local_x0_dict[s] for s in species_order]
    params = [params_dict[r] for r in reaction_order]
    rate_list = [get_rate_list_dict()[r] for r in reaction_order]

    return params, x0, rate_list, params_name, scenario_name


# ==========================================
# PERTURBATION UTILITIES
# ==========================================

def perturb_value(val, variance=0.10):
    sigma = abs(val) * variance
    sampled = np.random.normal(loc=val, scale=sigma)
    return max(sampled, 0.0)

def perturb_vector(vec, variance=0.10):
    return [perturb_value(v, variance) for v in vec]

def perturb_nested_params(params, variance=0.10):
    perturbed = []
    for group in params:
        new_group = [perturb_value(v, variance) for v in group]
        perturbed.append(new_group)
    return perturbed


# ==========================================
# STEADY STATE INDICATORS
# ==========================================

def calculate_amf_indicators(flux_vector, t_span, steady_frac=0.2):
    """Calculate mycorrhizal function indicators at steady state."""
    t_ss_init = t_span[1] - steady_frac * (t_span[1] - t_span[0])  # tiempo de inicio del SS
    time_arr  = flux_vector['Time'].values                           # array de tiempos
    n_crit    = np.searchsorted(time_arr, t_ss_init)                 # índice >= t_ss_init

    def ss(rxn):
        return float(flux_vector[rxn].iloc[n_crit:].mean())         # promedio en zona SS

    v1  = ss('r1')
    v6  = ss('r6')
    v13 = ss('r13')
    v15 = ss('r15')

    C_demand  = v13 / (v1  + 1e-10)
    P_benefit = v15 / (v6  + v15 + 1e-10)

    ind = {'C_supply': v1, 'C_demand': C_demand, 'P_benefit': P_benefit}
    return ind, n_crit


# ==========================================
# ECOLOGICAL / JOHNSON CONDITIONS
# ==========================================

def evaluate_johnson_conditions(indicators, scenario):
    C_supply  = indicators['C_supply']
    C_demand  = indicators['C_demand']
    P_benefit = indicators['P_benefit']
    L_high = 0.5
    S_high = 0.25
    C_high = 0.75
    P_high = 0.25

    if scenario == 'I':
        cond_order     = (P_benefit > C_supply >= C_demand)
        cond_p_benefit = (P_benefit >= L_high)
        cond_c_supply  = (C_supply < L_high)
        all_ok = (cond_order and cond_p_benefit and cond_c_supply)
        return {'order_ok': cond_order, 'p_benefit_ok': cond_p_benefit,
                'c_supply_ok': cond_c_supply, 'all_ok': all_ok}

    elif scenario == 'II':
        cond_c_demand  = (C_supply > C_demand > S_high)
        cond_p_benefit = (P_benefit > S_high)
        cond_order     = (C_supply > P_benefit > S_high)
        all_ok = (cond_c_demand and cond_p_benefit and cond_order)
        return {'c_demand_ok': cond_c_demand, 'p_benefit_ok': cond_p_benefit,
                'order_ok': cond_order, 'all_ok': all_ok}

    elif scenario == 'III':
        cond_c_supply  = (C_demand < C_supply <= C_high)
        cond_p_benefit = (C_demand < P_benefit <= C_high)
        all_ok = (cond_c_supply and cond_p_benefit)
        return {'c_supply_ok': cond_c_supply, 'p_benefit_ok': cond_p_benefit,
                'all_ok': all_ok}

    elif scenario == 'IV':
        cond_order = (P_benefit < P_high < C_supply < C_demand)
        return {'order_ok': cond_order, 'all_ok': cond_order}

    else:
        raise ValueError("Invalid scenario")


# ==========================================
# HELPERS PARA EXPANSIÓN EN EXCEL
# ==========================================

def expand_x0_to_dict(x0_mc, species_order):
    """
    Convierte el vector x0 perturbado en un dict con nombres de especie.
    Ejemplo: {'x0_Pact': 5.32, 'x0_Plim': 1.48, ...}
    """
    return {f'x0_{sp}': val for sp, val in zip(species_order, x0_mc)}


def expand_params_to_dict(params_mc, reaction_order):
    """
    Convierte la lista de listas de parámetros perturbados en un dict plano.
    Si una reacción tiene un solo valor  → columna 'p_r1'
    Si tiene múltiples valores (Vmax, Km) → columnas 'p_r1_0', 'p_r1_1', ...
    Ejemplo: {'p_r1_0': 2.31, 'p_r1_1': 9.87, 'p_r2': 9.8e-6, ...}
    """
    result = {}
    for rname, vals in zip(reaction_order, params_mc):
        if len(vals) == 1:
            result[f'p_{rname}'] = vals[0]
        else:
            for i, v in enumerate(vals):
                result[f'p_{rname}_{i}'] = v
    return result


# ==========================================
# RUN SINGLE SCENARIO
# ==========================================

def run_single_scenario(rn, scenario, variance, n_simulations,
                        t_span, n_steps, steady_frac):

    print("=" * 70)
    print(" MONTE CARLO ROBUSTNESS ANALYSIS ")
    print("=" * 70)
    print(f"Scenario:      {scenario}")
    print(f"Variance:      {variance*100:.1f}%")
    print(f"Simulations:   {n_simulations}")

    # Obtenemos el orden real de especies y reacciones para poder etiquetar
    species_order  = [s.name for s in rn.species()]
    reaction_order = [r.name() if callable(r.name) else r.name
                      for r in rn.reactions()]

    params_base, x0_base, rate_list, _, scenario_name = \
        get_scenario_config(rn, scenario)

    stats = []
    condition_counts = {}
    success_all = 0

    for k in range(n_simulations):
        print(f"Simulation {k+1}/{n_simulations}", end='\r')

        x0_mc     = perturb_vector(x0_base, variance)
        params_mc = perturb_nested_params(params_base, variance)

        try:
            ts, flux = simulation(
                rn, x0=x0_mc, rate=rate_list, spec_vector=params_mc,
                t_span=t_span, n_steps=n_steps, method='LSODA',
                rtol=1e-8, atol=1e-10
            )

            indicators, _ = calculate_amf_indicators(flux, t_span, steady_frac)
            ev             = evaluate_johnson_conditions(indicators, scenario)

            success_all += int(ev['all_ok'])
            for key in ev:
                condition_counts[key] = condition_counts.get(key, 0) + int(ev[key])

            row = {
                'simulation': k + 1,
                **expand_x0_to_dict(x0_mc, species_order),
                **expand_params_to_dict(params_mc, reaction_order),
                'C_supply':  indicators['C_supply'],
                'C_demand':  indicators['C_demand'],
                'P_benefit': indicators['P_benefit'],
            }
            for key, val in ev.items():
                row[key] = val

            stats.append(row)

        except Exception as e:
            print(f"\nSimulation failed: {e}")

            # En caso de error también expandimos para mantener columnas consistentes
            row = {
                'simulation': k + 1,
                **{f'x0_{sp}':  np.nan for sp in species_order},
                **{f'p_{rn_}':  np.nan for rn_ in reaction_order},
                'C_supply':  np.nan,
                'C_demand':  np.nan,
                'P_benefit': np.nan,
                'all_ok':    False,
            }
            stats.append(row)

    df = pd.DataFrame(stats)

    print("\n")
    print("=" * 70)
    print(f" ROBUSTNESS RESULTS – {scenario_name}")
    print("=" * 70)

    for key, count in condition_counts.items():
        pct = 100 * count / n_simulations
        print(f"{key:<20}: {pct:.2f}%")

    print("-" * 70)
    pct_all = 100 * success_all / n_simulations
    print(f"{'ALL CONDITIONS':<20}: {pct_all:.2f}%")
    print("=" * 70)

    print("\nIndicator statistics:\n")
    print(df[['C_supply', 'C_demand', 'P_benefit']].describe())

    return df


# ==========================================
# VISUALIZATION CONFIGURATION
# ==========================================

INDICATOR_BASE_COLORS = {
    'c_supply':  '#2ECC71',
    'c_demand':  '#E67E22',
    'p_benefit': '#3498DB'
}

INDICATOR_NAMES = {
    'c_supply':  'C supply',
    'c_demand':  'C demand',
    'p_benefit': 'P benefit'
}

SCENARIO_COLORS = {
    'I':   '#2E86C1',
    'II':  '#28B463',
    'III': '#E67E22',
    'IV':  '#E74C3C'
}

SCENARIO_LABELS = {
    'I':   'I – C-limited mutualism (N↓P↓)',
    'II':  'II – Strong mutualism (N↑P↓)',
    'III': 'III – Commensalism (N↓P↑)',
    'IV':  'IV – Parasitism (N↑P↑)',
}


# ==========================================
# COLLECT TRAJECTORIES WITH SUCCESS/FAILURE
# ==========================================

def collect_monte_carlo_trajectories_with_labels(
    rn, scenario, variance=0.10, n_simulations=100,
    t_span=(0, 600), n_steps=601, steady_frac=0.2
):
    params_base, x0_base, rate_list, _, _ = get_scenario_config(rn, scenario)

    trajectories = {
        'time': None,
        'c_supply':  {'success': [], 'failure': []},
        'c_demand':  {'success': [], 'failure': []},
        'p_benefit': {'success': [], 'failure': []}
    }

    successful_sims = 0
    failed_sims     = 0
    all_ok_count    = 0

    for k in range(n_simulations):
        if (k + 1) % 20 == 0:
            print(f"    Simulation {k+1}/{n_simulations} "
                  f"(success_cond: {all_ok_count}, "
                  f"success_sim: {successful_sims}, "
                  f"failed_sim: {failed_sims})")

        x0_mc     = perturb_vector(x0_base, variance)
        params_mc = perturb_nested_params(params_base, variance)

        try:
            ts, flux = simulation(
                rn, x0=x0_mc, rate=rate_list, spec_vector=params_mc,
                t_span=t_span, n_steps=n_steps, method='LSODA',
                rtol=1e-8, atol=1e-10
            )

            t         = flux['Time'].values
            c_supply  = flux['r1'].values
            v13       = flux['r13'].values
            v15       = flux['r15'].values
            v6        = flux['r6'].values
            c_demand  = v13 / (c_supply + 1e-12)
            p_benefit = v15 / (v15 + v6 + 1e-12)

            indicators_ss, _ = calculate_amf_indicators(flux, t_span, steady_frac)
            ev               = evaluate_johnson_conditions(indicators_ss, scenario)
            is_success    = ev['all_ok']

            if trajectories['time'] is None:
                trajectories['time'] = t

            label = 'success' if is_success else 'failure'
            trajectories['c_supply'][label].append(c_supply)
            trajectories['c_demand'][label].append(c_demand)
            trajectories['p_benefit'][label].append(p_benefit)

            if is_success:
                all_ok_count += 1
            successful_sims += 1

        except Exception as e:
            failed_sims += 1
            continue

    print(f"    Final: {all_ok_count} ALL_CONDITIONS_OK, "
          f"{successful_sims - all_ok_count} FAILED_CONDITIONS, "
          f"{failed_sims} simulation failures")

    for indicator in ['c_supply', 'c_demand', 'p_benefit']:
        for label in ['success', 'failure']:
            arr = trajectories[indicator][label]
            trajectories[indicator][label] = np.array(arr) if arr else np.array([])

    return trajectories, successful_sims, failed_sims, all_ok_count


# ==========================================
# PLOT TRAJECTORIES
# ==========================================

def plot_monte_carlo_trajectories_by_indicator(
    rn, all_trajectories, show_success=True, show_failure=True,
    show_legend=True, show_median=False, save_path=None
):
    scenarios = ['I', 'II', 'III', 'IV']
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    layout = {'III': (0, 0), 'IV': (0, 1), 'I': (1, 0), 'II': (1, 1)}

    for scen in scenarios:
        row, col = layout[scen]
        ax = axes[row][col]

        if scen not in all_trajectories:
            ax.text(0.5, 0.5, f'No data for scenario {scen}',
                    ha='center', va='center', transform=ax.transAxes)
            continue

        trajectories    = all_trajectories[scen]['trajectories']
        successful_cond = all_trajectories[scen]['all_ok_count']
        total_sims      = all_trajectories[scen]['successful_sims']
        failed_cond     = total_sims - successful_cond

        if total_sims == 0:
            ax.text(0.5, 0.5, f'No successful simulations\nfor scenario {scen}',
                    ha='center', va='center', transform=ax.transAxes)
            continue

        time = trajectories['time']

        for indicator, color in INDICATOR_BASE_COLORS.items():
            if show_success and successful_cond > 0:
                n_traj = len(trajectories[indicator]['success'])
                if n_traj > 0:
                    for i in range(n_traj):
                        ax.plot(time, trajectories[indicator]['success'][i],
                                color=color, alpha=0.4, linewidth=0.6)

            if show_median and successful_cond > 0:
                n_traj = len(trajectories[indicator]['success'])
                if n_traj > 0:
                    median_s = np.median(trajectories[indicator]['success'], axis=0)
                    ax.plot(time, median_s, color=color, linewidth=3.0, alpha=1.0)

            if show_failure and failed_cond > 0:
                n_traj = len(trajectories[indicator]['failure'])
                if n_traj > 0:
                    for i in range(n_traj):
                        ax.plot(time, trajectories[indicator]['failure'][i],
                                color=color, alpha=0.2, linewidth=0.6, linestyle='--')

        if show_legend:
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color=c, linewidth=2.5, label=INDICATOR_NAMES[i])
                for i, c in INDICATOR_BASE_COLORS.items()
            ]
            if show_median and successful_cond > 0:
                legend_elements.append(
                    Line2D([0], [0], color='black', linewidth=2.5, label='Median (Success)'))
            if show_failure and failed_cond > 0:
                legend_elements.append(
                    Line2D([0], [0], color='#E74C3C', linewidth=2, linestyle='--',
                           label=f'Failure (n={failed_cond})', alpha=0.6))
            ax.legend(handles=legend_elements, fontsize=8, loc='best')

        success_rate = 100 * successful_cond / total_sims if total_sims > 0 else 0
        ax.set_title(
            f"{SCENARIO_LABELS[scen]}\n(Success rate: {success_rate:.1f}%)",
            fontweight='bold', fontsize=11, color=SCENARIO_COLORS[scen])
        ax.set_xlabel("Time (days)", fontsize=10)
        ax.set_ylabel("Indicator values", fontsize=10)

        all_values = []
        if show_failure:
            for ind in ['c_supply', 'c_demand', 'p_benefit']:
                if len(trajectories[ind]['failure']) > 0:
                    all_values.extend(trajectories[ind]['failure'].flatten())
        if show_median and successful_cond > 0:
            for ind in ['c_supply', 'c_demand', 'p_benefit']:
                if len(trajectories[ind]['success']) > 0:
                    all_values.extend(np.median(trajectories[ind]['success'], axis=0))

        if all_values:
            max_val = np.percentile(all_values, 99)
            ax.set_ylim(0, max_val * 1.1)

        ax.grid(True, alpha=0.25)
        stats_text = (f"Total: {total_sims}\n"
                      f"✓ Success: {successful_cond} ({success_rate:.1f}%)\n"
                      f"✗ Failure: {failed_cond}")
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=8,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    fig.text(0.53, 0.02,
             r"Limitation $\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad"
             r"\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad$ Luxury",
             ha='center', fontsize=14, color='#555')
    fig.text(0.53, 0.001, r"N supply", ha='center', fontsize=18,
             fontweight='bold', color='#555')
    fig.text(0.001, 0.53, r"P supply", va='center', rotation='vertical',
             fontsize=18, fontweight='bold', color='#555')
    fig.text(0.02, 0.52,
             r"Limitation $\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad"
             r"\qquad\qquad$ Luxury",
             va='center', rotation='vertical', fontsize=14, color='#555')

    fig.suptitle(
        'Monte Carlo Trajectories - Color Coding:\n'
        'Green = C_supply | Orange = C_demand | Blue = P_benefit\n'
        'Solid = Median of Success | Dashed = Failure trajectories',
        fontsize=12, y=0.98, fontweight='bold')

    plt.tight_layout(rect=[0.04, 0.04, 1, 0.96])

    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"\nFigure saved: {save_path}")

    plt.show()
    return fig


# ==========================================
# EXCEL SAVING HELPER
# ==========================================

def save_results_to_excel(all_results, scenarios_list, variance, n_simulations,
                           t_span, n_steps, steady_frac, random_seed, save_excel_path):
    """
    Guarda los resultados en Excel de forma estructurada:
      - Una hoja 'Summary'  con el porcentaje de éxito por escenario.
      - Una hoja por escenario con:
          · Bloque de configuración (parámetros de la corrida).
          · Tabla de datos con columnas legibles para x0 (x0_Pact, x0_Plim, ...)
            y para los parámetros cinéticos (p_r1_0, p_r1_1, p_r2, ...).
    """
    print(f"\nSaving results to Excel: {save_excel_path} ...")
    os.makedirs(os.path.dirname(save_excel_path), exist_ok=True)

    try:
        with pd.ExcelWriter(save_excel_path, engine='openpyxl') as writer:

            # --------------------------------------------------
            # Hoja resumen
            # --------------------------------------------------
            summary_rows = []
            for sc in scenarios_list:
                df_sc   = all_results[sc]
                n_total = len(df_sc)
                n_ok    = int(df_sc['all_ok'].sum()) if 'all_ok' in df_sc.columns else 0
                n_nan   = int(df_sc['C_supply'].isna().sum())
                summary_rows.append({
                    'Scenario':          sc,
                    'Total_simulations': n_total,
                    'Successful_runs':   n_total - n_nan,
                    'Failed_runs':       n_nan,
                    'All_conditions_ok': n_ok,
                    'Success_rate_%':    round(100 * n_ok / n_total, 2) if n_total > 0 else 0.0,
                })

            pd.DataFrame(summary_rows).to_excel(writer, sheet_name='Summary', index=False)

            # --------------------------------------------------
            # Una hoja por escenario
            # --------------------------------------------------
            for sc in scenarios_list:
                df_sc = all_results[sc]

                # Bloque de configuración en las primeras filas
                config_df = pd.DataFrame({
                    'Parameter': ['Scenario', 'Variance', 'Simulations',
                                  't_span', 'n_steps', 'steady_frac', 'random_seed'],
                    'Value':     [sc, variance, n_simulations,
                                  str(t_span), n_steps, steady_frac, random_seed]
                })
                config_df.to_excel(writer, sheet_name=f'Scenario_{sc}',
                                   index=False, startrow=0)

                # Tabla de datos debajo del bloque de configuración
                # (ya contiene columnas x0_* y p_* expandidas desde run_single_scenario)
                data_start_row = len(config_df) + 2  # fila vacía de separación
                df_sc.to_excel(writer, sheet_name=f'Scenario_{sc}',
                               index=False, startrow=data_start_row)

        print("  ✓ Excel saved successfully.")

    except Exception as e:
        print(f"  ✗ Error saving Excel: {e}")


# ==========================================
# MAIN MONTE CARLO FUNCTION
# ==========================================

def monte_carlo_amf_robustness(
    rn, scenario='all', variance=0.10, n_simulations=100,
    t_span=(0, 600), n_steps=601, steady_frac=0.2, random_seed=123,
    plot_trajectories=True, show_success=True, show_failure=True,
    show_median=False, save_figure_path=None, save_excel_path=None
):
    np.random.seed(random_seed)

    if scenario == 'all':
        scenarios_to_run = ['I', 'II', 'III', 'IV']
        print("=" * 70)
        print(" MONTE CARLO ROBUSTNESS ANALYSIS - ALL SCENARIOS ")
        print("=" * 70)

        all_results      = {}
        all_trajectories = {}

        for sc in scenarios_to_run:
            print(f"\n\n{'='*70}\n EXECUTING SCENARIO {sc} \n{'='*70}")
            df_scenario       = run_single_scenario(rn, sc, variance, n_simulations,
                                                    t_span, n_steps, steady_frac)
            all_results[sc]   = df_scenario

            if plot_trajectories:
                trajectories, n_success, n_failed, n_ok = \
                    collect_monte_carlo_trajectories_with_labels(
                        rn, sc, variance, n_simulations, t_span, n_steps, steady_frac)
                all_trajectories[sc] = {
                    'trajectories':   trajectories,
                    'successful_sims': n_success,
                    'failed_sims':     n_failed,
                    'all_ok_count':    n_ok
                }

        # -----------------------------------------------
        # CORRECCIÓN: se usa la función dedicada al Excel
        # -----------------------------------------------
        if save_excel_path:
            save_results_to_excel(
                all_results, scenarios_to_run,
                variance, n_simulations, t_span, n_steps, steady_frac, random_seed,
                save_excel_path
            )

        if plot_trajectories:
            fig = plot_monte_carlo_trajectories_by_indicator(
                rn, all_trajectories,
                show_success=show_success, show_failure=show_failure,
                show_legend=True, show_median=show_median,
                save_path=save_figure_path
            )

        return all_results, all_trajectories

    else:
        df_result = run_single_scenario(rn, scenario, variance, n_simulations,
                                        t_span, n_steps, steady_frac)

        if save_excel_path:
            save_results_to_excel(
                {scenario: df_result}, [scenario],
                variance, n_simulations, t_span, n_steps, steady_frac, random_seed,
                save_excel_path
            )

        if plot_trajectories:
            trajectories, n_success, n_failed, n_ok = \
                collect_monte_carlo_trajectories_with_labels(
                    rn, scenario, variance, n_simulations, t_span, n_steps, steady_frac)

            fig, ax = plt.subplots(1, 1, figsize=(12, 8))
            time        = trajectories['time']
            failed_cond = n_success - n_ok

            for indicator, color in INDICATOR_BASE_COLORS.items():
                if n_ok > 0 and len(trajectories[indicator]['success']) > 0:
                    alpha_s = max(0.08, min(0.4, 20.0 / n_ok))
                    for i in range(len(trajectories[indicator]['success'])):
                        ax.plot(time, trajectories[indicator]['success'][i],
                                color=color, alpha=alpha_s, linewidth=0.6)
                if failed_cond > 0 and len(trajectories[indicator]['failure']) > 0:
                    for i in range(len(trajectories[indicator]['failure'])):
                        ax.plot(time, trajectories[indicator]['failure'][i],
                                color=color, alpha=0.15, linewidth=0.4)

            success_rate = 100 * n_ok / n_success if n_success > 0 else 0
            ax.set_title(
                f'Monte Carlo Trajectories - Scenario {scenario}\n'
                f'Success rate: {success_rate:.1f}% ({n_ok}) | '
                f'Failure: {failed_cond}')
            ax.set_xlabel("Time (days)")
            ax.set_ylabel("Indicator values")

            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], color=c, lw=2.5, label=INDICATOR_NAMES[i])
                for i, c in INDICATOR_BASE_COLORS.items()
            ]
            ax.legend(handles=legend_elements, loc='best')

            if save_figure_path:
                plt.savefig(save_figure_path, dpi=300, bbox_inches='tight')

            plt.show()
            return df_result, trajectories

        return df_result


# ==========================================
# MAIN EXECUTION
# ==========================================

if __name__ == '__main__':

    from pyCOT.io.functions import read_txt
    from pyCOT.simulations.ode import simulation

    FILE_PATH = 'data/Ecological_models/AMF_2.txt'
    rn = read_txt(FILE_PATH)

    print("\n" + "="*80)
    print(" MONTE CARLO ROBUSTNESS ANALYSIS WITH INDICATOR-COLORED TRAJECTORIES ")
    print(" GREEN = C_supply | ORANGE = C_demand | BLUE = P_benefit")
    print(" SATURATED = Success (ALL CONDITIONS) | FADED = Failure")
    print("="*80 + "\n")

    results_all, trajectories_all = monte_carlo_amf_robustness(
        rn,
        scenario='all',
        variance=0.10,
        n_simulations=100,
        n_steps=601,
        plot_trajectories=True,
        show_success=True,
        show_failure=False,
        show_median=True,
        save_figure_path='projects/AMF_3/outputs_3/monte_carlo_trajectories_by_indicator_exitos.png',
        save_excel_path='projects/AMF_3/outputs_3/monte_carlo_results_all_scenarios.xlsx'
    )

    print("\n" + "="*80)
    print(" ANALYSIS COMPLETED SUCCESSFULLY ")
    print("="*80)