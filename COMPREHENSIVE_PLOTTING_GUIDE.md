# Comprehensive Plotting Guide
## All Scenario Combinations for Article Figures

**Status**: ✅ IMPLEMENTED - All 45 combinations now visualized
**Backend**: Non-interactive (Agg) - Plots saved directly, not displayed

---

## Overview

The script now generates **comprehensive visualizations** for ALL combinations:
- **3** conflict scenarios (low, medium, severe)
- **5** climate parameters (0.01, 0.02, 0.03, 0.04, 0.05)
- **3** budget renewal rates (0.33, 0.66, 1.0)
- **Total**: 3 × 5 × 3 = **45 combinations**

Each combination is tested with:
- 4 government strategies (development, security, balanced, adaptive)
- 3 armed group strategies (recruitment, displacement, balanced)
- **Total simulations**: 45 × 4 × 3 = **540 simulations**

---

## Output Directory Structure

```
visualizations/conflict_simulation_results_adaptive/
├── 01_budget_comparison.png
├── 02_net_budget_efficiency.png
├── 03_prosperity_baselines.png
├── 04_strong_resilience_baseline.png
├── 05_armed_groups_baseline.png
├── 06_adaptive_performance.png
├── 07_heatmap_effectiveness.png
├── 0A_budget_renewal_comparison_all.png
├── 0B_best_strategy_by_renewal.png
├── 0C_detailed_metrics_grid.png
├── 0D_overall_winner_analysis.png
├── 0E_win_rate_analysis.png
├── 0F_timeseries_all_budgets.png
├── 0G_timeseries_all_scenarios.png
├── 0H_budget_evolution_timeseries.png
├── 0I_process_mode_evolution.png
├── 0J_multiscale_mode_distribution.png
├── 0K_cognitive_domain_fraction.png
├── 0L_seasonal_variation.png
├── 0M_steady_state_mechanism.png
├── COMPREHENSIVE_REPORT.txt
├── individual_combinations/          ← NEW: 45 detailed combo plots
│   ├── combo_low_0.010_0.33.png
│   ├── combo_low_0.010_0.66.png
│   ├── combo_low_0.010_1.00.png
│   ├── ... (42 more)
│   └── combo_severe_0.050_1.00.png
└── time_series_all/                  ← NEW: 45 time series plots
    ├── ts_low_0.010_0.33.png
    ├── ts_low_0.010_0.66.png
    ├── ts_low_0.010_1.00.png
    ├── ... (42 more)
    └── ts_severe_0.050_1.00.png
```

---

## Plot Categories

### 1. Summary Plots (19 plots)

These show aggregated results across scenarios:

| File | Description | Scenarios Shown |
|------|-------------|-----------------|
| **01** | Budget comparison across all climates | Medium, middle renewal, all climates |
| **02** | Net budget efficiency by renewal rate | All scenarios, middle climate |
| **03** | Prosperity (E, T) with baselines | Medium, middle renewal, all climates |
| **04** | Strong resilience with baseline | All scenarios, middle climate |
| **05** | Armed groups with baseline | All scenarios, middle climate |
| **06** | Adaptive strategy performance | Medium, middle climate/renewal |
| **07** | Strategy effectiveness heatmap | All scenarios, middle climate/renewal |
| **0A** | Budget renewal comparison (comprehensive) | All scenarios, all renewals, middle climate |
| **0B** | Best strategy by renewal rate | All scenarios, all renewals |
| **0C** | Detailed metrics grid | All combinations |
| **0D** | Overall winner analysis | All combinations |
| **0E** | Win rate analysis | All combinations |
| **0F** | Time series across ALL budgets | Medium, middle climate, all renewals |
| **0G** | Time series across ALL scenarios | All scenarios, middle climate/renewal |
| **0H** | Budget evolution time series | Medium, middle climate/renewal |
| **0I** | Process mode evolution | Representative (middle combo) |
| **0J** | Multi-scale mode distribution | Representative (middle combo) |
| **0K** | Cognitive domain fraction | Representative (middle combo) |
| **0L** | Seasonal variation effects | Representative (middle combo) |
| **0M** | Steady state mechanism | Representative (middle combo) |

### 2. Individual Combination Plots (45 plots) - **NEW**

**Directory**: `individual_combinations/`

Each plot shows a **3×3 grid** of metrics for one specific combination:

**Metrics shown**:
1. Final Armed Groups (AG) - red bars
2. Final Strong Resilient (SR) - green bars
3. Final Economy (E) - gold bars
4. Final Trust (T) - teal bars
5. Final Governance - blue bars
6. Final Violence (V) - dark red bars
7. Total Budget Consumed - purple bars
8. Total Budget Generated - cyan bars
9. Net Budget (Gen - Cons) - navy bars

**For each metric**:
- All 12 strategy combinations shown (4 gov × 3 AG)
- Adaptive strategies highlighted in dark green
- Best performer: lime border
- Worst performer: red border
- Values labeled on each bar

**Filename format**: `combo_{scenario}_{climate}_{renewal}.png`

**Example**: `combo_medium_0.030_0.66.png`
- Scenario: medium conflict
- Climate amplitude: 0.030
- Budget renewal: 0.66

### 3. Time Series Plots (45 plots) - **NEW**

**Directory**: `time_series_all/`

Each plot shows a **2×3 grid** of temporal evolution for one specific combination:

**Variables shown**:
1. Armed Groups (AG) over time - red lines
2. Strong Resilient (SR) over time - green lines
3. Economy (E) over time - gold lines
4. Trust (T) over time - teal lines
5. Violence (V) over time - dark red lines
6. Governance over time - blue lines

**For each variable**:
- All 4 government strategies plotted
- Adaptive strategy: thicker line (3px vs 2px)
- Initial baseline: black dashed line
- Averaged over 3 AG strategies for clarity

**Filename format**: `ts_{scenario}_{climate}_{renewal}.png`

**Example**: `ts_severe_0.050_1.00.png`
- Scenario: severe conflict
- Climate amplitude: 0.050
- Budget renewal: 1.00

---

## Complete List of Individual Plots

### Low Conflict (15 plots)

```
individual_combinations/
├── combo_low_0.010_0.33.png    time_series_all/├── ts_low_0.010_0.33.png
├── combo_low_0.010_0.66.png    ├── ts_low_0.010_0.66.png
├── combo_low_0.010_1.00.png    ├── ts_low_0.010_1.00.png
├── combo_low_0.020_0.33.png    ├── ts_low_0.020_0.33.png
├── combo_low_0.020_0.66.png    ├── ts_low_0.020_0.66.png
├── combo_low_0.020_1.00.png    ├── ts_low_0.020_1.00.png
├── combo_low_0.030_0.33.png    ├── ts_low_0.030_0.33.png
├── combo_low_0.030_0.66.png    ├── ts_low_0.030_0.66.png
├── combo_low_0.030_1.00.png    ├── ts_low_0.030_1.00.png
├── combo_low_0.040_0.33.png    ├── ts_low_0.040_0.33.png
├── combo_low_0.040_0.66.png    ├── ts_low_0.040_0.66.png
├── combo_low_0.040_1.00.png    ├── ts_low_0.040_1.00.png
├── combo_low_0.050_0.33.png    ├── ts_low_0.050_0.33.png
├── combo_low_0.050_0.66.png    ├── ts_low_0.050_0.66.png
└── combo_low_0.050_1.00.png    └── ts_low_0.050_1.00.png
```

### Medium Conflict (15 plots)

```
individual_combinations/
├── combo_medium_0.010_0.33.png    time_series_all/├── ts_medium_0.010_0.33.png
├── combo_medium_0.010_0.66.png    ├── ts_medium_0.010_0.66.png
├── combo_medium_0.010_1.00.png    ├── ts_medium_0.010_1.00.png
├── combo_medium_0.020_0.33.png    ├── ts_medium_0.020_0.33.png
├── combo_medium_0.020_0.66.png    ├── ts_medium_0.020_0.66.png
├── combo_medium_0.020_1.00.png    ├── ts_medium_0.020_1.00.png
├── combo_medium_0.030_0.33.png    ├── ts_medium_0.030_0.33.png
├── combo_medium_0.030_0.66.png    ├── ts_medium_0.030_0.66.png
├── combo_medium_0.030_1.00.png    ├── ts_medium_0.030_1.00.png
├── combo_medium_0.040_0.33.png    ├── ts_medium_0.040_0.33.png
├── combo_medium_0.040_0.66.png    ├── ts_medium_0.040_0.66.png
├── combo_medium_0.040_1.00.png    ├── ts_medium_0.040_1.00.png
├── combo_medium_0.050_0.33.png    ├── ts_medium_0.050_0.33.png
├── combo_medium_0.050_0.66.png    ├── ts_medium_0.050_0.66.png
└── combo_medium_0.050_1.00.png    └── ts_medium_0.050_1.00.png
```

### Severe Conflict (15 plots)

```
individual_combinations/
├── combo_severe_0.010_0.33.png    time_series_all/├── ts_severe_0.010_0.33.png
├── combo_severe_0.010_0.66.png    ├── ts_severe_0.010_0.66.png
├── combo_severe_0.010_1.00.png    ├── ts_severe_0.010_1.00.png
├── combo_severe_0.020_0.33.png    ├── ts_severe_0.020_0.33.png
├── combo_severe_0.020_0.66.png    ├── ts_severe_0.020_0.66.png
├── combo_severe_0.020_1.00.png    ├── ts_severe_0.020_1.00.png
├── combo_severe_0.030_0.33.png    ├── ts_severe_0.030_0.33.png
├── combo_severe_0.030_0.66.png    ├── ts_severe_0.030_0.66.png
├── combo_severe_0.030_1.00.png    ├── ts_severe_0.030_1.00.png
├── combo_severe_0.040_0.33.png    ├── ts_severe_0.040_0.33.png
├── combo_severe_0.040_0.66.png    ├── ts_severe_0.040_0.66.png
├── combo_severe_0.040_1.00.png    ├── ts_severe_0.040_1.00.png
├── combo_severe_0.050_0.33.png    ├── ts_severe_0.050_0.33.png
├── combo_severe_0.050_0.66.png    ├── ts_severe_0.050_0.66.png
└── combo_severe_0.050_1.00.png    └── ts_severe_0.050_1.00.png
```

---

## Total Plot Count

| Category | Count | Description |
|----------|-------|-------------|
| Summary plots | 19 | Aggregated analysis across scenarios |
| Individual combo plots | 45 | Bar charts for each combination |
| Time series plots | 45 | Temporal evolution for each combination |
| **TOTAL** | **109** | **Complete visualization suite** |

---

## Usage for Article

### Figure Selection Strategy

1. **Main Results**: Use summary plots (01-0E) for overall findings
2. **Specific Comparisons**: Use individual combo plots for detailed analysis
3. **Temporal Dynamics**: Use time series plots to show evolution
4. **Supplementary Materials**: Include all 45 combo plots as appendix

### Recommended Figure Layout for Article

**Main Text Figures**:
- Figure 1: Plot 0E (Win rate analysis) - Shows adaptive superiority
- Figure 2: Plot 0A (Budget renewal comparison) - All scenarios
- Figure 3: Plot 0G (Time series all scenarios) - Temporal dynamics
- Figure 4: Plot 0L (Seasonal variation) - Environmental effects
- Figure 5: Plot 0I (Process mode evolution) - Behavioral analysis

**Supplementary Figures**:
- S1-S15: Low conflict individual combos
- S16-S30: Medium conflict individual combos
- S31-S45: Severe conflict individual combos
- S46-S90: All time series plots

### How to Identify Best Scenarios

From individual combo plots:
- **Green bars with lime border**: Adaptive strategy winning
- **Highest SR, lowest AG**: Best outcomes
- **Lowest budget consumption**: Most efficient

From time series plots:
- **Green line (adaptive) above others for SR**: Sustained success
- **Red line (adaptive) below others for AG**: Conflict reduction
- **Stable vs oscillating**: Dynamic stability analysis

---

## Key Improvements from Previous Version

### Before
- Only ~20 plots generated
- Many combinations not visualized
- Had to infer missing data
- Selected "representative" scenarios only
- Some plots showed only middle climate/renewal

### After
- **109 plots** generated
- **Every combination visualized**
- Complete data coverage
- All scenarios, climates, renewals shown
- No gaps in analysis

### Specific Enhancements

1. **Matplotlib Backend**: Set to 'Agg' (non-interactive)
   - Plots save directly without display windows
   - Faster execution
   - No manual closing required

2. **Individual Combination Plots**: New feature
   - 3×3 metric grid for each combo
   - All strategy combinations visible
   - Color-coded for quick identification

3. **Comprehensive Time Series**: New feature
   - 2×3 temporal grid for each combo
   - All government strategies compared
   - Baseline reference lines included

4. **Model File**: Updated to model3.txt
   - Consistent with diagnostic scripts
   - Uses latest model specification

---

## Verification Checklist

After running the script, verify:

✅ **Summary directory exists**:
```bash
ls visualizations/conflict_simulation_results_adaptive/
```

✅ **19 summary plots created**:
```bash
ls visualizations/conflict_simulation_results_adaptive/*.png | wc -l
# Should output: 19
```

✅ **45 individual combo plots created**:
```bash
ls visualizations/conflict_simulation_results_adaptive/individual_combinations/*.png | wc -l
# Should output: 45
```

✅ **45 time series plots created**:
```bash
ls visualizations/conflict_simulation_results_adaptive/time_series_all/*.png | wc -l
# Should output: 45
```

✅ **Total plot count**:
```bash
find visualizations/conflict_simulation_results_adaptive/ -name "*.png" | wc -l
# Should output: 109
```

✅ **Report exists**:
```bash
ls visualizations/conflict_simulation_results_adaptive/COMPREHENSIVE_REPORT.txt
```

---

## Execution Time Estimate

Based on 540 simulations + plotting:
- Simulations: ~15-20 minutes (depends on CPU)
- Summary plots: ~2-3 minutes
- Individual combo plots: ~3-5 minutes
- Time series plots: ~3-5 minutes
- **Total**: ~25-35 minutes

---

## Troubleshooting

### Issue: "No plots generated"
**Solution**: Check save_dir parameter and permissions

### Issue: "Missing combinations"
**Solution**: Verify CLIMATE_PARAMETERS and BUDGET_RENEWAL_RATES lists

### Issue: "Plots displayed during run"
**Solution**: Ensure `matplotlib.use('Agg')` is before `import matplotlib.pyplot`

### Issue: "Memory error"
**Solution**: Reduce number of climate parameters or run in batches

---

## Summary

✅ **All 45 scenario combinations fully visualized**
✅ **Non-interactive plotting (no display windows)**
✅ **Comprehensive coverage for article figures**
✅ **Ready for publication-quality analysis**

Run the script and you'll have a complete visualization suite covering every possible combination of parameters!

---

**End of Guide**
