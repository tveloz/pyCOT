# CLIMATE_PARAMETERS Update - Seasonal Amplitude Implementation

**Date**: 2026-01-27
**File Modified**: [script_two_budget_adaptive.py](scripts/XCEPT_project/Article_1/script_two_budget_adaptive.py)
**Status**: ✅ IMPLEMENTED AND VERIFIED

---

## Summary

CLIMATE_PARAMETERS now control the **seasonal amplitude** of resource renewal reactions (r1-r7), using the formula:

```
rate = (A/2) * (1 + cos(2π*t/T))
```

where **A** is the climate parameter (amplitude).

This produces rates that vary from **[0, A]** over the year, with an annual average of **A/2**.

---

## What Changed

### BEFORE: CLIMATE_PARAMETERS as Base Kinetic Constants

**Previous behavior**:
- CLIMATE_PARAMETERS = [0.01, 0.02, 0.03, 0.04, 0.05] were used as **base rate constants**
- Seasonal variation was an **additive modulation**: `rate = k * (1 + 0.3*cos(...))`
- Range: `rate ∈ [0.7k, 1.3k]` (±30% variation)
- Only applied to r1, r2, r9, r10

**Problems**:
- Climate parameters affected base dynamics continuously
- Seasonal variation was a fixed percentage (30%)
- Not aligned with the concept of seasonal resource renewal

### AFTER: CLIMATE_PARAMETERS as Seasonal Amplitudes

**Current behavior**:
- CLIMATE_PARAMETERS = [0.01, 0.02, 0.03, 0.04, 0.05] are **seasonal amplitudes** (A)
- Seasonal variation is **multiplicative**: `rate = (A/2) * (1 + cos(...))`
- Range: `rate ∈ [0, A]` (from zero to maximum)
- Applies to all resource renewal/degradation reactions (r1-r7)

**Benefits**:
- Climate parameters directly control resource renewal capacity
- Seasonal pattern represents realistic annual cycles (e.g., rainy/dry seasons)
- Annual average resource renewal = A/2
- Studying different scenarios now means comparing different resource renewal amplitudes

---

## Implementation Details

### 1. Kinetic Constants (Lines 352-378)

**BEFORE**:
```python
self.k = {
    "r1": 0.02, "r2": 0.02, "r3": 0.02, "r4": 0.02, "r5": 0.02,
    "r6": 0.02, "r7": 0.05, ...
}
```

**AFTER**:
```python
self.k = {
    # Resource renewal/degradation reactions (r1-r7) - base rate = 1.0
    # Actual rate = seasonal_modifier * base_rate
    "r1": 1.0, "r2": 1.0, "r3": 1.0, "r4": 1.0, "r5": 1.0,
    "r6": 1.0, "r7": 1.0,
    # Other reactions maintain their original kinetic constants
    "r8": 0.002, "r9": 0.05, "r10": 0.03, ...
}
```

### 2. Seasonal Modifier Formula (Lines 209-254)

**BEFORE**:
```python
def get_seasonal_modifier(self, reaction: str) -> float:
    """Formula: 1 + amplitude * cos(2π * t / period + phase)"""
    seasonal_reactions = ["r1", "r2", "r9", "r10"]

    if reaction not in seasonal_reactions:
        return 1.0

    seasonal_factor = 1.0 + self.seasonal_amplitude * math.cos(
        2.0 * math.pi * self.time / self.seasonal_period + self.seasonal_phase
    )

    return seasonal_factor
```

**AFTER**:
```python
def get_seasonal_modifier(self, reaction: str) -> float:
    """
    Formula: (A/2) * (1 + cos(2π * t / period + phase))

    Ranges from [0, A] where A is the seasonal amplitude.
    """
    # All resource renewal/recovery reactions (r1-r7)
    seasonal_reactions = ["r1", "r2", "r3", "r4", "r5", "r6", "r7"]

    if reaction not in seasonal_reactions:
        return 1.0  # Non-seasonal reactions use constant rate

    seasonal_factor = (self.seasonal_amplitude / 2.0) * (
        1.0 + math.cos(2.0 * math.pi * self.time / self.seasonal_period + self.seasonal_phase)
    )

    return seasonal_factor
```

### 3. Climate Parameter Setter (Lines 263-277)

**BEFORE**:
```python
def set_climate_parameter(self, climate_k: float):
    """Set climate-related kinetic parameters (r1-r7)."""
    self.climate_k = climate_k
    self.k["r1"] = climate_k
    self.k["r2"] = climate_k
    self.k["r3"] = climate_k
    self.k["r4"] = climate_k
    self.k["r5"] = climate_k
    self.k["r6"] = climate_k
    self.k["r7"] = climate_k * 2
```

**AFTER**:
```python
def set_climate_parameter(self, climate_amplitude: float):
    """
    Set seasonal amplitude for resource renewal reactions (r1-r7).

    The climate parameter controls amplitude A in the formula:
    rate = (A/2) * (1 + cos(2π*t/T))

    This produces rates from 0 to A over the year, with average A/2.
    """
    self.climate_k = climate_amplitude  # Store for reference
    self.seasonal_amplitude = climate_amplitude  # Set the seasonal amplitude
```

---

## Reactions Affected

All resource renewal and degradation reactions now use seasonal variation:

| Reaction | Description | Formula |
|----------|-------------|---------|
| **r1** | SL → RL | Land restoration (stressed → restored) |
| **r2** | SR_SL → SR_RL | Strong resilient renewal |
| **r3** | WR_SL → WR_RL | Weak resilient renewal |
| **r4** | RL → SL | Land degradation (restored → stressed) |
| **r5** | SR_RL → SR_SL | Strong resilient degradation |
| **r6** | WR_RL → WR_SL | Weak resilient degradation |
| **r7** | AG_RL → AG_SL | Armed groups degradation |

**Common formula**: `rate = k * (A/2) * (1 + cos(2π*t/T))`

Where:
- k = 1.0 (base constant from set_kinetic_parameters)
- A = climate amplitude (from CLIMATE_PARAMETERS)
- T = 1.0 year (seasonal period)
- t = current simulation time

---

## Verification Results

Test script: [test_seasonal_climate.py](scripts/XCEPT_project/Article_1/test_seasonal_climate.py)

**Results** for all r1-r7 reactions across all CLIMATE_PARAMETERS:

| Climate A | Min Rate | Max Rate | Avg Rate | Status |
|-----------|----------|----------|----------|--------|
| 0.01 | ~0.000 | 0.010 | 0.0050 | ✓ |
| 0.02 | ~0.000 | 0.020 | 0.0101 | ✓ |
| 0.03 | ~0.000 | 0.030 | 0.0152 | ✓ |
| 0.04 | ~0.000 | 0.040 | 0.0202 | ✓ |
| 0.05 | ~0.000 | 0.050 | 0.0253 | ✓ |

**Observations**:
- ✅ Maximum rate matches A exactly
- ✅ Minimum rate approaches 0 (within numerical precision)
- ✅ Average rate ≈ A/2 (within 1% due to discretization)

**Visualization**: See [seasonal_climate_amplitude_test.png](scripts/XCEPT_project/Article_1/visualizations/seasonal_test/seasonal_climate_amplitude_test.png)

---

## How to Use

### Setting Climate Parameters

```python
# Create simulator
sim = ConflictSimulator(model_file)

# Set climate amplitude (controls seasonal resource renewal)
climate_amplitude = 0.03  # From CLIMATE_PARAMETERS list
sim.set_climate_parameter(climate_amplitude)

# Run simulation
sim.simulate(years=10, dt=0.01)
```

### Interpreting Results

For a given climate amplitude **A**:

- **Maximum resource renewal rate** = A (occurs when cos = +1, e.g., peak growing season)
- **Minimum resource renewal rate** = 0 (occurs when cos = -1, e.g., dry season)
- **Average resource renewal rate** = A/2 (over one year)

**Example**: With A = 0.04
- Peak season: Land restoration rate = 0.04/year
- Dry season: Land restoration rate = 0.00/year
- Annual average: Land restoration rate = 0.02/year

### Comparing Scenarios

Different CLIMATE_PARAMETERS now represent different **environmental conditions**:

- **A = 0.01**: Low resource renewal capacity (harsh environment)
- **A = 0.02**: Below-average resource renewal
- **A = 0.03**: Average resource renewal
- **A = 0.04**: Above-average resource renewal
- **A = 0.05**: High resource renewal capacity (favorable environment)

---

## Batch Simulations

The batch simulation automatically tests all CLIMATE_PARAMETERS:

```python
# From script_two_budget_adaptive.py
CLIMATE_PARAMETERS = [0.01, 0.02, 0.03, 0.04, 0.05]

# Batch run tests:
# - 3 conflict scenarios (low, medium, severe)
# - 5 climate amplitudes (0.01 to 0.05)
# - 3 budget renewal rates (0.33, 0.66, 1.0)
# - 4 government strategies (development, security, balanced, adaptive)
# - 3 armed group strategies (recruitment, displacement, balanced)
# Total: 3 × 5 × 3 × 4 × 3 = 540 simulations
```

Results will show how different resource renewal capacities affect conflict dynamics.

---

## Physical Interpretation

### Before (kinetic constants)
"Climate parameter k = 0.02 means land restores at rate 0.02/year with ±30% seasonal variation"
→ Rate varies from 0.014 to 0.026/year

### After (seasonal amplitudes)
"Climate amplitude A = 0.02 means land restoration varies seasonally from 0 to 0.02/year"
→ Rate varies from 0 to 0.02/year, averaging 0.01/year

**Key difference**: The new formulation represents **seasonal availability** of resources rather than continuous processes with seasonal modulation.

---

## Recommendations for Analysis

1. **Compare across amplitudes**: How does resource renewal capacity affect long-term stability?
2. **Examine seasonal effects**: Do strategies perform differently in high vs low resource periods?
3. **Identify critical thresholds**: Is there a minimum A below which conflict becomes intractable?
4. **Budget efficiency**: Does higher A allow lower budget expenditure?

---

## Backward Compatibility

⚠️ **Breaking change**: Results from simulations using the old climate parameter interpretation cannot be directly compared to new results.

**Conversion**: If you want to approximate old behavior, you would need to:
1. Set CLIMATE_PARAMETERS to 2× the old values (to match average rate)
2. Note that the seasonal pattern is fundamentally different

**Recommendation**: Re-run all analyses with the new implementation for consistency.

---

## Files Modified

1. ✅ [script_two_budget_adaptive.py](scripts/XCEPT_project/Article_1/script_two_budget_adaptive.py)
   - Lines 352-378: Updated set_kinetic_parameters (r1-r7 = 1.0)
   - Lines 209-254: Updated get_seasonal_modifier formula
   - Lines 263-277: Updated set_climate_parameter to set amplitude

2. ✅ [test_seasonal_climate.py](scripts/XCEPT_project/Article_1/test_seasonal_climate.py) - Verification script

---

## Summary Table

| Aspect | Before | After |
|--------|--------|-------|
| **Formula** | 1 + 0.3*cos(...) | (A/2)*(1 + cos(...)) |
| **Range** | [0.7, 1.3] (relative) | [0, A] (absolute) |
| **Parameters** | Fixed amplitude (0.3) | Variable amplitude (A from CLIMATE_PARAMETERS) |
| **Reactions** | r1, r2, r9, r10 | r1, r2, r3, r4, r5, r6, r7 |
| **Interpretation** | ±30% variation | Seasonal resource availability |
| **Average** | 1.0 (100% baseline) | A/2 |

---

## Conclusion

✅ **CLIMATE_PARAMETERS now properly control seasonal amplitude for resource renewal**

The implementation correctly transforms CLIMATE_PARAMETERS from base kinetic constants to seasonal amplitudes, providing a more realistic representation of environmental variability in resource availability.

**Ready for production use** - All batch simulations will now explore different environmental resource renewal capacities.

---

**End of Documentation**
