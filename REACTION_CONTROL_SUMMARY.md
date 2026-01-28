# Reaction Control Mechanism - Summary

## Your Question Was Correct

You asked why r31 (V decay), r26 (Gov decay), and r25 (trust destruction) become inactive when they should have autonomous kinetics. **You are absolutely right to be confused** - the implementation has a critical bug.

## Evidence from Simulation Output

From the steady state analysis of the security strategy:

```
V (STEADY STATE):
  Mean: 15.770, Std: 0.000
  Production: 0.0000
  Consumption: 0.0000
  Inactive consumers: r15, r16, r25, r31  ← r31 SHOULD BE ACTIVE!

Gov (STEADY STATE):
  Mean: 19.098, Std: 0.000
  Production: 0.0000
  Consumption: 0.0000
  Inactive consumers: r8, r21, r24, r29, r30, r26, r36  ← r26 SHOULD BE ACTIVE!
```

With V = 15.77 and k31 = 0.1, the rate of r31 should be:
**rate = 0.1 × 15.77 = 1.577** (definitely not zero!)

With Gov = 19.10 and k26 = 0.02, the rate of r26 should be:
**rate = 0.02 × 19.10 = 0.382** (definitely not zero!)

## Complete Classification of Reactions

### CATEGORY A: Pure Mass Action (Autonomous - Should Always Be Active)

| Reaction | Formula | Implementation | Rate Formula | Status |
|----------|---------|----------------|--------------|--------|
| r1  | SL => RL | `rate *= SL` | k × SL | ✓ Correct |
| r2  | SR_SL => SR_RL | `rate *= SR_SL` | k × SR_SL | ✓ Correct |
| r3  | WR_SL => WR_RL | `rate *= WR_SL` | k × WR_SL | ✓ Correct |
| r4  | RL => SL | `rate *= RL` | k × RL | ✓ Correct |
| r5  | SR_RL => SR_SL | `rate *= SR_RL` | k × SR_RL | ✓ Correct |
| r6  | WR_RL => WR_SL | `rate *= WR_RL` | k × WR_RL | ✓ Correct |
| r7  | AG_RL => AG_SL | `rate *= AG_RL` | k × AG_RL | ✓ Correct |
| r9  | SR_RL => SR_SL + E | `rate *= SR_RL` | k × SR_RL | ✓ Correct |
| r10 | WR_RL => WR_SL + E | `rate *= WR_RL` | k × WR_RL | ✓ Correct |
| r12 | => Gov | `pass` (rate = k) | k (constant) | ✓ Correct |
| r13 | SR_SL + RL => ... | `rate *= SR_SL × RL` | k × SR_SL × RL | ✓ Correct |
| r14 | WR_SL + RL => ... | `rate *= WR_SL × RL` | k × WR_SL × RL | ✓ Correct |
| r17 | AG_SL + RL => ... | `rate *= AG_SL × RL` | k × AG_SL × RL | ✓ Correct |
| r20 | SR_SL => WR_SL | `rate *= SR_SL` | k × SR_SL | ✓ Correct |
| **r26** | **Gov => ;** | `rate *= Gov` | k × Gov | ✓ Correct implementation, **but inactive in security strategy!** |
| r27 | WR_SL => WR_SL + V | `rate *= WR_SL` | k × WR_SL | ✓ Correct |
| r28 | AG_SL + SR_RL => ... | `rate *= AG_SL × SR_RL` | k × AG_SL × SR_RL | ✓ Correct |
| **r31** | **V => ;** | `rate *= V` | k × V | ✓ Correct implementation, **but inactive in security strategy!** |
| r33 | AG_SL + WR_RL => ... | `rate *= AG_SL × WR_RL` | k × AG_SL × WR_RL | ✓ Correct |
| r34 | AG_SL + SR_RL => ... | `rate *= AG_SL × SR_RL` | k × AG_SL × SR_RL | ✓ Correct |
| r35 | AG_SL + E => AG_RL | `rate *= AG_SL × E` | k × AG_SL × E | ✓ Correct (E consumed manually) |
| r36 | AG_RL + Gov => AG_RL | `rate *= AG_RL × Gov` | k × AG_RL × Gov | ✓ Correct (Gov consumed manually) |

### CATEGORY B: Budget-Controlled by Government (E or Gov Allocation)

| Reaction | Formula | Implementation | Rate Formula | Status |
|----------|---------|----------------|--------------|--------|
| r8  | Gov + SR_SL => SR_RL | `rate *= Gov_alloc × SR_SL` | k × Gov_alloc["r8"] × SR_SL | ✓ Correct |
| r11 | E => Gov | `rate *= E_alloc` | k × E_alloc["r11"] | ✓ Correct |
| r18 | E + WR_RL => SR_RL | `rate *= E_alloc × WR_RL` | k × E_alloc["r18"] × WR_RL | ✓ Correct |
| r19 | E + WR_SL => SR_SL | `rate *= E_alloc × WR_SL` | k × E_alloc["r19"] × WR_SL | ✓ Correct |
| r21 | E + Gov + AG_SL => WR_SL | `rate *= sqrt(E_alloc × Gov_alloc) × AG_SL` | Geometric mean | ✓ Correct |
| r22 | SR_RL + E => SR_RL + T | `rate *= SR_RL × E_alloc` | k × SR_RL × E_alloc["r22"] | ✓ Correct |
| r23 | Gov + SR_RL => ... + T | `rate *= Gov_alloc × SR_RL` | k × Gov_alloc["r23"] × SR_RL | ✓ Correct |
| r24 | Gov + WR_RL => SR_RL | `rate *= Gov_alloc × WR_RL` | k × Gov_alloc["r24"] × WR_RL | ✓ Correct |
| r29 | AG_SL + Gov => ... + V | `rate *= AG_SL × Gov_alloc` | k × AG_SL × Gov_alloc["r29"] | ✓ Correct |
| r30 | AG_RL + Gov => ... + V | `rate *= AG_RL × Gov_alloc` | k × AG_RL × Gov_alloc["r30"] | ✓ Correct |

### CATEGORY C: Budget-Controlled by Armed Groups (V Allocation)

| Reaction | Formula | Implementation | Rate Formula | Status |
|----------|---------|----------------|--------------|--------|
| r15 | SR_RL + SL + V => ... | `rate *= SR_RL × SL × V_alloc` | k × SR_RL × SL × V_alloc["r15"] | ✓ Correct |
| r16 | WR_RL + SL + V => ... | `rate *= WR_RL × SL × V_alloc` | k × WR_RL × SL × V_alloc["r16"] | ✓ Correct |
| **r25** | **V + T => ;** | `rate *= V_alloc × T` | k × V_alloc["r25"] × T | **✗ BUG! Should be pure mass action!** |
| r32 | WR_SL + AG_SL + V => ... | `rate *= WR_SL × AG_SL × V_alloc` | k × WR_SL × AG_SL × V_alloc["r32"] | ✓ Correct |

## The Critical Bug: r25

**MODEL SPECIFICATION**: `r25: V + T => ;`

This is a standard bimolecular reaction showing automatic trust destruction by violence.

**EXPECTED IMPLEMENTATION** (pure mass action):
```python
elif reaction == "r25":
    rate *= state.get("V", 0) * state.get("T", 0)
    # No budget allocation needed
```

**CURRENT IMPLEMENTATION** (incorrectly budget-controlled):
```python
elif reaction == "r25":
    v_used = V_alloc.get("r25", 0)
    rate *= v_used * state.get("T", 0)
    consumption["V"] = rate
```

This makes r25 dependent on armed groups' strategic allocation, which is conceptually wrong.

## Saturation Mechanism: Rate Limiter

The only saturation mechanism is the rate limiter (lines 725-739):

```python
# Rate limiter: prevent excessive consumption
for i, species in enumerate(self.SPECIES):
    current = self.state.get(species, 0)
    if delta_raw[i] < 0 and current > 0:
        max_consumption = current * 0.5  # Max 50% per timestep
        if abs(delta_raw[i]) > max_consumption:
            species_scale = max_consumption / abs(delta_raw[i])
            scale_factor = min(scale_factor, species_scale)

if scale_factor < 1.0:
    process_vector = process_vector * scale_factor
```

This prevents any species from losing more than 50% of its stock in a single timestep. It scales down ALL reactions proportionally when this limit is exceeded.

## Why r31 and r26 Are Inactive (The Mystery)

### The Evidence

1. r31 and r26 are correctly implemented as pure mass action
2. Both have non-zero kinetic constants (k31=0.1, k26=0.02)
3. Both have non-zero reactants (V=15.77, Gov=19.10)
4. Neither is budget-controlled
5. Yet both show rate = 0 in the security strategy steady state

### The Investigation

Looking at the state update code (lines 759-762):

```python
v_production = sum(process_vector[i] * dt for i, rxn in enumerate(self.REACTIONS)
                   if rxn in ["r27", "r28", "r29", "r30"])
r31_idx = self.REACTIONS.index("r31")
v_decay = process_vector[r31_idx] * dt
new_state["V"] = max(0.0, old_state["V"] + v_production - v_decay - total_consumption["V"] * dt)
```

The V update formula is:
**V(new) = V(old) + production - decay - consumption**

For V to be exactly constant (Std = 0.000):
**production = decay + consumption**

But the analysis shows: production = 0, decay = 0, consumption = 0

### Possible Explanations

1. **Numerical Freezing**: After AG_SL reaches exactly 0.0, a cascade of reactions become inactive. Through some numerical artifact, this propagates to make ALL reactions inactive, including r31 and r26.

2. **Implementation Bug**: There may be additional code not visible in the sections I examined that explicitly sets reactions to zero under certain conditions.

3. **Rate History Recording**: The rates might be computed correctly but recorded incorrectly in rate_history.

4. **Perfect Balance Roundoff**: The system reaches such a perfect balance that floating-point roundoff causes rates to be stored as exactly zero.

### What We Know For Sure

- **r26 implementation is correct**: It should decay governance at rate 0.02 × Gov
- **r31 implementation is correct**: It should decay violence at rate 0.1 × V
- **r25 implementation is WRONG**: It should use mass action, not budget allocation
- **Something else is making r26 and r31 inactive** in the security strategy

## Recommended Actions

### 1. Fix r25 (Confirmed Bug)

**File**: [script_two_budget_adaptive.py](script_two_budget_adaptive.py:638-641)

```python
# BEFORE (lines 638-641)
elif reaction == "r25":
    v_used = V_alloc.get("r25", 0)
    rate *= v_used * state.get("T", 0)
    consumption["V"] = rate

# AFTER
elif reaction == "r25":
    rate *= state.get("V", 0) * state.get("T", 0)
    # No consumption tracking - V consumed through stoichiometry only
```

Also update V_weights to remove r25 and renormalize (lines 527-531):

```python
# BEFORE
if strategy == "recruitment":
    self.V_weights = {"r15": 0.1, "r16": 0.1, "r25": 0.2, "r32": 0.6}
elif strategy == "displacement":
    self.V_weights = {"r15": 0.4, "r16": 0.4, "r25": 0.1, "r32": 0.1}
elif strategy == "balanced":
    self.V_weights = {"r15": 0.25, "r16": 0.25, "r25": 0.25, "r32": 0.25}

# AFTER (renormalized)
if strategy == "recruitment":
    self.V_weights = {"r15": 0.125, "r16": 0.125, "r32": 0.75}
elif strategy == "displacement":
    self.V_weights = {"r15": 0.5, "r16": 0.5, "r32": 0.0}  # Or r32: 0.0, redistributed
elif strategy == "balanced":
    self.V_weights = {"r15": 0.33, "r16": 0.33, "r32": 0.34}
```

### 2. Investigate r31 and r26 Inactivity

Add diagnostic output in [apply_timestep](script_two_budget_adaptive.py:718) to track:
- process_vector values for r31 and r26
- scale_factor from rate limiter
- Actual state values when rates become zero

### 3. Check for Hidden Code

Search for any code that might explicitly zero out reactions:
- Conditional statements that set rates to zero
- Threshold checks (e.g., "if rate < 1e-10: rate = 0")
- Special handling for "steady state" detection

## Summary Table

| Reaction | Type | Implemented Correctly? | Actually Works? | Notes |
|----------|------|----------------------|-----------------|-------|
| r1-r20 | Mass Action | ✓ Yes | ✓ Yes | - |
| r21 | Budget (geometric mean) | ✓ Yes | ✓ Yes | Special handling |
| r22-r24 | Budget (Gov/E) | ✓ Yes | ✓ Yes | - |
| **r25** | **Should be Mass Action** | **✗ No** | **✗ No** | **BUG: Uses V allocation** |
| **r26** | **Mass Action** | **✓ Yes** | **✗ No** | **Mystery: Becomes inactive** |
| r27-r30 | Mass Action / Budget | ✓ Yes | ✓ Yes | - |
| **r31** | **Mass Action** | **✓ Yes** | **✗ No** | **Mystery: Becomes inactive** |
| r32-r36 | Budget / Mass Action | ✓ Yes | ✓ Yes | - |

## Kinetic Constants

```python
k25 = 0.002  # Trust destruction rate
k26 = 0.02   # Governance decay rate
k31 = 0.1    # Violence decay rate (largest of the three!)
```

With V=15.77 and k31=0.1, r31 should have rate = 1.577 per time unit, which with dt=0.01 gives a change of 0.01577 per timestep. This is 0.1% of V, well within the 50% rate limiter threshold.

**There is no reason r31 should be inactive based on the code we can see.**
