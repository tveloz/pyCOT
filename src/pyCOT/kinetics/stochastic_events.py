"""
Stochastic Event-Based Kinetics

This module provides mechanisms for modeling random, discrete events in
reaction network simulations. Unlike deterministic kinetics that give continuous
rates, these methods generate discrete stochastic events with specified statistics.

Includes:
- Poisson events (rare shocks like floods, droughts, conflicts)
- Gillespie SSA (exact stochastic simulation)
- Tau-leaping (approximate, faster stochastic)
- Langevin dynamics (continuous noise, SDE)
- Markov regime switching

Use cases: Climate extremes, sudden policy changes, random interventions,
discrete resource shocks, catastrophic events.
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Callable, Any


class PoissonEventGenerator:
    """
    Generate rare, random events following a Poisson process.
    
    Events occur with average frequency λ (events per time unit).
    Waiting times between events follow exponential distribution.
    
    Use for modeling:
    - Climate extremes (floods, droughts, hurricanes)
    - Conflict shocks (sudden violence spikes, terrorist attacks)
    - Economic crises (market crashes, supply disruptions)
    - Policy interventions (sudden regulations, aid packages)
    
    Attributes:
    -----------
    event_rate : float
        Average number of events per time unit (λ)
    event_impact : callable or dict
        Function or dict defining how events affect species concentrations
        Signature: event_impact(concentrations, t) -> delta_concentrations
    target_species : list, optional
        Species affected by events (None = all species)
    magnitude_distribution : callable, optional
        Distribution for event magnitudes (default: constant)
    """
    
    def __init__(self, 
                 event_rate: float,
                 event_impact: Callable,
                 target_species: Optional[List[str]] = None,
                 magnitude_distribution: Optional[Callable] = None,
                 seed: Optional[int] = None):
        """
        Initialize Poisson event generator.
        
        Parameters:
        -----------
        event_rate : float
            Mean events per time unit (higher = more frequent)
            Examples: 0.1 (rare), 1.0 (moderate), 10.0 (frequent)
        event_impact : callable
            Function(concentrations, species_idx, magnitude) -> delta_x
            Returns change in concentrations due to event
        target_species : list, optional
            Names of species affected (None = all)
        magnitude_distribution : callable, optional
            Function() -> float for event intensity
            Default: constant magnitude 1.0
        seed : int, optional
            Random seed for reproducibility
        """
        self.event_rate = event_rate
        self.event_impact = event_impact
        self.target_species = target_species
        self.magnitude_distribution = magnitude_distribution or (lambda: 1.0)
        self.rng = np.random.default_rng(seed)
        
        # Event schedule (generated lazily)
        self.next_event_time = None
        self.event_history = []
    
    def generate_next_event_time(self, current_time: float) -> float:
        """
        Generate time of next event using exponential distribution.
        
        Waiting time ~ Exponential(λ), where λ = event_rate
        Mean waiting time = 1/λ
        """
        waiting_time = self.rng.exponential(1.0 / self.event_rate)
        return current_time + waiting_time
    
    def apply_event(self, 
                   concentrations: np.ndarray,
                   species_idx: Dict[str, int],
                   current_time: float) -> np.ndarray:
        """
        Apply event impact to current concentrations.
        
        Returns:
        --------
        delta_x : array
            Change in concentrations due to event
        """
        magnitude = self.magnitude_distribution()
        delta_x = self.event_impact(concentrations, species_idx, magnitude)
        
        # Record event
        self.event_history.append({
            'time': current_time,
            'magnitude': magnitude,
            'delta_x': delta_x.copy()
        })
        
        return delta_x
    
    def get_event_schedule(self, t_start: float, t_end: float) -> List[float]:
        """
        Pre-generate all event times in simulation window.
        
        Useful for event-driven simulation approach (Option C from strategy).
        """
        event_times = []
        current_time = t_start
        
        while current_time < t_end:
            current_time = self.generate_next_event_time(current_time)
            if current_time < t_end:
                event_times.append(current_time)
        
        return event_times


# Pre-defined event impact templates
def concentration_shock(magnitude: float, 
                       target_species_indices: List[int],
                       shock_type: str = 'additive'):
    """
    Factory for creating concentration shock impacts.
    
    Parameters:
    -----------
    magnitude : float
        Strength of shock
    target_species_indices : list
        Indices of affected species
    shock_type : str
        'additive': Δx = magnitude
        'multiplicative': Δx = magnitude * x
        'depletion': Δx = -magnitude * x
    
    Returns:
    --------
    callable
        Impact function compatible with PoissonEventGenerator
    """
    def impact_function(concentrations, species_idx, event_magnitude):
        delta_x = np.zeros_like(concentrations)
        
        for idx in target_species_indices:
            if shock_type == 'additive':
                delta_x[idx] = magnitude * event_magnitude
            elif shock_type == 'multiplicative':
                delta_x[idx] = magnitude * event_magnitude * concentrations[idx]
            elif shock_type == 'depletion':
                delta_x[idx] = -magnitude * event_magnitude * concentrations[idx]
        
        return delta_x
    
    return impact_function


def rate_perturbation_event(base_rate_params: Dict,
                            perturbation_factor: float,
                            duration: float):
    """
    Factory for temporary rate parameter perturbations.
    
    Events temporarily modify kinetic parameters (e.g., rate constants).
    Returns to baseline after duration.
    
    Parameters:
    -----------
    base_rate_params : dict
        Original kinetic parameters
    perturbation_factor : float
        Multiplier for parameters (2.0 = double rates)
    duration : float
        How long perturbation lasts
    
    Use case: Policy interventions that temporarily boost/suppress reactions.
    """
    # This requires integration with simulation loop
    # Returns event specification for simulation wrapper
    return {
        'type': 'rate_perturbation',
        'factor': perturbation_factor,
        'duration': duration,
        'base_params': base_rate_params
    }


class GillespieSSA:
    """
    Gillespie Stochastic Simulation Algorithm (exact method).
    
    Treats reactions as discrete, stochastic events rather than continuous flows.
    Exact for systems with discrete molecule counts and stochastic reaction events.
    
    More computationally expensive than ODE integration but captures
    fluctuations important in small population systems.
    
    Reference: Gillespie, D. T. (1977). J. Phys. Chem. 81(25), 2340–2361.
    """
    
    def __init__(self, rn, kinetic_params, seed=None):
        """
        Initialize Gillespie SSA simulator.
        
        Parameters:
        -----------
        rn : ReactionNetwork
            Network structure
        kinetic_params : dict
            Rate constants for each reaction
        seed : int, optional
            Random seed
        """
        self.rn = rn
        self.kinetic_params = kinetic_params
        self.rng = np.random.default_rng(seed)
        
    def propensity(self, x, reaction_idx):
        """
        Calculate propensity (instantaneous probability) of reaction firing.
        
        For mass action: a_j = k_j * Π(x_i choose α_ij)
        """
        # Simplified mass action propensity
        # Full implementation requires stoichiometry
        pass
    
    def step(self, x, t):
        """
        Single Gillespie step: choose reaction and advance time.
        
        Returns:
        --------
        x_new : array
            Updated state after one reaction
        t_new : float
            Time of next reaction
        reaction_fired : int
            Index of reaction that occurred
        """
        # Calculate all propensities
        propensities = np.array([self.propensity(x, i) for i in range(len(self.rn.reactions()))])
        a0 = propensities.sum()
        
        if a0 == 0:
            return x, np.inf, None  # No reactions possible
        
        # Time to next reaction (exponential)
        tau = self.rng.exponential(1.0 / a0)
        
        # Choose reaction (weighted by propensities)
        reaction_idx = self.rng.choice(len(propensities), p=propensities/a0)
        
        # Apply stoichiometry
        S = self.rn.stoichiometry_matrix()
        x_new = x + S[:, reaction_idx]
        
        return x_new, t + tau, reaction_idx


class TauLeaping:
    """
    Tau-leaping algorithm (approximate stochastic simulation).
    
    Faster than Gillespie SSA by allowing multiple reactions per time step.
    Good compromise between accuracy and speed.
    
    Reference: Gillespie, D. T. (2001). J. Chem. Phys. 115(4), 1716–1733.
    """
    pass  # Implementation follows similar structure to GillespieSSA


class LangevinDynamics:
    """
    Stochastic Langevin dynamics (continuous noise).
    
    Adds white noise terms to deterministic ODEs:
    dx/dt = f(x) + σ(x) * ξ(t)
    
    where ξ(t) is Gaussian white noise.
    
    Good for systems with many molecules where noise is continuous.
    """
    pass  # Requires specialized integrators (Euler-Maruyama, Milstein)


class BudgetAllocationSimulator:
    """
    Discrete-time stochastic simulator using species budget allocation.

    Conceptual model
    ----------------
    At each time step dt every species has a budget (its current amount).
    Each unit of a species can commit to exactly one reaction that consumes it,
    or stay idle.  The allocation is stochastic: for consumed species X and
    reaction r the per-unit commitment probability is

        p(X -> r) = stoich_X * lambda_r * dt / [X]

    where lambda_r is the reaction propensity (events / time unit).  Because
    individual probabilities must sum to <= 1 the remainder is the idle
    fraction -- no free parameter required.  Species that appear only as
    catalysts (net_change == 0) contribute to the propensity formula but draw
    no budget.

    Propensity formula
    ------------------
    lambda_r = k_r
               * seasonal_factor(t)          (if reaction is seasonal)
               * product over consumed  species:  g(X_i)
               * product over enabling  species:  g(X_i)

    where  g(X) = X / (X + K_X)  if species X has a K entry in kinetics,
                = X               otherwise (linear).
    Enabling species appearing multiple times in the list contribute that
    many multiplicative factors (e.g. twice for a squared catalytic effect).

    reaction_defs format
    --------------------
    {
      'r7': {
        'consumed' : {'SR_RL': 1},        # draws from budget; enter propensity
        'enabling' : ['Gov'],              # catalyst; enters propensity only
        'net_change': {'SR_RL': -1, 'SR_SL': +1, 'E_peace': +1},
        'seasonal'  : None,               # None | 'recovery' | 'degradation'
      }, ...
    }

    kinetics format
    ---------------
    {
      'r7': {'k': 0.1, 'K': {'SR_RL': 20.0}},
      ...
    }
    K entries select the saturating form for that species.  Absent = linear.

    Parameters
    ----------
    reaction_defs : dict
    kinetics      : dict
    dt            : float   Time step in years  (default 0.05)
    climate       : float   [0,1]  0=dry, 0.5=neutral, 1=wet.
                            Scales seasonal recovery/degradation amplitudes.
    rng_seed      : int, optional
    """

    def __init__(self, reaction_defs, kinetics, dt=0.05, climate=0.5,
                 rng_seed=None):
        self.reaction_defs = reaction_defs
        self.kinetics      = kinetics
        self.dt            = dt
        self.climate       = float(np.clip(climate, 0.0, 1.0))
        self.rng           = np.random.default_rng(rng_seed)
        self._rxn_ids      = list(reaction_defs.keys())

    # ── internal helpers ────────────────────────────────────────────────────

    def _seasonal_factor(self, t: float, phase: str) -> float:
        """
        Seasonal modulation in [0, 1], scaled by climate.

        phase='recovery'    peaks at t=0,1,2... (wet season).
            amplitude = climate  (wetter -> stronger recovery)
        phase='degradation' peaks at t=0.5,1.5... (dry season).
            amplitude = 1-climate (drier -> stronger degradation)
        The factor is multiplied by 2 so the mean over the year equals the
        amplitude, keeping the overall rate comparable across climate values.
        """
        base = (1.0 + np.cos(2.0 * np.pi * t)) / 2.0  # in [0, 1]
        if phase == 'degradation':
            base = 1.0 - base
        amplitude = self.climate if phase == 'recovery' else (1.0/self.climate)
        return base * amplitude * 2.0

    # ── propensity ──────────────────────────────────────────────────────────

    def compute_propensity(self, rxn_id: str, state: dict, t: float) -> float:
        """
        Compute lambda_r (events per unit time) for one reaction at time t.
        """
        rxn = self.reaction_defs[rxn_id]
        kin = self.kinetics[rxn_id]
        K   = kin.get('K', {})
        lam = float(kin['k'])

        # Seasonal modulation
        seasonal = rxn.get('seasonal')
        if seasonal in ('recovery', 'degradation'):
            lam *= self._seasonal_factor(t, seasonal)

        # Consumed and enabling species both enter the propensity
        def _g(sp):
            x = max(0.0, state.get(sp, 0.0))
            return x / (x + K[sp]) if sp in K else x

        for sp in rxn.get('consumed', {}):
            lam *= _g(sp)
        for sp in rxn.get('enabling', []):
            lam *= _g(sp)

        return lam

    def compute_all_propensities(self, state: dict, t: float) -> dict:
        """Return {rxn_id: lambda_r} for all reactions."""
        return {r: self.compute_propensity(r, state, t) for r in self._rxn_ids}

    # ── allocation ──────────────────────────────────────────────────────────

    def _compute_alloc_probs(self, propensities: dict, state: dict) -> dict:
        """
        For each consumed species compute per-reaction allocation probability.

        p(X_i -> r) = stoich_i * lambda_r * dt / max([X_i], eps)

        If the sum over all reactions for a species exceeds 1 the probabilities
        are renormalised and a RuntimeWarning is issued (dt too large).
        Returns {species: {rxn_id: prob}}.
        """
        import warnings
        alloc: Dict[str, Dict[str, float]] = {}

        for rxn_id, lam in propensities.items():
            consumed = self.reaction_defs[rxn_id].get('consumed', {})
            for sp, stoich in consumed.items():
                x = max(state.get(sp, 0.0), 1e-12)
                p = stoich * lam * self.dt / x
                bucket = alloc.setdefault(sp, {})
                bucket[rxn_id] = bucket.get(rxn_id, 0.0) + p

        for sp, rxn_probs in alloc.items():
            total = sum(rxn_probs.values())
            if total > 1.0:
                warnings.warn(
                    f"BudgetAllocationSimulator: species '{sp}' allocation "
                    f"sum = {total:.3f} > 1.  Reduce dt or k values. "
                    f"Clipping to 1.",
                    RuntimeWarning, stacklevel=3,
                )
                scale = 1.0 / total
                for r in rxn_probs:
                    rxn_probs[r] *= scale

        return alloc

    def _stochastic_draw(self, alloc_probs: dict, state: dict) -> dict:
        """
        For each (species, reaction) pair draw committed units from
        Poisson(state[X] * p(X->r)).  Works for continuous amounts.
        Returns {rxn_id: {species: committed_float}}.
        """
        committed: Dict[str, Dict[str, float]] = {r: {} for r in self._rxn_ids}
        for sp, rxn_probs in alloc_probs.items():
            x = max(state.get(sp, 0.0), 0.0)
            for rxn_id, p in rxn_probs.items():
                expected = x * p
                drawn = float(self.rng.poisson(max(expected, 0.0)))
                committed[rxn_id][sp] = drawn
        return committed

    def _resolve_firings(self, committed: dict, propensities: dict) -> dict:
        """
        For reactions with consumed species:
            firings = floor( min over consumed of committed[sp] / stoich )
        For pure catalytic reactions (no consumed species):
            firings ~ Poisson( lambda_r * dt )
        Returns {rxn_id: n_firings}.
        """
        firings: Dict[str, int] = {}
        for rxn_id in self._rxn_ids:
            consumed = self.reaction_defs[rxn_id].get('consumed', {})
            if consumed:
                slots = [
                    committed[rxn_id].get(sp, 0.0) / stoich
                    for sp, stoich in consumed.items()
                ]
                firings[rxn_id] = int(min(slots)) if slots else 0
            else:
                # Pure catalytic: Poisson draw from expected rate
                lam = propensities[rxn_id]
                firings[rxn_id] = int(self.rng.poisson(lam * self.dt))
        return firings

    def _apply_firings(self, state: dict, firings: dict) -> dict:
        """
        Apply net_change * n_firings to state.  All values clipped to >= 0.
        Returns new state dict.
        """
        new_state = dict(state)
        for rxn_id, n in firings.items():
            if n == 0:
                continue
            for sp, delta in self.reaction_defs[rxn_id].get('net_change', {}).items():
                new_state[sp] = max(0.0, new_state.get(sp, 0.0) + delta * n)
        return new_state

    # ── public simulation API ───────────────────────────────────────────────

    def step(self, state: dict, t: float) -> Tuple[dict, dict]:
        """
        Advance simulation by one time step dt from time t.

        Returns
        -------
        new_state : dict   {species: amount}
        flux      : dict   {rxn_id: firings / dt}
        """
        propensities = self.compute_all_propensities(state, t)
        alloc_probs  = self._compute_alloc_probs(propensities, state)
        committed    = self._stochastic_draw(alloc_probs, state)
        firings      = self._resolve_firings(committed, propensities)
        new_state    = self._apply_firings(state, firings)
        flux         = {r: firings[r] / self.dt for r in self._rxn_ids}
        return new_state, flux

    def simulate(self, ic: dict, T_max: float):
        """
        Run simulation from initial condition ic for T_max years.

        Parameters
        ----------
        ic    : dict  {species: initial_amount}
        T_max : float  simulation horizon in years

        Returns
        -------
        trajectory_df : pd.DataFrame  index=Time, columns=species
        flux_df       : pd.DataFrame  index=Time, columns=rxn_ids
                        (time points t=dt, 2dt, ..., T_max)
        """
        import pandas as pd
        state   = {sp: float(v) for sp, v in ic.items()}
        n_steps = int(round(T_max / self.dt))

        states  = [dict(state)]
        fluxes  = []
        times   = [0.0]

        for i in range(n_steps):
            t = i * self.dt
            state, flux = self.step(state, t)
            states.append(dict(state))
            fluxes.append(dict(flux))
            times.append(t + self.dt)

        traj_df = pd.DataFrame(states, index=times)
        traj_df.index.name = 'Time'

        flux_df = pd.DataFrame(fluxes, index=times[1:])
        flux_df.index.name = 'Time'

        return traj_df, flux_df


# Registry of stochastic methods
STOCHASTIC_METHODS = {
    'poisson_events': PoissonEventGenerator,
    'gillespie': GillespieSSA,
    'tau_leaping': TauLeaping,
    'langevin': LangevinDynamics
}


__all__ = [
    'PoissonEventGenerator',
    'concentration_shock',
    'rate_perturbation_event',
    'GillespieSSA',
    'TauLeaping',
    'LangevinDynamics',
    'BudgetAllocationSimulator',
    'STOCHASTIC_METHODS',
]
