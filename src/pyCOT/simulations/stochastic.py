from __future__ import annotations

import math
import warnings

from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from pyCOT.kinetics import KINETIC_REGISTRY
from pyCOT.simulations.core import (
    build_reaction_dict,
    parse_parameters,
    validate_rate_list,
    generate_default_parameters,
    print_differential_equations,
    print_velocity_expressions,
    time_series_dataframe,
    flux_vector_dataframe,
    generate_random_vector,
)


def _propensity_mak_exact(
    state: np.ndarray,
    reactant_coeffs: np.ndarray,
    k: float,
) -> float:
    """
    Exact stochastic mass-action propensity:

        a_j(x) = k * prod_i C(x_i, n_ij)

    where n_ij is the stoichiometric coefficient of species i
    as a reactant in reaction j.
    """
    a = float(k)

    for x_i, n_i in zip(state, reactant_coeffs):

        if n_i == 0:
            continue

        x_i = int(round(x_i))
        n_i = int(round(n_i))

        if x_i < n_i:
            return 0.0

        a *= math.comb(x_i, n_i)

    return max(a, 0.0)


def _prepare_transport_rates(
    species: List[str],
    num_patches: int,
    D_dict: Optional[Dict[str, float]],
    connectivity_matrix: Optional[np.ndarray],
    transport_rates: Optional[Dict[str, np.ndarray]],
) -> Dict[str, np.ndarray]:
    """
    Construct species-specific transport-rate matrices.

    Convention:

        T[s][p, q] = rate at which one molecule of species s
                     moves from patch p to patch q.

    Diagonal entries must be zero.
    """

    if transport_rates is not None and (
        D_dict is not None or connectivity_matrix is not None
    ):
        raise ValueError(
            "Use either transport_rates or "
            "(D_dict + connectivity_matrix), not both."
        )

    T = {}

    # ------------------------------------------------------------
    # General formulation
    # ------------------------------------------------------------
    if transport_rates is not None:

        for s in species:

            if s not in transport_rates:
                T[s] = np.zeros((num_patches, num_patches))
                continue

            M = np.asarray(
                transport_rates[s],
                dtype=float,
            )

            if M.shape != (num_patches, num_patches):
                raise ValueError(
                    f"transport_rates['{s}'] must have shape "
                    f"({num_patches}, {num_patches}), "
                    f"received {M.shape}."
                )

            if np.any(M < 0):
                raise ValueError(
                    f"transport_rates['{s}'] cannot contain negative rates."
                )

            M = M.copy()
            np.fill_diagonal(M, 0.0)

            T[s] = M

        return T

    # ------------------------------------------------------------
    # Special case T^(s) = D_s C
    # ------------------------------------------------------------

    if connectivity_matrix is None:

        connectivity_matrix = np.zeros(
            (num_patches, num_patches),
            dtype=float,
        )

    C = np.asarray(
        connectivity_matrix,
        dtype=float,
    )

    if C.shape != (num_patches, num_patches):
        raise ValueError(
            f"connectivity_matrix must have shape "
            f"({num_patches}, {num_patches}), "
            f"received {C.shape}."
        )

    if np.any(C < 0):
        raise ValueError(
            "connectivity_matrix cannot contain negative values."
        )

    C = C.copy()
    np.fill_diagonal(C, 0.0)

    D_dict = D_dict or {}

    for s in species:

        D = float(D_dict.get(s, 0.0))

        if D < 0:
            raise ValueError(
                f"D_dict['{s}'] cannot be negative."
            )

        T[s] = D * C

    return T


def gillespie(
    rn,
    rate: Union[str, List[str]] = 'mak',
    spec_vector: Optional[List] = None,
    x0: Optional[Dict[str, List]] = None,
    t_span: Tuple[float, float] = (0, 100),
    n_steps: Optional[int] = 101,
    additional_laws: Optional[Dict] = None,
    max_iter: int = 100_000,
    exact_mass_action: bool = True,
    seed: Optional[int] = None,
    stop_condition: Optional[Callable] = None,
    verbose: bool = True,

    # ============================================================
    # Metapopulation
    # ============================================================
    num_patches: int = 1,
    D_dict: Optional[Dict[str, float]] = None,
    connectivity_matrix: Optional[np.ndarray] = None,
    transport_rates: Optional[Dict[str, np.ndarray]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Gillespie Direct Method for a reaction-network metapopulation.

    The stochastic system contains TWO types of events:

        1. Local reaction events
        2. Inter-patch dispersal events

    Reaction event:

        (j, p):
            X^(p) -> X^(p) + nu_j

    with propensity

        a_R[j,p] = k_j prod_i C(X_i^(p), n_ij)

    for exact mass-action kinetics.

    Dispersal event:

        (s, p -> q):
            X_s^(p) -= 1
            X_s^(q) += 1

    with propensity

        a_D[s,p->q] = T_s[p,q] * X_s^(p)

    Parameters
    ----------
    rn : ReactionNetwork
        pyCOT reaction network.

    rate : str or list[str]
        Kinetic law(s).

    spec_vector : list
        Reaction parameters.

    x0 : dict
        Initial populations.

        Single patch:

            {
                'l':  [10],
                's1': [15],
                's2': [5]
            }

        Multiple patches:

            {
                'l':  [100, 20, 0],
                's1': [10, 50, 20],
                's2': [0, 20, 80]
            }

    t_span : tuple
        Simulation interval.

    n_steps : int, optional
        Number of points in the output grid.

    additional_laws : dict, optional
        Additional kinetic laws.

    max_iter : int
        Maximum number of Gillespie events.

    exact_mass_action : bool
        If True, exact stochastic mass-action propensities are used.

    seed : int, optional
        Random seed.

    stop_condition : callable, optional
        Function receiving (X, species) that can stop the simulation.

    verbose : bool
        If True, prints model information and simulation summary.
        Velocity expressions and differential equations are printed
        through the standard core helpers, so the SSA shows exactly
        the same model description as the ODE simulator.

    num_patches : int
        Number of patches.

    D_dict : dict, optional
        Species-specific dispersal coefficients.

    connectivity_matrix : ndarray, optional
        Connectivity matrix.

    transport_rates : dict, optional
        Explicit species-specific transport matrices.

    Returns
    -------
    time_series_df : pandas.DataFrame

        Population trajectories.

    flux_vector_df : pandas.DataFrame

        Reaction propensities evaluated at output states.
    """

    # ============================================================
    # 0. Random generator
    # ============================================================

    rng = np.random.default_rng(seed)

    # ============================================================
    # 1. Reaction network
    # ============================================================

    sm = rn.stoichiometry_matrix()

    mat_species = list(sm.species)
    mat_reactions = list(sm.reactions)

    out_species = [
        s.name for s in rn.species()
    ]

    out_reactions = [
        r.name() for r in rn.reactions()
    ]

    if sorted(mat_species) != sorted(out_species):

        raise ValueError(
            "Mismatch between species() and "
            "stoichiometry_matrix().species"
        )

    if sorted(mat_reactions) != sorted(out_reactions):

        raise ValueError(
            "Mismatch between reactions() and "
            "stoichiometry_matrix().reactions"
        )

    n_species = len(mat_species)
    n_reactions = len(mat_reactions)

    reactant_matrix = np.asarray(
        rn.reactants_matrix(),
        dtype=float,
    )

    stoich_matrix = np.asarray(
        sm,
        dtype=float,
    )

    rn_dict = build_reaction_dict(rn)

    # ============================================================
    # 2. Kinetic laws
    # ============================================================

    rate = validate_rate_list(
        rate,
        n_reactions,
    )

    rate_laws = dict(KINETIC_REGISTRY)

    if additional_laws:

        for name, fn in additional_laws.items():

            rate_laws[name] = fn

    # ============================================================
    # 3. Exactness
    # ============================================================

    if exact_mass_action:

        non_mak = [
            r
            for r in rate
            if r != 'mak'
        ]

        if non_mak:

            raise ValueError(
                "Exact metapopulation SSA requires 'mak'. "
                f"Received kinetic laws: {non_mak}. "
                "Use exact_mass_action=False only if you "
                "explicitly want phenomenological propensities."
            )

    # ============================================================
    # 4. Parameters
    # ============================================================

    if spec_vector is None:

        spec_vector = generate_default_parameters(
            rate,
            n_reactions,
            additional_laws,
        )

    if len(spec_vector) != n_reactions:

        raise ValueError(
            f"spec_vector must have {n_reactions} entries, "
            f"received {len(spec_vector)}."
        )

    parameters = parse_parameters(
        rn,
        rn_dict,
        rate,
        spec_vector,
    )

    # ============================================================
    # 5. Initial state
    # ============================================================

    X = np.zeros(
        (n_species, num_patches),
        dtype=float,
    )

    if x0 is None:

        for i, s in enumerate(mat_species):

            values = generate_random_vector(
                num_patches,
                seed=None,
                min_value=50,
                max_value=100,
            ).astype(int)

            X[i, :] = values

    else:

        for s, values in x0.items():

            if s not in mat_species:

                raise ValueError(
                    f"Species '{s}' does not exist "
                    "in the network."
                )

            values = np.asarray(
                values,
                dtype=float,
            )

            if values.shape != (num_patches,):

                raise ValueError(
                    f"x0['{s}'] must have "
                    f"{num_patches} values."
                )

            if np.any(values < 0):

                raise ValueError(
                    f"x0['{s}'] cannot contain "
                    "negative populations."
                )

            X[
                mat_species.index(s),
                :
            ] = np.round(values)

    # ============================================================
    # 6. Transport matrices
    # ============================================================

    T = _prepare_transport_rates(
        species=mat_species,
        num_patches=num_patches,
        D_dict=D_dict,
        connectivity_matrix=connectivity_matrix,
        transport_rates=transport_rates,
    )

    # ============================================================
    # 7. Species indices
    # ============================================================

    species_idx = {
        s: i
        for i, s in enumerate(mat_species)
    }

    # ============================================================
    # 8. Verbose model information
    # ============================================================
    #
    # The description of the model uses the same core helpers as the
    # ODE simulator, so that both backends report the reactions,
    # velocity expressions and differential equations identically.
    #
    #   print_velocity_expressions(rn, rn_dict, rate, rate_laws,
    #                              additional_laws)
    #   print_differential_equations(rn, rn_dict)
    #
    # Only the SSA-specific information (x0 per patch, spec_vector,
    # parsed parameters, patches, seed) is printed here.

    if verbose:

        print("\n" + "=" * 70)
        print("GILLESPIE DIRECT METHOD")
        print("=" * 70)

        print("\n[1] Reaction network")
        print("-" * 70)

        print(f"Species   : {mat_species}")
        print(f"Reactions : {mat_reactions}")
        print(f"Rate laws : {rate}")
        print(f"Patches   : {num_patches}") 
        print("Initial condition x0 =", x0)
        print("spec_vector =", spec_vector)
        print("-" * 70)

        # # --------------------------------------------------------
        # # Parsed parameters
        # # --------------------------------------------------------

        # print("\n[4] Parsed parameters")
        # print("-" * 70)

        # for reaction in mat_reactions:

        #     print(f"  {reaction}: {parameters[reaction]}")

        # --------------------------------------------------------
        # Velocity expressions (core helper)
        # --------------------------------------------------------

        print("\n[2] Velocity expressions")
        print("-" * 70)

        print_velocity_expressions(
            rn,
            rn_dict,
            rate,
            rate_laws,
            additional_laws,
        )

        # --------------------------------------------------------
        # Differential equations (core helper)
        # --------------------------------------------------------

        print("\n[3] Differential equations")
        print("-" * 70)

        print_differential_equations(
            rn,
            rn_dict,
        )

        print("\n" + "=" * 70)

    # ============================================================
    # 9. Reaction propensities
    # ============================================================

    def reaction_propensities(
        X,
        t,
    ):

        a = np.zeros(
            (n_reactions, num_patches),
            dtype=float,
        )

        for p in range(num_patches):

            x_local = X[:, p]

            for j, reaction in enumerate(
                mat_reactions
            ):

                reactants = rn_dict[
                    reaction
                ][0]

                # =================================================
                # Exact stochastic mass action
                # =================================================

                if (
                    rate[j] == 'mak'
                    and exact_mass_action
                ):

                    a[j, p] = (
                        _propensity_mak_exact(
                            x_local,
                            reactant_matrix[:, j],
                            float(
                                parameters[
                                    reaction
                                ][0]
                            ),
                        )
                    )

                # =================================================
                # Phenomenological propensity
                # =================================================

                else:

                    if np.any(
                        x_local
                        < reactant_matrix[:, j]
                    ):

                        continue

                    fn = rate_laws[
                        rate[j]
                    ]

                    a[j, p] = max(
                        float(
                            fn(
                                reactants,
                                x_local,
                                species_idx,
                                parameters[
                                    reaction
                                ],
                            )
                        ),
                        0.0,
                    )

        return a

    # ============================================================
    # 10. Dispersal propensities
    # ============================================================

    def dispersal_events(X):

        events = []

        propensities = []

        for i, s in enumerate(
            mat_species
        ):

            T_s = T[s]

            for p in range(
                num_patches
            ):

                N = X[i, p]

                if N <= 0:
                    continue

                for q in range(
                    num_patches
                ):

                    if p == q:
                        continue

                    rate_pq = T_s[p, q]

                    if rate_pq <= 0:
                        continue

                    a = rate_pq * N

                    events.append(
                        (i, p, q)
                    )

                    propensities.append(
                        a
                    )

        return (
            events,
            np.asarray(
                propensities,
                dtype=float,
            ),
        )

    # ============================================================
    # 11. Global SSA loop
    # ============================================================

    t0, t_end = map(
        float,
        t_span,
    )

    t = t0

    times = [t]

    states = [
        X.copy()
    ]

    reaction_event_counts = {

        r: np.zeros(
            num_patches,
            dtype=int,
        )

        for r in mat_reactions
    }

    dispersal_event_counts = {

        s: np.zeros(
            (
                num_patches,
                num_patches,
            ),
            dtype=int,
        )

        for s in mat_species
    }

    n_events_executed = 0

    # ============================================================
    # SSA
    # ============================================================

    for iteration in range(
        max_iter
    ):

        # --------------------------------------------------------
        # Reaction propensities
        # --------------------------------------------------------

        A_R = reaction_propensities(
            X,
            t,
        )

        reaction_events = []

        reaction_prop = []

        for j in range(
            n_reactions
        ):

            for p in range(
                num_patches
            ):

                if A_R[j, p] > 0:

                    reaction_events.append(
                        (j, p)
                    )

                    reaction_prop.append(
                        A_R[j, p]
                    )

        reaction_prop = np.asarray(
            reaction_prop,
            dtype=float,
        )

        # --------------------------------------------------------
        # Dispersal propensities
        # --------------------------------------------------------

        (
            dispersion_events,
            dispersion_prop,
        ) = dispersal_events(X)

        # --------------------------------------------------------
        # Combine all events
        # --------------------------------------------------------

        all_events = (

            [
                ('reaction', e)
                for e in reaction_events
            ]

            +

            [
                ('dispersal', e)
                for e in dispersion_events
            ]
        )

        all_prop = np.concatenate(
            [
                reaction_prop,
                dispersion_prop,
            ]
        )

        a0 = all_prop.sum()

        # --------------------------------------------------------
        # Absorbing state
        # --------------------------------------------------------

        if a0 <= 0:
            break

        # --------------------------------------------------------
        # Gillespie waiting time
        # --------------------------------------------------------

        u1 = 1.0 - rng.random()

        tau = (
            -math.log(u1)
            / a0
        )

        if t + tau > t_end:
            break

        # --------------------------------------------------------
        # Select event
        # --------------------------------------------------------

        u2 = rng.random()

        threshold = (
            u2 * a0
        )

        cumulative = np.cumsum(
            all_prop
        )

        event_idx = int(
            np.searchsorted(
                cumulative,
                threshold,
                side='right',
            )
        )

        event_idx = min(
            event_idx,
            len(all_events) - 1,
        )

        event_type, event = (
            all_events[event_idx]
        )

        # --------------------------------------------------------
        # Execute reaction event
        # --------------------------------------------------------

        if event_type == 'reaction':

            j, p = event

            X[:, p] += (
                stoich_matrix[:, j]
            )

            reaction_event_counts[
                mat_reactions[j]
            ][p] += 1

        # --------------------------------------------------------
        # Execute dispersal event
        # --------------------------------------------------------

        else:

            i, p, q = event

            X[i, p] -= 1.0
            X[i, q] += 1.0

            dispersal_event_counts[
                mat_species[i]
            ][p, q] += 1

        # --------------------------------------------------------
        # Numerical consistency check
        # --------------------------------------------------------

        if np.any(
            X < -1e-12
        ):

            raise RuntimeError(
                "Negative population generated. "
                "This indicates an inconsistency "
                "in the propensity/event definition."
            )

        X[
            np.abs(X) < 1e-12
        ] = 0.0

        # --------------------------------------------------------
        # Advance time
        # --------------------------------------------------------

        t += tau

        times.append(t)

        states.append(
            X.copy()
        )

        n_events_executed += 1

        # --------------------------------------------------------
        # Optional stopping condition
        # --------------------------------------------------------

        if stop_condition is not None:

            if stop_condition(
                X,
                mat_species,
            ):

                break

    # ============================================================
    # 12. Convert trajectory
    # ============================================================

    times = np.asarray(
        times,
        dtype=float,
    )

    states = np.asarray(
        states,
        dtype=float,
    )

    # ============================================================
    # 13. Output grid
    # ============================================================

    if n_steps is not None:

        grid = np.linspace(
            t0,
            t_end,
            n_steps,
        )

        idx = (
            np.searchsorted(
                times,
                grid,
                side='right',
            )
            - 1
        )

        idx = np.clip(
            idx,
            0,
            len(times) - 1,
        )

        states = states[idx]

        times = grid

    # ============================================================
    # 14. Population DataFrame
    # ============================================================

    data = {
        'Time': times
    }

    if num_patches == 1:

        for i, s in enumerate(
            mat_species
        ):

            data[s] = (
                states[:, i, 0]
            )

    else:

        for i, s in enumerate(
            mat_species
        ):

            for p in range(
                num_patches
            ):

                data[
                    f'{s}_{p + 1}'
                ] = (
                    states[:, i, p]
                )

    time_series_df = pd.DataFrame(
        data
    )

    # ============================================================
    # 15. Reaction propensities at output states
    # ============================================================

    flux_data = {
        'Time': times
    }

    propensity_states = [

        reaction_propensities(
            state,
            0.0,
        )

        for state in states
    ]

    if num_patches == 1:

        for j, r in enumerate(
            mat_reactions
        ):

            flux_data[r] = [

                A[j, 0]

                for A in propensity_states
            ]

    else:

        for j, r in enumerate(
            mat_reactions
        ):

            for p in range(
                num_patches
            ):

                flux_data[
                    f'{r}_{p + 1}'
                ] = [

                    A[j, p]

                    for A in propensity_states
                ]

    flux_vector_df = pd.DataFrame(
        flux_data
    )

    # ============================================================
    # 16. Verbose simulation summary
    # ============================================================

    if verbose:

        print(
            "\n[4] Simulation summary"
        )

        print("-" * 70)

        print(
            f"Time interval : "
            f"[{t0:g}, {t_end:g}]"
        )

        print(
            f"Output points : "
            f"{len(times)}"
        )

        print(
            f"SSA events    : "
            f"{n_events_executed}"
        )

        print(
            f"Maximum events: "
            f"{max_iter}"
        )

        print(
            f"Seed          : "
            f"{seed}"
        )

        # --------------------------------------------------------
        # Reaction events
        # --------------------------------------------------------

        print(
            "\nReaction events:"
        )

        for r, counts in (
            reaction_event_counts.items()
        ):

            print(
                f"  {r}: "
                f"{counts.tolist()}"
            )

        # --------------------------------------------------------
        # Dispersal events
        # --------------------------------------------------------

        print(
            "\nDispersal events:"
        )

        any_dispersal = False

        for s, counts in (
            dispersal_event_counts.items()
        ):

            total = counts.sum()

            if total > 0:

                any_dispersal = True

                print(
                    f"  {s}: "
                    f"{total}"
                )

                if num_patches > 1:

                    for p in range(
                        num_patches
                    ):

                        for q in range(
                            num_patches
                        ):

                            if p == q:
                                continue

                            n = counts[p, q]

                            if n > 0:

                                print(
                                    f"      "
                                    f"{p + 1} -> "
                                    f"{q + 1}: "
                                    f"{n}"
                                )

        if not any_dispersal:

            print(
                "  None"
            )

        print(
            "\n" + "=" * 70
        )

    # ============================================================
    # 17. Return
    # ============================================================

    return (
        time_series_df,
        flux_vector_df,
    )


__all__ = ['gillespie']





################ Gillespie v.0 #######################
# """
# Stochastic simulation (Gillespie SSA) for pyCOT reaction networks.

# Mirrors the interface of `pyCOT.simulations.ode.simulation`: it accepts `rate`,
# `spec_vector`, `x0` and `additional_laws`, generates random initial conditions
# and parameters when they are not supplied, and returns the same two DataFrames
# (time series and flux vector) using the utilities in
# `pyCOT.simulations.core`.

# Essential difference with the ODE version: here time advances through discrete
# events (Gillespie's direct method), not through numerical integration.

# Warnings about kinetic laws
# ---------------------------
# The SSA is exact only under mass action ('mak'), where the propensity is
#     a_j = k_j * prod_i C(x_i, n_ij)
# with C the binomial coefficient (Gillespie 1977). The remaining laws in the
# registry ('mmk', 'hill', 'saturated') are phenomenological macroscopic rate
# laws: using them as propensities is a heuristic approximation with no grounding
# in the chemical master equation. The module accepts them for compatibility, but
# issues a warning.

# Time-dependent laws ('cosine', 'threshold_memory') violate the assumption of
# constant propensity between events on which tau = -ln(u)/a0 is derived. With
# them the simulation is NOT exact and strictly requires a thinning scheme or the
# non-Markovian Gillespie algorithm. An explicit warning is issued.

# Requires: numpy, pandas, scipy (via core). Optional: matplotlib.
# """

# from __future__ import annotations

# import math
# import warnings
# from typing import Callable, Dict, List, Optional, Tuple, Union

# import numpy as np
# import pandas as pd

# from pyCOT.kinetics import KINETIC_REGISTRY
# from pyCOT.simulations.core import (
#     validate_rate_list,
#     build_reaction_dict,
#     parse_parameters,
#     generate_default_parameters,
#     print_differential_equations,
#     print_velocity_expressions,
#     time_series_dataframe,
#     flux_vector_dataframe,
#     generate_random_vector,
# )

# # Laws whose propensity depends explicitly on t
# _TIME_DEPENDENT = ('cosine', 'threshold_memory')

# # The only law with an exact grounding in the chemical master equation
# _EXACT_LAWS = ('mak',)


# # ---------------------------------------------------------------------------
# # Reproducible generation of default parameters
# # ---------------------------------------------------------------------------

# def generate_default_parameters_seeded(rate_list: List[str],
#                                        n_reactions: int,
#                                        additional_laws: Optional[Dict] = None,
#                                        seed: Optional[int] = None) -> List[List[float]]:
#     """Seeded version of `core.generate_default_parameters`.

#     `core.generate_default_parameters` does not propagate a seed (it calls
#     `generate_random_vector` without `seed`), so two runs with the same
#     Gillespie `seed` would draw different parameters. This function replicates
#     exactly the same sampling ranges as core, but with a single seeded
#     generator, so that the whole simulation is reproducible.

#     If `seed` is None, it delegates to core's implementation.
#     """
#     if seed is None:
#         return generate_default_parameters(rate_list, n_reactions, additional_laws)

#     rng = np.random.default_rng(seed)

#     def u(lo, hi):
#         return float(np.round(rng.uniform(lo, hi), 2))

#     spec_vector = []
#     for kinetic in rate_list:
#         if kinetic == 'mak':
#             params = [u(0.01, 1.0)]
#         elif kinetic == 'mmk':
#             params = [u(1, 1.5), u(5, 10)]
#         elif kinetic == 'hill':
#             params = [u(1, 1.5), u(5, 10), u(1, 4)]
#         elif kinetic == 'cosine':
#             params = [u(0.5, 2.0), u(0.1, 1.0)]
#         elif kinetic == 'saturated':
#             params = [u(1, 1.5), u(5, 10)]
#         elif kinetic == 'threshold_memory':
#             params = [u(0.5, 2.0), u(0.5, 1.5), u(0.1, 0.8), u(0.1, 0.5), 5.0]
#         elif kinetic in (additional_laws or {}):
#             params = [float(np.round(v, 3)) for v in rng.uniform(0.1, 1.0, 3)]
#         else:
#             raise ValueError(f"Unknown kinetic law: {kinetic}")
#         spec_vector.append(params)

#     return spec_vector


# # ---------------------------------------------------------------------------
# # Propensities
# # ---------------------------------------------------------------------------

# def _propensity_mak_exact(state: np.ndarray, reactant_coeffs: np.ndarray, k: float) -> float:
#     """Exact stochastic propensity: k * prod_i C(x_i, n_ij)."""
#     a = k
#     for x_i, n_i in zip(state, reactant_coeffs):
#         if n_i == 0:
#             continue
#         x_i_int = int(round(x_i))
#         n_i_int = int(round(n_i))
#         if x_i_int < n_i_int:
#             return 0.0
#         a *= math.comb(x_i_int, n_i_int)
#     return a


# # ---------------------------------------------------------------------------
# # Simulator
# # ---------------------------------------------------------------------------

# def gillespie(rn,
#               rate: Union[str, List[str]] = 'mak',
#               spec_vector: Optional[List] = None,
#               x0: Optional[List] = None,
#               t_span: Tuple[float, float] = (0, 10),
#               n_steps: Optional[int] = None,
#               additional_laws: Optional[Dict] = None,
#               max_iter: int = 100_000,
#               exact_mass_action: bool = True,
#               seed: Optional[int] = None,
#               stop_condition: Optional[Callable] = None,
#               verbose: bool = True) -> Tuple[pd.DataFrame, pd.DataFrame]:
#     """
#     Simulate a pyCOT reaction network with the Gillespie algorithm (SSA).

#     Interface parallel to `pyCOT.simulations.ode.simulation`: if `x0` or
#     `spec_vector` are not supplied, they are generated at random following the
#     same conventions as `core` (integer x0 in [50, 100]; parameters according
#     to each reaction's kinetic law).

#     Parameters
#     ----------
#     rn : ReactionNetwork
#         pyCOT reaction network object.
#     rate : str or list
#         Kinetic law name(s), a single one or one per reaction.
#         Available: 'mak', 'mmk', 'hill', 'cosine', 'saturated',
#         'threshold_memory'. Only 'mak' yields an exact SSA (see the module
#         docstring).
#     spec_vector : list, optional
#         Parameters per reaction, format [[params_r1], [params_r2], ...].
#         If None, defaults are generated. For 'mak' the only parameter is the
#         rate constant k_j.
#     x0 : list, optional
#         Initial populations (integer counts). If None, random integer values
#         in [50, 100] are generated. A dict {species_name: count} is also
#         accepted.
#     t_span : tuple
#         Time interval (t_start, t_end).
#     n_steps : int, optional
#         If given, the output is resampled onto a uniform grid of `n_steps`
#         points using zero-order hold (the state is constant between events,
#         which is exactly the SSA trajectory). Useful for averaging replicates
#         or comparing against the ODE. If None, the original event times are
#         returned.
#     additional_laws : dict, optional
#         Custom kinetic laws {name: function}.
#     max_iter : int
#         Maximum number of events, to bound the run.
#     exact_mass_action : bool
#         If True, reactions with the 'mak' law use the exact combinatorial
#         propensity k*prod C(x_i, n_ij) instead of the macroscopic form
#         k*prod x_i^n_ij from the registry. Recommended: both coincide when all
#         coefficients are 1, but only the combinatorial form is correct with
#         larger coefficients.
#     seed : int, optional
#         Random seed. Controls the random x0, the random spec_vector and the
#         stochastic trajectory, so that the run is reproducible.
#     stop_condition : callable, optional
#         f(state, species_names) -> bool. Stops the simulation if True.
#     verbose : bool
#         Print equations, velocity expressions and parameters.

#     Returns
#     -------
#     time_series_df : DataFrame
#         Columns: Time, species_1, species_2, ...
#     flux_vector_df : DataFrame
#         Columns: Time, reaction_1, reaction_2, ...

#     Examples
#     --------
#     >>> ts, flux = gillespie(rn)                       # everything random
#     >>> ts, flux = gillespie(rn, seed=42)              # random but reproducible
#     >>> ts, flux = gillespie(rn, spec_vector=[[0.7], [0.5], [1.0]],
#     ...                      x0={'l': 10, 's1': 5, 's2': 0})
#     """
#     rng = np.random.default_rng(seed)

#     # --- Canonical ordering ----------------------------------------------
#     # core.build_reaction_dict and core.parse_parameters index by the order of
#     # stoichiometry_matrix(); time_series_dataframe uses the order of
#     # rn.species(). Both orderings and the mapping between them are built here
#     # to prevent a silent misalignment from assigning parameters to the wrong
#     # reaction.
#     sm = rn.stoichiometry_matrix()
#     mat_species = list(sm.species)
#     mat_reactions = list(sm.reactions)

#     out_species = [s.name for s in rn.species()]
#     out_reactions = [r.name() for r in rn.reactions()]

#     if sorted(mat_species) != sorted(out_species):
#         raise ValueError("Mismatch between species() and stoichiometry_matrix().species")
#     if sorted(mat_reactions) != sorted(out_reactions):
#         raise ValueError("Mismatch between reactions() and stoichiometry_matrix().reactions")

#     n_species = len(mat_species)
#     n_reactions = len(mat_reactions)

#     reactant_matrix = np.asarray(rn.reactants_matrix(), dtype=float)
#     stoich_matrix = np.asarray(sm, dtype=float)

#     rn_dict = build_reaction_dict(rn)
#     rate = validate_rate_list(rate, n_reactions)

#     # --- Warnings about SSA validity --------------------------------------
#     non_exact = sorted({r for r in rate if r not in _EXACT_LAWS})
#     time_dep = sorted({r for r in rate if r in _TIME_DEPENDENT})
#     if non_exact:
#         warnings.warn(
#             f"Kinetic laws that are not exact for the SSA: {non_exact}. The "
#             f"Gillespie algorithm is exact only with 'mak'; the rest are "
#             f"evaluated as heuristic propensities.",
#             RuntimeWarning, stacklevel=2)
#     if time_dep:
#         warnings.warn(
#             f"Time-dependent laws: {time_dep}. The derivation of "
#             f"tau = -ln(u)/a0 assumes a constant propensity between events; "
#             f"the resulting trajectory is not exact.",
#             RuntimeWarning, stacklevel=2)

#     # --- Initial conditions ------------------------------------------------
#     state = np.zeros(n_species, dtype=float)
#     if x0 is None:
#         x0_vals = generate_random_vector(
#             n_species, seed=seed, min_value=50, max_value=100.0).astype(int).tolist()
#         state = np.array(x0_vals, dtype=float)          # mat_species ordering
#         if verbose: 
#             print(f"\nx0 = {dict(zip(mat_species, x0_vals))}")
#     elif isinstance(x0, dict):
#         for name, qty in x0.items():
#             if name not in mat_species:
#                 raise ValueError(f"Species '{name}' does not exist in the network (species: {mat_species})")
#             state[mat_species.index(name)] = qty
#     else:
#         x0_arr = np.asarray(x0, dtype=float)
#         if x0_arr.shape[0] != n_species:
#             raise ValueError(
#                 f"x0 must have {n_species} elements, received {x0_arr.shape[0]}")
#         state = x0_arr.copy()

#     if np.any(state < 0):
#         raise ValueError("x0 cannot contain negative populations")
#     state = np.round(state)     # the SSA operates on integer counts

#     # --- Parameters --------------------------------------------------------
#     if spec_vector is None:
#         spec_vector = generate_default_parameters_seeded(
#             rate, n_reactions, additional_laws, seed=seed)
#     if len(spec_vector) != n_reactions:
#         raise ValueError(
#             f"spec_vector must have {n_reactions} entries, received {len(spec_vector)}")
#     parameters = parse_parameters(rn, rn_dict, rate, spec_vector)
#     if verbose:
#         print("spec_vector =", spec_vector)

#     # --- Kinetic law registry ----------------------------------------------
#     rate_laws = dict(KINETIC_REGISTRY)
#     if additional_laws:
#         for name in rate:
#             if name not in rate_laws:
#                 if name in additional_laws:
#                     rate_laws[name] = additional_laws[name]
#                 else:
#                     raise NotImplementedError(
#                         f"Kinetic law '{name}' not defined or registered")
#     for name in rate:
#         if name not in rate_laws:
#             raise NotImplementedError(f"Kinetic law '{name}' not defined or registered")

#     species_idx = {s: i for i, s in enumerate(mat_species)}

#     if verbose:
#         print_differential_equations(rn, rn_dict)
#         print_velocity_expressions(rn, rn_dict, rate, rate_laws, additional_laws)

#     # --- Propensities ------------------------------------------------------
#     def propensities(x: np.ndarray, t: float) -> np.ndarray:
#         a = np.zeros(n_reactions)
#         for j, reaction in enumerate(mat_reactions):
#             kinetic = rate[j]
#             reactants = rn_dict[reaction][0]

#             if kinetic == 'mak' and exact_mass_action:
#                 a[j] = _propensity_mak_exact(x, reactant_matrix[:, j],
#                                              float(parameters[reaction][0]))
#                 continue

#             if kinetic == 'mmk' and not reactants:
#                 continue

#             # Zero propensity if any reactant falls below its coefficient:
#             # without this, saturating laws can fire reactions with no
#             # substrate available.
#             if np.any(x < reactant_matrix[:, j]):
#                 continue

#             fn = rate_laws[kinetic]
#             if kinetic in _TIME_DEPENDENT:
#                 a[j] = fn(reactants, x, species_idx, parameters[reaction], t=t)
#             else:
#                 a[j] = fn(reactants, x, species_idx, parameters[reaction])
#         return np.maximum(a, 0.0)

#     # --- SSA loop (direct method) ------------------------------------------
#     t0, t_end = float(t_span[0]), float(t_span[1])
#     t = t0
#     times = [t]
#     populations = [state.copy()]
#     reaction_counts = {name: 0 for name in mat_reactions}

#     for _ in range(max_iter):
#         a = propensities(state, t)
#         a0 = a.sum()
#         if a0 <= 0:
#             break                                   # absorbing state

#         # u in (0, 1]: rng.random() returns [0, 1), which admits 0 and would
#         # blow up log(1/u).
#         u1 = 1.0 - rng.random()
#         u2 = rng.random()

#         tau = -math.log(u1) / a0
#         if t + tau > t_end:
#             break

#         j = int(np.searchsorted(np.cumsum(a), u2 * a0, side="right"))
#         j = min(j, n_reactions - 1)

#         state = state + stoich_matrix[:, j]
#         state = np.maximum(state, 0.0)
#         t += tau
#         reaction_counts[mat_reactions[j]] += 1

#         times.append(t)
#         populations.append(state.copy())

#         if stop_condition is not None and stop_condition(state, mat_species):
#             break

#     times = np.array(times)
#     result = np.array(populations)                  # (n_events+1, n_species)

#     # --- Optional resampling onto a uniform grid ---------------------------
#     if n_steps is not None:
#         grid = np.linspace(t0, t_end, n_steps)
#         # Zero-order hold: the state is constant between events, so this does
#         # not interpolate anything, it merely reindexes the same trajectory.
#         idx = np.searchsorted(times, grid, side="right") - 1
#         idx = np.clip(idx, 0, len(times) - 1)
#         result = result[idx]
#         times = grid

#     # --- Reorder to rn.species() ordering for the DataFrames ---------------
#     reorder = [mat_species.index(s) for s in out_species]
#     result_out = result[:, reorder]

#     time_series_df = time_series_dataframe(rn, result_out, times)
#     flux_vector_df = flux_vector_dataframe(
#         rn, rn_dict, rate, result_out, times, spec_vector, rate_laws)

#     return time_series_df, flux_vector_df


# __all__ = ['gillespie', 'generate_default_parameters_seeded']