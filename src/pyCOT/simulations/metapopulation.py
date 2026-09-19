"""
Metapopulation Simulation Module
=================================

Hybrid deterministic/stochastic metapopulation model with:

    local reactions + inter-patch transport

Mathematical model
------------------

For species s and patch p:

    dx_s^(p)/dt =
        reaction_s^(p)(x^(p))
        + sum_q (L_s)[p,q] x_s^(q)

where L_s is the species-specific transport Laplacian:

    (L_s)[p,q] = T_s[q -> p],        p != q

    (L_s)[p,p] = -sum_{q != p} T_s[p -> q]

Every column of L_s sums to zero, hence 1^T L_s = 0 and transport
conserves the total amount of each species. Only reactions create
or destroy mass.

Reporting (verbose=True)
------------------------
`print_simulation_setup` runs BEFORE the integration loop and prints the
complete configuration: patches, modes, kinetic laws and parameters,
initial conditions, connectivity, transport rates, Laplacians and
migration matrices, plus the effective time discretization.

`_print_full_system` (also before the loop) prints, per patch:

  [1] the local reaction equations and velocity expressions,
  [2] the transport coupling terms per species,
  [3] the full coupled equation (local + transport) per species.

Parameter NAMES are shown symbolically (k_R, Vmax_R, Km_R, Kd_R, n_R),
not their numerical values, so the report documents the structure of
the system rather than a particular parametrization.

`print_simulation_summary` runs after the loop and reports the final
state and the mass balance per species.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy.linalg import expm

from pyCOT.simulations.ode import simulation
from pyCOT.simulations.stochastic import gillespie
from pyCOT.simulations.core import (
    generate_default_parameters,
    validate_rate_list,
    build_reaction_dict,
    parse_parameters,
    print_differential_equations,
    print_velocity_expressions,
)

from pyCOT.kinetics import KINETIC_REGISTRY


# ============================================================================
# Connectivity templates (public helpers)
# ============================================================================

def uniform_connectivity(num_patches: int, p_stay: float = 0.5) -> np.ndarray:
    """Well-mixed template: every patch disperses equally to all others."""
    if num_patches < 1:
        raise ValueError("num_patches must be >= 1.")
    if not 0.0 <= p_stay <= 1.0:
        raise ValueError("p_stay must lie in [0, 1].")

    if num_patches == 1:
        return np.ones((1, 1))

    off = (1.0 - p_stay) / (num_patches - 1)
    C = np.full((num_patches, num_patches), off, dtype=float)
    np.fill_diagonal(C, p_stay)
    return C


def ring_connectivity(
    num_patches: int,
    p_stay: float = 0.5,
    directed: bool = False,
) -> np.ndarray:
    """Ring template. directed=True -> advective ring."""
    if num_patches < 2:
        raise ValueError("A ring needs at least 2 patches.")
    if not 0.0 <= p_stay <= 1.0:
        raise ValueError("p_stay must lie in [0, 1].")

    C = np.zeros((num_patches, num_patches), dtype=float)
    move = 1.0 - p_stay

    for p in range(num_patches):
        C[p, p] = p_stay
        if directed:
            C[p, (p + 1) % num_patches] += move
        else:
            C[p, (p + 1) % num_patches] += move / 2.0
            C[p, (p - 1) % num_patches] += move / 2.0

    return C


def chain_connectivity(num_patches: int, p_stay: float = 0.5) -> np.ndarray:
    """Open 1D chain."""
    if num_patches < 2:
        raise ValueError("A chain needs at least 2 patches.")

    C = np.zeros((num_patches, num_patches), dtype=float)
    move = 1.0 - p_stay

    for p in range(num_patches):
        neighbours = [q for q in (p - 1, p + 1) if 0 <= q < num_patches]
        C[p, p] = p_stay
        for q in neighbours:
            C[p, q] = move / len(neighbours)

    return C


# ============================================================================
# Internal helpers
# ============================================================================

def _as_patch_list(value, num_patches: int, name: str) -> List:
    if isinstance(value, (list, tuple)) and not isinstance(value, str):
        if len(value) != num_patches:
            raise ValueError(
                f"'{name}' has {len(value)} entries for "
                f"{num_patches} patches."
            )
        return list(value)
    return [value] * num_patches


def _rate_per_patch(rate, num_patches: int) -> List:
    if isinstance(rate, str):
        return [rate] * num_patches

    if isinstance(rate, (list, tuple)):
        nested = all(
            isinstance(r, (list, tuple)) and not isinstance(r, str)
            for r in rate
        )
        if nested and len(rate) == num_patches:
            return [list(r) for r in rate]
        return [list(rate)] * num_patches

    raise ValueError(
        "'rate' must be a string, a list of laws, or a list of "
        "per-patch law lists."
    )


def _normalize_connectivity(
    C: Optional[np.ndarray],
    num_patches: int,
    verbose: bool = True,
) -> np.ndarray:
    if C is None:
        return np.eye(num_patches)

    C = np.asarray(C, dtype=float)

    if C.shape != (num_patches, num_patches):
        raise ValueError(
            f"connectivity_matrix must have shape "
            f"({num_patches}, {num_patches}), got {C.shape}"
        )

    if np.any(C < 0):
        raise ValueError("connectivity_matrix cannot contain negative values.")

    row_sums = C.sum(axis=1)

    if np.any(row_sums <= 0):
        raise ValueError(
            "Every row of connectivity_matrix must have a positive sum."
        )

    if not np.allclose(row_sums, 1.0, atol=1e-10):
        if verbose:
            print("Warning: normalizing connectivity_matrix so rows sum to 1.")
        C = C / row_sums[:, None]

    return C


def _build_laplacian_from_connectivity(
    C: np.ndarray,
    D: float,
    present: np.ndarray,
) -> np.ndarray:
    n = C.shape[0]
    L = np.zeros((n, n), dtype=float)

    if D <= 0:
        return L

    idx = np.flatnonzero(present)

    if len(idx) < 2:
        return L

    Cs = D * C[np.ix_(idx, idx)]

    L_block = Cs.T.copy()

    diag = -Cs.sum(axis=1) + np.diag(Cs)
    L_block[np.diag_indices_from(L_block)] = diag

    L[np.ix_(idx, idx)] = L_block
    return L


def _build_laplacian_from_rates(
    rates: np.ndarray,
    present: np.ndarray,
) -> np.ndarray:
    rates = np.asarray(rates, dtype=float)
    n = rates.shape[0]

    if rates.shape != (n, n):
        raise ValueError("Transport-rate matrix must be square.")

    if np.any(rates < 0):
        raise ValueError("Transport rates cannot be negative.")

    if not np.allclose(np.diag(rates), 0.0):
        raise ValueError("The diagonal of a transport-rate matrix must be zero.")

    L = np.zeros_like(rates)
    idx = np.flatnonzero(present)

    if len(idx) < 2:
        return L

    Rs = rates[np.ix_(idx, idx)]

    L_block = Rs.T.copy()
    diag = -Rs.sum(axis=1) + np.diag(Rs)
    L_block[np.diag_indices_from(L_block)] = diag

    L[np.ix_(idx, idx)] = L_block
    return L


def _migration_matrix(L: np.ndarray, dt: float) -> np.ndarray:
    if dt < 0:
        raise ValueError("dt must be non-negative.")

    n = L.shape[0]

    if dt == 0.0 or np.allclose(L, 0.0):
        return np.eye(n)

    M = L.T * dt
    P = expm(M)

    P[np.abs(P) < 1e-14] = 0.0
    P = np.clip(P, 0.0, None)

    row_sums = P.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    P = P / row_sums

    return P


def _stochastic_round(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    floor = np.floor(x)
    remainder = x - floor
    rounded = floor + (rng.random(x.shape) < remainder)
    return rounded.astype(float)


def _check_invariants(
    L_species: Dict[str, np.ndarray],
    P_mig: Dict[str, np.ndarray],
    tol: float = 1e-10,
) -> None:
    for s, L in L_species.items():
        col = L.sum(axis=0)
        if not np.allclose(col, 0.0, atol=tol):
            raise AssertionError(f"L_{s} column sums deviate from 0: {col}")

    for s, P in P_mig.items():
        rs = P.sum(axis=1)
        if not np.allclose(rs, 1.0, atol=tol):
            raise AssertionError(f"P_{s} row sums deviate from 1: {rs}")
        if np.any(P < -tol):
            raise AssertionError(f"P_{s} has negative entries.")


# ============================================================================
# Symbolic helpers for the full-system report
# ============================================================================

def _velocity_string(kinetic: str, reactants, r_name: str, sub: int) -> str:
    """
    Return a symbolic string for the velocity of a single reaction,
    using the symbolic PARAMETER NAMES (k_R, Vmax_R, Km_R, Kd_R, n_R)
    instead of their numerical values.

    This matches the notation used by `print_velocity_expressions` from
    the ODE backend, so that the printed full system is expressed in
    terms of the parameters, not their particular values.

    `sub` is the 1-based patch index used as subscript for the species
    names, e.g. 'A_1' for patch 1.

    Unknown laws fall back to a placeholder 'v_<reaction>' so the
    report never fails.
    """
    def sname(sp):
        return f"{sp}_{sub}"

    if kinetic == "mak":
        k = f"k_{r_name}"
        if reactants:
            parts = []
            for sp, coef in reactants:
                if coef == 1:
                    parts.append(sname(sp))
                else:
                    parts.append(f"{sname(sp)}^{int(coef)}")
            return f"{k}*" + "*".join(parts)
        return k

    if kinetic == "mmk":
        Vmax = f"Vmax_{r_name}"
        Km = f"Km_{r_name}"
        sp = reactants[0][0] if reactants else "0"
        return f"{Vmax}*{sname(sp)}/({Km}+{sname(sp)})"

    if kinetic == "hill":
        Vmax = f"Vmax_{r_name}"
        Kd = f"Kd_{r_name}"
        n = f"n_{r_name}"
        sp = reactants[0][0] if reactants else "0"
        return (f"{Vmax}*{sname(sp)}^{n}/"
                f"({Kd}^{n}+{sname(sp)}^{n})")

    return f"v_{r_name}"


def _local_rhs_terms(rn, rn_dict, rate_list, spec_list, patch_idx: int):
    """
    Build the symbolic local reaction RHS for one patch.

    Returns a dict {species_name: list_of_signed_terms}, where each term
    is a string such as '-1.0*(A_1)' or '+0.3*(A_1*B_1)'. Species names
    carry the patch subscript (1-based) so the terms can be concatenated
    directly with the transport terms.

    Parameters (spec_list) are only used to keep the signature aligned
    with the setup reporting; the printed velocity uses symbolic
    parameter names, not their numerical values.
    """
    sub = patch_idx + 1
    terms = {s.name: [] for s in rn.species()}

    for r_idx, reaction in enumerate(rn.reactions()):
        r_name = reaction.name()
        kinetic = rate_list[r_idx]
        reactants = rn_dict[r_name][0]
        products = rn_dict[r_name][1]

        v = _velocity_string(kinetic, reactants, r_name, sub)

        for sp, coef in reactants:
            c = -int(coef)
            if c == -1:
                terms[sp].append(f"-({v})")
            else:
                terms[sp].append(f"{c}*({v})")

        for sp, coef in products:
            c = int(coef)
            if c == 1:
                terms[sp].append(f"+({v})")
            else:
                terms[sp].append(f"+{c}*({v})")

    return terms


def _transport_terms(L: np.ndarray, s: str, p: int):
    """
    Symbolic transport terms for species `s` in patch `p`:
        sum_q L[p, q] * s_q
    Returns a list of strings like ['-0.02*A_1', '+0.02*A_2'].
    """
    parts = []
    for q in range(L.shape[0]):
        coeff = L[p, q]
        if abs(coeff) > 0:
            sign = "+" if coeff >= 0 else "-"
            parts.append(f"{sign} {abs(coeff):.4g}*{s}_{q + 1}")
    return parts


def _clean_join(terms):
    """Join a list of signed terms, removing the leading '+ '."""
    if not terms:
        return "0"
    s = " ".join(terms)
    if s.startswith("+ "):
        s = s[2:]
    return s


# ============================================================================
# Reporting helpers
# ============================================================================

def _tabla_por_parche(nombres: List[str], M: np.ndarray, present=None,
                      ancho: int = 10, fmt: str = "{:.4g}") -> str:
    m = M.shape[1]
    lineas = ["    " + "".join(f"{'P' + str(p + 1):>{ancho}}" for p in range(m))]

    for i, nombre in enumerate(nombres):
        fila = f"  {nombre:<2}"
        for p in range(m):
            if present is not None and not present[i, p]:
                fila += f"{'-':>{ancho}}"
            else:
                fila += f"{fmt.format(M[i, p]):>{ancho}}"
        lineas.append(fila)

    return "\n".join(lineas)


def _matriz(M: np.ndarray, sangria: str = "    ", precision: int = 4) -> str:
    txt = np.array2string(M, precision=precision, suppress_small=True)
    return "\n".join(sangria + linea for linea in txt.splitlines())


def _rates_desde_laplaciano(L: np.ndarray) -> np.ndarray:
    T = L.T.copy()
    np.fill_diagonal(T, 0.0)
    return T


def _print_full_system(
    rns,
    rate_per_patch,
    spec_per_patch,
    additional_laws,
    species,
    patch_species,
    patch_reactions,
    L_species,
    moves,
    verbose: bool = True,
) -> None:
    """
    Print the complete metapopulation system in three sections:

      [1] Local reaction equations and velocity expressions per patch.
      [2] Transport coupling terms per species and patch.
      [3] Full coupled equation (local + transport) per species and patch.

    Local velocities are shown with symbolic parameter names, not values.
    """
    if not verbose:
        return

    ancho = 74
    num_patches = len(rns)

    print("\n" + "=" * ancho)
    print("SISTEMA COMPLETO: REACCIONES LOCALES + TRANSPORTE")
    print("=" * ancho)

    rate_laws = dict(KINETIC_REGISTRY)
    if additional_laws:
        rate_laws.update(additional_laws)

    homogeneo_local = all(
        rate_per_patch[p] == rate_per_patch[0]
        and [list(np.atleast_1d(v)) for v in spec_per_patch[p]]
        == [list(np.atleast_1d(v)) for v in spec_per_patch[0]]
        and patch_reactions[p] == patch_reactions[0]
        for p in range(num_patches)
    )

    # ==================================================================
    # [1] Local reactions
    # ==================================================================
    print("\n" + "-" * ancho)
    print("[1] REACCIONES LOCALES")
    print("-" * ancho)

    parches_locales = [0] if homogeneo_local else range(num_patches)

    if homogeneo_local and num_patches > 1:
        print(f"  (identicas en todos los {num_patches} parches)")

    for p in parches_locales:
        if not homogeneo_local:
            print(f"\n  --- Parche {p + 1} ---")
        rn_dict = build_reaction_dict(rns[p])
        _ = parse_parameters(
            rns[p], rn_dict, rate_per_patch[p], spec_per_patch[p]
        )
        print(f"\n  spec_vector (parche {p + 1}) = {spec_per_patch[p]}")
        print_differential_equations(rns[p], rn_dict)
        print()
        print_velocity_expressions(
            rns[p], rn_dict, rate_per_patch[p], rate_laws, additional_laws
        )

    # ==================================================================
    # [2] Transport coupling
    # ==================================================================
    print("\n" + "-" * ancho)
    print("[2] TRANSPORTE ENTRE PARCHES")
    print("-" * ancho)
    print("  Para la especie s y el parche p:")
    print("    dx_s^(p)/dt |_transporte = sum_q L_s[p,q] * x_s^(q)")
    print()

    for p in range(num_patches):
        print(f"  Parche {p + 1}:")
        found_any = False
        for s in species:
            if not moves.get(s, False):
                continue
            if s not in patch_species[p]:
                continue
            L = L_species[s]
            parts = _transport_terms(L, s, p)
            print(f"    d{s}_{p + 1}/dt = {_clean_join(parts)}")
            found_any = True
        if not found_any:
            print("    (sin transporte)")

    # ==================================================================
    # [3] Full coupled system
    # ==================================================================
    print("\n" + "-" * ancho)
    print("[3] SISTEMA COMPLETO POR PARCHE")
    print("-" * ancho)
    print("  Cada ecuacion es la suma de la parte local [1] mas el")
    print("  transporte [2], con la especie subindicada por parche.")
    print()

    for p in range(num_patches):
        print(f"  Parche {p + 1}:")

        rn_dict = build_reaction_dict(rns[p])
        local_terms = _local_rhs_terms(
            rns[p], rn_dict, rate_per_patch[p], spec_per_patch[p], p
        )

        for s in species:
            if s not in patch_species[p]:
                continue

            lhs = f"d{s}_{p + 1}/dt"

            if local_terms.get(s):
                local_str = _clean_join(local_terms[s])
            else:
                local_str = "0"

            L = L_species[s]
            transport_str = _clean_join(_transport_terms(L, s, p))

            if transport_str == "0":
                rhs = local_str
            elif local_str == "0":
                rhs = transport_str
            else:
                rhs = f"{local_str} + {transport_str}"

            print(f"    {lhs} = {rhs}")

    print("\n" + "=" * ancho)


def print_simulation_setup(
    *, num_patches, modes, rn_is_list, species, reactions, present,
    patch_species, patch_reactions, rate_per_patch, spec_per_patch,
    X0, t_span, n_steps, dt_out, dt_couple, dt, n_sub,
    C, connectivity_given, D_dict, D_default, transport_rates,
    L_species, P_mig, moves, method, rtol, atol, seed,
    max_patches_matrices: int = 8,
):
    """Print the full configuration before integrating."""

    ancho = 74
    print("\n" + "=" * ancho)
    print("METAPOPULATION SIMULATION - CONFIGURACION")
    print("=" * ancho)

    print(f"Parches ............ {num_patches}")
    print(f"Modos .............. {modes}")
    print(f"Redes .............. {'heterogeneas' if rn_is_list else 'compartida'}")
    print(f"Especies ........... {species}")
    print(f"Reacciones ......... {reactions}")
    print(f"Semilla ............ {seed}")
    if any(mo == "ode" for mo in modes):
        print(f"Integrador ODE ..... {method}  (rtol={rtol:g}, atol={atol:g})")

    print("\n" + "-" * ancho)
    print("DISCRETIZACION TEMPORAL")
    print("-" * ancho)
    print(f"  t_span ........... ({t_span[0]:g}, {t_span[1]:g})")
    print(f"  n_steps .......... {n_steps}     dt_out = {dt_out:g}")
    print(f"  dt_couple pedido . {dt_couple:g}")
    print(f"  dt efectivo ...... {dt:g}     ({n_sub} subpaso(s) por intervalo)")

    if dt_couple > dt_out * (1 + 1e-12):
        print(f"  AVISO: dt_couple ({dt_couple:g}) es mayor que dt_out "
              f"({dt_out:g}).")
        print(f"         El paso efectivo es min(dt_couple, dt_out) = {dt:g}, "
              f"porque")
        print(f"         n_sub = ceil(dt_out/dt_couple) nunca baja de 1. "
              f"Subir dt_couple")
        print(f"         por encima de dt_out no cambia nada.")

    print("\n" + "-" * ancho)
    print("LEYES CINETICAS Y PARAMETROS")
    print("-" * ancho)

    homogeneo = all(
        rate_per_patch[p] == rate_per_patch[0]
        and [list(np.atleast_1d(v)) for v in spec_per_patch[p]]
        == [list(np.atleast_1d(v)) for v in spec_per_patch[0]]
        and patch_reactions[p] == patch_reactions[0]
        for p in range(num_patches)
    )

    parches_a_mostrar = [0] if homogeneo else range(num_patches)

    if homogeneo and num_patches > 1:
        print("  (identicos en todos los parches)")

    for p in parches_a_mostrar:
        if not homogeneo:
            print(f"\n  Parche {p + 1}:")
        for j, r in enumerate(patch_reactions[p]):
            params = ", ".join(
                f"{v:g}" for v in np.atleast_1d(spec_per_patch[p][j])
            )
            print(f"    {r:<6} ley={rate_per_patch[p][j]:<6} params=[{params}]")

    print("\n" + "-" * ancho)
    print("CONDICIONES INICIALES  x0  (filas = especies, columnas = parches)")
    print("-" * ancho)
    print(_tabla_por_parche(species, X0, present))

    totales = np.array([
        np.nansum(np.where(present[i], X0[i], np.nan))
        for i in range(len(species))
    ])
    print("\n  Total por especie: "
          + ", ".join(f"{s}={t:.4g}" for s, t in zip(species, totales)))

    if not present.all():
        print("\n  Matriz de presencia (1 = la especie existe en ese parche):")
        print(_tabla_por_parche(species, present.astype(int), fmt="{:d}"))

    print("\n" + "-" * ancho)
    print("TRANSPORTE")
    print("-" * ancho)

    if transport_rates:
        print(f"  Especies por transport_rates: {sorted(transport_rates)}")

    otras = [s for s in species
             if not (transport_rates and s in transport_rates)]

    if otras:
        print(f"  Especies por D_dict + connectivity_matrix: {otras}")
        print(f"  D_dict = {D_dict}   D_default = {D_default:g}")

        if not connectivity_given:
            print("  connectivity_matrix = None -> se usa la identidad, "
                  "que NO produce transporte.")
        else:
            print("\n  Matriz de conectividad C (filas suman 1, "
                  "diagonal = retencion):")
            print(_matriz(C))
            for s in otras:
                D = float(D_dict.get(s, D_default))
                if D > 0:
                    tasas = D * (1.0 - np.diag(C))
                    print(f"    {s}: D = {D:g}  ->  tasa de emigracion "
                          f"D(1 - C_pp) = "
                          f"{np.array2string(tasas, precision=4, suppress_small=True)}")

    for s in species:
        L = L_species[s]
        print(f"\n  --- especie {s} ---")

        if not moves[s]:
            print("    sin transporte (L = 0)")
            continue

        if num_patches <= max_patches_matrices:
            print("    T[p -> q] (fila = origen, columna = destino):")
            print(_matriz(_rates_desde_laplaciano(L), sangria="      "))
            print("    L (columna q = desde q; fila p = hacia p):")
            print(_matriz(L, sangria="      "))
        else:
            print(f"    ({num_patches} parches: matrices omitidas)")

        print("    tasa de salida por parche -diag(L): "
              f"{np.array2string(-np.diag(L), precision=4, suppress_small=True)}")
        print("    suma de columnas de L: "
              f"{np.array2string(L.sum(axis=0), precision=2, suppress_small=True)}"
              "   (debe ser 0: conservacion)")
        print("    tipo: " + ("simetrico -> difusion pura"
                              if np.allclose(L, L.T)
                              else "asimetrico -> transporte dirigido (deriva)"))

        if not np.allclose(L.sum(axis=0), 0.0, atol=1e-10):
            print("    ATENCION: el transporte NO conserva masa.")

        if num_patches <= max_patches_matrices:
            P = P_mig[s]
            print(f"    P = expm(L^T dt) con dt = {dt:g} "
                  "(fila p = destino de lo que sale de p):")
            print(_matriz(P, sangria="      "))
            print("    suma de filas de P: "
                  f"{np.array2string(P.sum(axis=1), precision=6)}"
                  "   (debe ser 1)")

    print("=" * ancho)


def print_simulation_summary(species, X0, X_final, present):
    ancho = 74
    print("\n" + "=" * ancho)
    print("METAPOPULATION SIMULATION - RESUMEN")
    print("=" * ancho)

    print("Estado final (filas = especies, columnas = parches):")
    print(_tabla_por_parche(species, X_final, present))

    print("\nBalance de masa por especie "
          "(variacion = reacciones + redondeo hibrido):")
    print(f"  {'especie':<8}{'inicial':>12}{'final':>12}{'variacion':>14}")

    for i, s in enumerate(species):
        mi = np.nansum(np.where(present[i], X0[i], np.nan))
        mf = np.nansum(np.where(present[i], X_final[i], np.nan))
        print(f"  {s:<8}{mi:12.4g}{mf:12.4g}{mf - mi:+14.4g}")

    print("=" * ancho)


# ============================================================================
# Post-processing helper
# ============================================================================

def total_mass(X_out: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    return {s: np.nansum(v, axis=1) for s, v in X_out.items()}


# ============================================================================
# Main simulation
# ============================================================================

def simulate_metapopulation_dynamics(
    rn: Union[object, Sequence],
    rate: Union[str, List] = "mak",
    modes: Union[str, List[str]] = "ode",
    num_patches: Optional[int] = None,

    D_dict: Optional[Dict[str, float]] = None,
    D_default: float = 0.0,
    connectivity_matrix: Optional[np.ndarray] = None,

    transport_rates: Optional[Dict[str, np.ndarray]] = None,

    x0_dict: Optional[Dict[str, np.ndarray]] = None,
    spec_vector: Optional[Union[List, Sequence[List]]] = None,

    t_span: Tuple[float, float] = (0, 20),
    n_steps: int = 500,
    dt_couple: Optional[float] = None,

    additional_laws: Optional[Dict] = None,
    method: str = "LSODA",

    rtol: float = 1e-8,
    atol: float = 1e-10,

    seed: Optional[int] = None,
    verbose: bool = True,
    return_info: bool = False,
    check_invariants: bool = False,
):
    """
    Simulate a hybrid heterogeneous metapopulation.

    See the module docstring for the mathematical model and the
    reporting behaviour. When `verbose=True`, the function prints:

      * the configuration (print_simulation_setup),
      * the complete system per patch (_print_full_system),
      * the final state and mass balance (print_simulation_summary).

    Parameter names in the full-system report are symbolic
    (k_R, Vmax_R, Km_R, Kd_R, n_R), not their numerical values.
    """

    rng = np.random.default_rng(seed)

    # 1. Patches and reaction networks
    rn_is_list = isinstance(rn, (list, tuple))

    if num_patches is None:
        if rn_is_list:
            num_patches = len(rn)
        elif connectivity_matrix is not None:
            num_patches = np.asarray(connectivity_matrix).shape[0]
        elif transport_rates:
            first_matrix = next(iter(transport_rates.values()))
            num_patches = np.asarray(first_matrix).shape[0]
        elif isinstance(modes, (list, tuple)):
            num_patches = len(modes)
        else:
            num_patches = 3

    if num_patches < 1:
        raise ValueError("num_patches must be >= 1.")

    rns = list(rn) if rn_is_list else [rn] * num_patches

    if len(rns) != num_patches:
        raise ValueError(f"Got {len(rns)} networks for {num_patches} patches.")

    # 2. Simulation mode
    modes = _as_patch_list(modes, num_patches, "modes")

    for mode in modes:
        if mode not in ("ode", "ssa"):
            raise ValueError(f"mode must be 'ode' or 'ssa', got {mode!r}.")

    # 3. Species and reactions
    patch_species = [[s.name for s in network.species()] for network in rns]
    patch_reactions = [
        [reaction.name() for reaction in network.reactions()]
        for network in rns
    ]

    species: List[str] = []
    for species_list in patch_species:
        for s in species_list:
            if s not in species:
                species.append(s)

    reactions: List[str] = []
    for reaction_list in patch_reactions:
        for r in reaction_list:
            if r not in reactions:
                reactions.append(r)

    species_index = {s: i for i, s in enumerate(species)}
    patch_species_idx = [
        np.array([species_index[s] for s in patch_species[p]], dtype=int)
        for p in range(num_patches)
    ]

    # 4. Presence matrix
    present = np.zeros((len(species), num_patches), dtype=bool)

    for p, species_list in enumerate(patch_species):
        for s in species_list:
            present[species_index[s], p] = True

    # 5. Kinetic laws and parameters
    rate_per_patch = _rate_per_patch(rate, num_patches)

    rate_per_patch = [
        validate_rate_list(rate_per_patch[p], len(patch_reactions[p]))
        for p in range(num_patches)
    ]

    if spec_vector is None:
        spec_per_patch = [
            generate_default_parameters(
                rate_per_patch[p], len(patch_reactions[p]), additional_laws,
            )
            for p in range(num_patches)
        ]

    elif (
        isinstance(spec_vector, (list, tuple))
        and len(spec_vector) == num_patches
        and all(
            isinstance(sv, (list, tuple, np.ndarray))
            for sv in spec_vector
        )
        and all(
            len(sv) == len(patch_reactions[p])
            for p, sv in enumerate(spec_vector)
        )
    ):
        spec_per_patch = [list(sv) for sv in spec_vector]

    else:
        shared = list(spec_vector)
        for p in range(num_patches):
            if len(shared) != len(patch_reactions[p]):
                raise ValueError(
                    f"Shared spec_vector has {len(shared)} entries but "
                    f"patch {p} has {len(patch_reactions[p])} reactions. "
                    f"With heterogeneous networks, pass one spec_vector "
                    f"per patch."
                )
        spec_per_patch = [shared] * num_patches

    # 6. Connectivity
    C = _normalize_connectivity(connectivity_matrix, num_patches, verbose)

    if D_dict is None:
        D_dict = {}

    unknown = [s for s in D_dict if s not in species_index]
    if unknown:
        raise ValueError(
            f"D_dict refers to species absent from every network: {unknown}. "
            f"Known species: {species}"
        )

    if transport_rates:
        unknown = [s for s in transport_rates if s not in species_index]
        if unknown:
            raise ValueError(
                f"transport_rates refers to unknown species: {unknown}. "
                f"Known species: {species}"
            )

    uses_D = any(float(v) > 0 for v in D_dict.values()) or D_default > 0
    rates_cover_all = (
        transport_rates is not None
        and all(s in transport_rates for s in species)
    )

    if (
        uses_D
        and connectivity_matrix is None
        and num_patches > 1
        and not rates_cover_all
    ):
        if verbose:
            print(
                "Warning: dispersal coefficients were given but "
                "connectivity_matrix is None and transport_rates does not "
                "cover every species. The default identity template has no "
                "off-diagonal entries, so no transport will occur. Pass a "
                "connectivity_matrix (uniform_connectivity / ring_connectivity "
                "/ chain_connectivity) or transport_rates."
            )

    # 7. Initial conditions
    any_stochastic = any(mode == "ssa" for mode in modes)

    X = np.zeros((len(species), num_patches), dtype=float)

    if x0_dict is None:
        for i, s in enumerate(species):
            if any_stochastic:
                X[i] = rng.integers(1, 11, size=num_patches)
            else:
                X[i] = np.round(rng.uniform(0, 2.0, size=num_patches), 2)

    else:
        for s, values in x0_dict.items():
            if s not in species_index:
                raise ValueError(f"Species {s!r} is not present in any network.")

            values = np.atleast_1d(np.asarray(values, dtype=float))

            if values.size == 1:
                values = np.repeat(values, num_patches)

            if values.size != num_patches:
                raise ValueError(
                    f"x0 for {s!r} must have length {num_patches}, "
                    f"got {values.size}."
                )

            X[species_index[s]] = values

    if np.any(X < 0):
        raise ValueError("Initial conditions cannot be negative.")

    X[~present] = 0.0

    for p, mode in enumerate(modes):
        if mode == "ssa":
            X[:, p] = np.round(X[:, p])

    X0 = X.copy()

    # 8. Time discretization
    if n_steps < 2:
        raise ValueError("n_steps must be >= 2.")

    t_grid = np.linspace(t_span[0], t_span[1], n_steps)
    dt_out = (t_span[1] - t_span[0]) / max(n_steps - 1, 1)

    if dt_out <= 0:
        raise ValueError("t_span must be increasing.")

    if dt_couple is None:
        dt_couple = dt_out

    if dt_couple <= 0:
        raise ValueError("dt_couple must be positive.")

    n_sub = max(1, int(np.ceil(dt_out / dt_couple)))
    dt = dt_out / n_sub

    # 9. Species-specific transport Laplacians
    L_species: Dict[str, np.ndarray] = {}

    for i, s in enumerate(species):

        if transport_rates is not None and s in transport_rates:
            rates = np.asarray(transport_rates[s], dtype=float)

            if rates.shape != (num_patches, num_patches):
                raise ValueError(
                    f"Transport matrix for species {s!r} must have shape "
                    f"({num_patches}, {num_patches}), got {rates.shape}."
                )

            L_species[s] = _build_laplacian_from_rates(rates, present[i])

        else:
            D = float(D_dict.get(s, D_default))
            if D < 0:
                raise ValueError(
                    f"Dispersal coefficient for {s!r} is negative."
                )
            L_species[s] = _build_laplacian_from_connectivity(C, D, present[i])

    # 10. Exact migration matrices
    P_mig: Dict[str, np.ndarray] = {
        s: _migration_matrix(L_species[s], dt) for s in species
    }

    moves = {s: not np.allclose(L_species[s], 0.0) for s in species}

    if check_invariants:
        _check_invariants(L_species, P_mig)

    # Report configuration BEFORE integrating
    if verbose:
        print_simulation_setup(
            num_patches=num_patches, modes=modes, rn_is_list=rn_is_list,
            species=species, reactions=reactions, present=present,
            patch_species=patch_species, patch_reactions=patch_reactions,
            rate_per_patch=rate_per_patch, spec_per_patch=spec_per_patch,
            X0=X0, t_span=t_span, n_steps=n_steps, dt_out=dt_out,
            dt_couple=dt_couple, dt=dt, n_sub=n_sub,
            C=C, connectivity_given=connectivity_matrix is not None,
            D_dict=D_dict, D_default=D_default,
            transport_rates=transport_rates,
            L_species=L_species, P_mig=P_mig, moves=moves,
            method=method, rtol=rtol, atol=atol, seed=seed,
        )

        # Report the COMPLETE system (local + transport)
        _print_full_system(
            rns=rns,
            rate_per_patch=rate_per_patch,
            spec_per_patch=spec_per_patch,
            additional_laws=additional_laws,
            species=species,
            patch_species=patch_species,
            patch_reactions=patch_reactions,
            L_species=L_species,
            moves=moves,
            verbose=verbose,
        )

    # 11. Output arrays
    X_out = {s: np.full((n_steps, num_patches), np.nan) for s in species}
    flux_out = {r: np.full((n_steps, num_patches), np.nan) for r in reactions}

    # 12. Record state
    def record(step_idx, fluxes):
        for i, s in enumerate(species):
            for p in range(num_patches):
                if present[i, p]:
                    X_out[s][step_idx, p] = X[i, p]

        for p in range(num_patches):
            for j, reaction_name in enumerate(patch_reactions[p]):
                flux_out[reaction_name][step_idx, p] = fluxes[p][j]

    # 13. Advance one reaction step
    def advance_patch(p: int, delta_t: float, step_seed: Optional[int]):
        idx = patch_species_idx[p]
        n_rxn = len(patch_reactions[p])

        if n_rxn == 0:
            return X[idx, p].copy(), np.zeros(0)

        if delta_t <= 0:
            return X[idx, p].copy(), np.zeros(n_rxn)

        if modes[p] == "ode":
            x0_arg = [float(X[species_index[s], p]) for s in patch_species[p]]

            ts, fx = simulation(
                rns[p],
                rate=rate_per_patch[p],
                spec_vector=spec_per_patch[p],
                x0=x0_arg,
                t_span=(0.0, delta_t),
                n_steps=2,
                additional_laws=additional_laws,
                method=method,
                rtol=rtol,
                atol=atol,
                verbose=False,
            )

        else:
            x0_arg = {
                s: [float(X[species_index[s], p])]
                for s in patch_species[p]
            }

            ts, fx = gillespie(
                rns[p],
                rate=rate_per_patch[p],
                spec_vector=spec_per_patch[p],
                x0=x0_arg,
                t_span=(0.0, delta_t),
                n_steps=2,
                additional_laws=additional_laws,
                seed=step_seed,
                verbose=False,
            )

        x_new = ts[patch_species[p]].iloc[-1].values.astype(float)
        flux = fx[patch_reactions[p]].iloc[-1].values.astype(float)

        return x_new, flux

    # 14. Dispersal
    def disperse():
        for i, s in enumerate(species):

            if not moves[s]:
                continue

            P = P_mig[s]
            old = X[i].copy()
            new = np.zeros(num_patches, dtype=float)

            for p in range(num_patches):

                if not present[i, p] or old[p] <= 0:
                    continue

                if modes[p] == "ssa":
                    count = int(round(old[p]))
                    if count > 0:
                        new += rng.multinomial(count, P[p])
                else:
                    new += old[p] * P[p]

            for p in range(num_patches):

                if not present[i, p]:
                    new[p] = 0.0
                    continue

                if modes[p] == "ssa":
                    new[p] = _stochastic_round(np.array([new[p]]), rng)[0]

            X[i] = new

    # 15. Initial reaction flux
    fluxes = [np.zeros(len(patch_reactions[p])) for p in range(num_patches)]

    probe_dt = min(dt, 1e-9) if dt > 0 else 1e-9

    for p in range(num_patches):

        if len(patch_reactions[p]) == 0:
            continue

        try:
            _, fluxes[p] = advance_patch(p, probe_dt, None)
        except Exception as exc:                                # noqa: BLE001
            if verbose:
                print(
                    f"Warning: initial flux for patch {p} could not be "
                    f"evaluated ({type(exc).__name__}); recorded as zero."
                )

    record(0, fluxes)

    # 16. Main Lie-Trotter loop
    for step in range(1, n_steps):

        for k in range(n_sub):

            for p in range(num_patches):

                if seed is None:
                    step_seed = None
                else:
                    step_seed = int(
                        np.random.SeedSequence([seed, step, k, p])
                        .generate_state(1)[0]
                    )

                x_new, fluxes_p = advance_patch(p, dt, step_seed)

                fluxes[p] = fluxes_p
                X[patch_species_idx[p], p] = x_new

            disperse()

        record(step, fluxes)

    # 17. Summary
    if verbose:
        print_simulation_summary(species, X0, X, present)

    if return_info:
        info = {
            "species": species,
            "reactions": reactions,
            "modes": modes,
            "present": present,
            "L_species": L_species,
            "P_mig": P_mig,
            "dt": dt,
            "n_sub": n_sub,
            "num_patches": num_patches,
            "connectivity": C,
            "x0": X0,
        }
        return t_grid, X_out, flux_out, info

    return t_grid, X_out, flux_out


__all__ = [
    "simulate_metapopulation_dynamics",
    "uniform_connectivity",
    "ring_connectivity",
    "chain_connectivity",
    "total_mass",
    "print_simulation_setup",
    "print_simulation_summary",
]
















# """
# Metapopulation Simulation Module (hybrid, heterogeneous)
# ========================================================

# Discrete patches coupled by dispersal, where each patch may carry

#   - the SAME reaction network, or a DIFFERENT one per patch, and
#   - a DIFFERENT simulation mode: deterministic ('ode', via
#     `pyCOT.simulations.ode.simulation`) or stochastic ('ssa', via
#     `pyCOT.simulations.stochastic.gillespie`).

# Both existing simulators are reused as-is per patch, rather than
# reimplementing kinetics: whatever `simulation` and `gillespie` support
# (kinetic laws, spec_vector conventions, exact combinatorial propensities)
# is supported here.

# Numerical scheme
# ----------------
# Lie-Trotter operator splitting on a fixed coupling step `dt_couple`:

#     for each window [t, t+dt]:
#         1. REACTION  each patch advances independently with its own method
#         2. DISPERSAL the migration operator is applied across patches

# The splitting error is O(dt_couple), so `dt_couple` must be small compared
# with both the fastest reaction timescale and 1/D. It is NOT an exact
# simulation of the coupled process; there is no exact scheme for a system
# that is continuous in some patches and discrete in others.

# The dispersal sub-step, taken alone, IS exact, and implements the transport
# layer of metapoblaciones.pdf (secs. 3-4): for species s, transport is driven
# by a graph Laplacian L_s in R^(num_patches x num_patches) with

#     (L_s)_pq = D_s * C[q, p]                for p != q   (rate T_(q->p))
#     (L_s)_pp = -sum_{q != p} D_s * C[p, q]

# so that dx_s/dt|_transport = L_s @ x_s, matching the PDF's T_(q->p) = rate at
# which species s moves FROM patch q TO patch p. This module builds the
# transpose M = L_s^T instead (M[p, q] = D_s * C[p, q] for p != q), because M
# is the natural continuous-time-Markov-chain *generator* for one molecule:
# its rows sum to zero, and row p gives the hazard of a molecule at p jumping
# to every other patch (PDF sec. 6). Over a window dt the transition matrix
# P = expm(M dt) is row-stochastic, and P^T = expm(L_s dt) is exactly the
# PDF's exact linear solution for that window. Then:

#   - 'ode' source patch: x <- P^T x = expm(L_s dt) x    (exact linear solution)
#   - 'ssa' source patch: each molecule independently follows the CTMC with
#     generator M, so the counts leaving patch p are Multinomial(x_p, P[p, :]).
#     (exact -- PDF sec. 6 models a hop q->p as a unimolecular reaction with
#     propensity T_(q->p) * N_i^(q) -- and it cannot produce negative counts,
#     unlike an Euler dispersal step)

# Special case vs. the PDF's fully general L_i (sec. 3): every species here
# shares the SAME connectivity template `C`, only the scalar D_s varies per
# species (the PDF's own "L_i = d_i * L" special case, sec. 4). `present`
# still restricts each species to the patches whose network actually has it,
# which is what lets species have effectively different migration graphs.

# Mixed coupling: a fractional amount arriving from an 'ode' patch into an
# 'ssa' patch is rounded stochastically (floor plus Bernoulli on the remainder),
# which preserves the mean but is an approximation of the true hybrid coupling.

# Heterogeneous networks
# ----------------------
# Species are indexed on the UNION of all patch species. A species migrates
# only between patches whose networks actually contain it; migration into a
# patch lacking that species is blocked and the mass stays where it is. In the
# output, a species absent from a patch is reported as NaN, not as 0, to keep
# "not part of this patch's chemistry" distinct from "present at zero".
# """

# from __future__ import annotations

# from typing import Dict, List, Optional, Sequence, Tuple, Union

# import numpy as np
# from scipy.linalg import expm

# from pyCOT.simulations.ode import simulation
# from pyCOT.simulations.stochastic import gillespie
# from pyCOT.simulations.core import generate_default_parameters, validate_rate_list


# # ---------------------------------------------------------------------------
# # Helpers
# # ---------------------------------------------------------------------------

# def _as_patch_list(value, num_patches, name):
#     """Broadcast a single value to one entry per patch, or validate a list."""
#     if isinstance(value, (list, tuple)) and len(value) == num_patches \
#             and not isinstance(value, str):
#         return list(value)
#     return [value] * num_patches


# def _normalize_connectivity(C, num_patches, verbose=True):
#     if C is None:
#         C = np.eye(num_patches)
#     C = np.asarray(C, dtype=float)
#     if C.shape != (num_patches, num_patches):
#         raise ValueError(f"connectivity_matrix must be ({num_patches}, {num_patches}), got {C.shape}")
#     rows = C.sum(axis=1, keepdims=True)
#     if np.any(rows <= 0):
#         raise ValueError("Every row of connectivity_matrix must have a positive sum")
#     if not np.allclose(rows, 1.0, atol=1e-6):
#         if verbose:
#             print("Warning: normalizing connectivity_matrix so rows sum to 1")
#         C = C / rows
#     return C


# def _migration_matrix(C, D, dt, present):
#     """Row-stochastic transition matrix over dt for one species.

#     `present` is a boolean mask of the patches whose network contains the
#     species. Rows and columns of absent patches are removed from the
#     generator and reinserted as identity, so mass never migrates into a patch
#     that has no such species.

#     M below is the CTMC generator for one molecule (M[p, q] = D * C[p, q] =
#     T_(p->q) for p != q, rows sum to 0). Its transpose is the population-level
#     Laplacian of metapoblaciones.pdf sec. 3, L = M^T (L[p, q] = T_(q->p)), so
#     P = expm(M * dt) is row-stochastic and P^T = expm(L * dt) is the PDF's
#     exact linear transport solution -- see the module docstring for how the
#     'ode' and 'ssa' branches of `disperse()` each use this P.
#     """
#     n = C.shape[0]
#     P = np.eye(n)
#     idx = np.flatnonzero(present)
#     if len(idx) < 2 or D <= 0:
#         return P                                   # nowhere to go
#     Csub = C[np.ix_(idx, idx)]
#     # Renormalize: the mass that would have gone to absent patches stays put.
#     np.fill_diagonal(Csub, Csub.diagonal() + (C[idx].sum(axis=1) - Csub.sum(axis=1)))
#     M = D * (Csub - np.diag(Csub.sum(axis=1)))     # generator, rows sum to 0
#     P[np.ix_(idx, idx)] = expm(M * dt)
#     return P


# def _stochastic_round(x, rng):
#     """Round to integers preserving the expected value."""
#     floor = np.floor(x)
#     return floor + (rng.random(x.shape) < (x - floor))


# # ---------------------------------------------------------------------------
# # Main entry point
# # ---------------------------------------------------------------------------

# def simulate_metapopulation_dynamics(
#         rn: Union[object, Sequence],
#         rate: Union[str, List] = 'mak',
#         modes: Union[str, List[str]] = 'ode',
#         num_patches: Optional[int] = None,
#         D_dict: Optional[Dict[str, float]] = None,
#         x0_dict: Optional[Dict[str, np.ndarray]] = None,
#         spec_vector: Optional[Union[List, Sequence[List]]] = None,
#         t_span: Tuple[float, float] = (0, 20),
#         n_steps: int = 500,
#         dt_couple: Optional[float] = None,
#         connectivity_matrix: Optional[np.ndarray] = None,
#         additional_laws: Optional[Dict] = None,
#         method: str = 'LSODA',
#         rtol: float = 1e-8,
#         atol: float = 1e-10,
#         seed: Optional[int] = None,
#         verbose: bool = True):
#     """
#     Simulate metapopulation dynamics with per-patch networks and per-patch
#     simulation modes.

#     Parameters
#     ----------
#     rn : ReactionNetwork or list of ReactionNetwork
#         A single network used in every patch, or one network per patch.
#         If a list is given and `num_patches` is None, its length sets the
#         number of patches.
#     rate : str or list
#         Kinetic law(s). A single string applies everywhere. A list of strings
#         of length num_patches gives one law per patch; a list of lists gives
#         one law per reaction per patch.
#     modes : {'ode', 'ssa'} or list
#         Simulation method per patch. A single value applies to all patches.
#         'ode' calls `simulation`, 'ssa' calls `gillespie`.
#     num_patches : int, optional
#         Number of patches. Inferred from `rn`, `connectivity_matrix` or
#         `modes` when possible; otherwise defaults to 3.
#     D_dict : dict, optional
#         Dispersal rate per species {species_name: D}. Missing species get 0
#         (no dispersal), which is a deliberate default: silently inventing a
#         dispersal rate would change the model.
#     x0_dict : dict, optional
#         {species_name: array of length num_patches}. Missing species start at
#         0. If None, random values are drawn (integers 1..10 if any patch is
#         stochastic, uniform [0, 2] otherwise).
#     spec_vector : list, optional
#         Kinetic parameters. Either one spec_vector shared by all patches, or
#         a list of length num_patches with one per patch (required when the
#         patches carry different networks).
#     t_span, n_steps : tuple, int
#         Output time grid.
#     dt_couple : float, optional
#         Splitting step. Defaults to the output grid spacing. Smaller values
#         reduce the O(dt) splitting error at proportional cost.
#     connectivity_matrix : array, optional
#         (num_patches x num_patches), row-normalized. C[i, j] is the share of
#         the dispersing flux of patch i directed to patch j. Default: identity
#         (isolated patches).
#     seed : int, optional
#         Seed for the stochastic patches, the dispersal sampling and the
#         random defaults.
#     verbose : bool
#         Print the configuration summary.

#     Returns
#     -------
#     t : array (n_steps,)
#     X_out : dict {species: array (n_steps, num_patches)}
#         NaN where a species does not belong to that patch's network.
#     flux_out : dict {reaction: array (n_steps, num_patches)}
#         NaN where a reaction does not belong to that patch's network.
#         For 'ssa' patches this is the rate law evaluated at the current
#         state, which is the propensity only under mass action.
#     """
#     rng = np.random.default_rng(seed)

#     # --- Patches and networks ---------------------------------------------
#     rn_list_given = isinstance(rn, (list, tuple))
#     if num_patches is None:
#         if rn_list_given:
#             num_patches = len(rn)
#         elif connectivity_matrix is not None:
#             num_patches = np.asarray(connectivity_matrix).shape[0]
#         elif isinstance(modes, (list, tuple)):
#             num_patches = len(modes)
#         else:
#             num_patches = 3

#     rns = list(rn) if rn_list_given else [rn] * num_patches
#     if len(rns) != num_patches:
#         raise ValueError(f"Got {len(rns)} networks for {num_patches} patches")

#     modes = _as_patch_list(modes, num_patches, 'modes')
#     for m in modes:
#         if m not in ('ode', 'ssa'):
#             raise ValueError(f"mode must be 'ode' or 'ssa', got '{m}'")

#     patch_species = [[s.name for s in r.species()] for r in rns]
#     patch_reactions = [[x.name() for x in r.reactions()] for r in rns]

#     # Union, keeping first-appearance order
#     species: List[str] = []
#     for sl in patch_species:
#         for s in sl:
#             if s not in species:
#                 species.append(s)
#     reactions: List[str] = []
#     for rl in patch_reactions:
#         for x in rl:
#             if x not in reactions:
#                 reactions.append(x)

#     present = np.zeros((len(species), num_patches), dtype=bool)
#     for p, sl in enumerate(patch_species):
#         for s in sl:
#             present[species.index(s), p] = True

#     # --- Kinetic laws and parameters per patch ------------------------------
#     rate_per_patch = _as_patch_list(rate, num_patches, 'rate')
#     rate_per_patch = [validate_rate_list(rt, len(patch_reactions[p]))
#                       for p, rt in enumerate(rate_per_patch)]

#     if spec_vector is None:
#         spec_per_patch = [generate_default_parameters(
#             rate_per_patch[p], len(patch_reactions[p]), additional_laws)
#             for p in range(num_patches)]
#     elif len(spec_vector) == num_patches and spec_vector \
#             and isinstance(spec_vector[0], (list, tuple)) \
#             and spec_vector[0] and isinstance(spec_vector[0][0], (list, tuple)):
#         spec_per_patch = [list(sv) for sv in spec_vector]     # one per patch
#     else:
#         spec_per_patch = [list(spec_vector)] * num_patches    # shared
#     for p in range(num_patches):
#         if len(spec_per_patch[p]) != len(patch_reactions[p]):
#             raise ValueError(
#                 f"Patch {p}: spec_vector has {len(spec_per_patch[p])} entries "
#                 f"for {len(patch_reactions[p])} reactions. With heterogeneous "
#                 f"networks you must pass one spec_vector per patch.")

#     any_stochastic = any(m == 'ssa' for m in modes)

#     # --- Dispersal ----------------------------------------------------------
#     C = _normalize_connectivity(connectivity_matrix, num_patches, verbose)
#     if D_dict is None:
#         D_dict = {}
#     D_vec = np.array([float(D_dict.get(s, 0.0)) for s in species])

#     # --- Initial conditions -------------------------------------------------
#     X = np.zeros((len(species), num_patches))
#     if x0_dict is None:
#         for i, s in enumerate(species):
#             X[i] = (rng.integers(1, 11, num_patches) if any_stochastic
#                     else np.round(rng.uniform(0, 2.0, num_patches), 2))
#     else:
#         for s, vals in x0_dict.items():
#             if s not in species:
#                 raise ValueError(f"Species '{s}' is in no patch network (species: {species})")
#             vals = np.atleast_1d(np.asarray(vals, dtype=float))
#             if vals.size == 1:
#                 vals = np.repeat(vals, num_patches)
#             if vals.size != num_patches:
#                 raise ValueError(f"x0 for '{s}' must have length {num_patches}")
#             X[species.index(s)] = vals
#     X[~present] = 0.0
#     for p, m in enumerate(modes):
#         if m == 'ssa':
#             X[:, p] = np.round(X[:, p])
#     if verbose:
#         print("Initial conditions X0:\n", dict(zip(species, X.T.tolist())))
#     # --- Time grid ----------------------------------------------------------
#     t_grid = np.linspace(t_span[0], t_span[1], n_steps)
#     dt_out = (t_span[1] - t_span[0]) / max(n_steps - 1, 1)
#     if dt_couple is None:
#         dt_couple = dt_out
#     sub = max(1, int(round(dt_out / dt_couple)))
#     dt = dt_out / sub

#     if verbose:
#         print(f"\nPatches: {num_patches}   modes: {modes}")
#         print(f"Networks: {'heterogeneous' if rn_list_given else 'shared'}")
#         print(f"Species (union): {species}")
#         print(f"Dispersal D: {dict(zip(species, D_vec))}")
#         print(f"Connectivity matrix:\n{C}\nRow sums: {C.sum(axis=1)}")
#         print(f"Coupling step dt = {dt:g} ({sub} sub-step(s) per output step)")

#     # Migration transition matrix per species, precomputed for the fixed dt
#     P_mig = [_migration_matrix(C, D_vec[i], dt, present[i]) for i in range(len(species))]

#     # --- Output buffers -----------------------------------------------------
#     X_out = {s: np.full((n_steps, num_patches), np.nan) for s in species}
#     flux_out = {r: np.full((n_steps, num_patches), np.nan) for r in reactions}

#     def record(step_idx, fluxes):
#         for i, s in enumerate(species):
#             for p in range(num_patches):
#                 if present[i, p]:
#                     X_out[s][step_idx, p] = X[i, p]
#         for p in range(num_patches):
#             for j, r in enumerate(patch_reactions[p]):
#                 flux_out[r][step_idx, p] = fluxes[p][j]

#     # --- Per-patch reaction sub-step ---------------------------------------
#     def advance_patch(p, dt, step_seed):
#         """Advance patch p over dt with its own method. Returns (x_new, flux)."""
#         sp_p = patch_species[p]
#         x_local = [X[species.index(s), p] for s in sp_p]

#         if modes[p] == 'ode':
#             ts, fx = simulation(rns[p], rate=rate_per_patch[p],
#                                 spec_vector=spec_per_patch[p], x0=x_local,
#                                 t_span=(0.0, dt), n_steps=2,
#                                 additional_laws=additional_laws, method=method,
#                                 rtol=rtol, atol=atol, verbose=False)
#         else:
#             ts, fx = gillespie(rns[p], rate=rate_per_patch[p],
#                                spec_vector=spec_per_patch[p], x0=x_local,
#                                t_span=(0.0, dt), n_steps=2,
#                                additional_laws=additional_laws,
#                                seed=step_seed, verbose=False)

#         x_new = ts[sp_p].iloc[-1].values.astype(float)
#         flux = fx[patch_reactions[p]].iloc[-1].values.astype(float)
#         return x_new, flux

#     # --- Dispersal sub-step -------------------------------------------------
#     def disperse():
#         for i, s in enumerate(species):
#             if D_vec[i] <= 0:
#                 continue
#             P = P_mig[i]
#             row = X[i]
#             new = np.zeros(num_patches)
#             for p in range(num_patches):
#                 if row[p] <= 0:
#                     continue
#                 if modes[p] == 'ssa':
#                     # Each molecule migrates independently: exact, non-negative
#                     new += rng.multinomial(int(round(row[p])), P[p])
#                 else:
#                     new += row[p] * P[p]
#             # Stochastic rounding of fractional inflow into discrete patches
#             for p in range(num_patches):
#                 if modes[p] == 'ssa' and abs(new[p] - round(new[p])) > 1e-12:
#                     new[p] = _stochastic_round(np.array([new[p]]), rng)[0]
#             X[i] = new
#         X[~present] = 0.0

#     # --- Main loop ----------------------------------------------------------
#     fluxes = [np.zeros(len(patch_reactions[p])) for p in range(num_patches)]
#     for p in range(num_patches):
#         _, fluxes[p] = advance_patch(p, 1e-12, None if seed is None else seed)
#     record(0, fluxes)

#     for step in range(1, n_steps):
#         for k in range(sub):
#             for p in range(num_patches):
#                 sd = None if seed is None else int(seed + 1_000_003 * (step * sub + k) + p)
#                 x_new, fluxes[p] = advance_patch(p, dt, sd)
#                 for j, s in enumerate(patch_species[p]):
#                     X[species.index(s), p] = x_new[j]
#             disperse()
#         record(step, fluxes)

#     return t_grid, X_out, flux_out


# __all__ = ['simulate_metapopulation_dynamics']


















# """
# Metapopulation Simulation Module

# Implements discrete patch dynamics with local reactions and inter-patch dispersal.

# Features:
# - Multiple discrete patches/populations
# - Local reaction dynamics in each patch
# - Dispersal between patches via connectivity matrix
# - Flexible connectivity (migration networks, metapopulations)

# Mathematical form:
# dx_i^(p)/dt = f_i(x^(p)) + Σ_q [D_i * C_qp * x_i^(q) - D_i * C_pq * x_i^(p)]

# where:
# - p, q are patch indices
# - f_i(x^(p)) are local reactions in patch p
# - D_i is dispersal rate for species i
# - C_pq is connectivity from patch p to q
# """

# import numpy as np
# from scipy.integrate import solve_ivp
# from typing import List, Dict, Optional, Tuple
# import sys
# import os

# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# from kinetics import KINETIC_REGISTRY

# from .core import (
#     validate_rate_list,
#     generate_default_parameters,
#     generate_random_vector
# )


# def get_reaction_components(reaction, species):
#     """
#     Extract reactants, products, and stoichiometry from pyCOT Reaction object.
    
#     Parameters:
#     -----------
#     reaction : Reaction
#         pyCOT reaction object
#     species : list
#         List of species names
        
#     Returns:
#     --------
#     reactants : list of tuple
#         [(species_idx, coefficient), ...]
#     products : list of tuple
#         [(species_idx, coefficient), ...]
#     stoichiometry : list
#         Net stoichiometric coefficients for each species
#     """
#     try:
#         reactants = [(species.index(edge.species_name), edge.coefficient) 
#                      for edge in reaction.support_edges()]
#         products = [(species.index(edge.species_name), edge.coefficient) 
#                     for edge in reaction.products_edges()]
        
#         stoichiometry = [0] * len(species)
#         for sp_idx, coeff in reactants:
#             stoichiometry[sp_idx] -= coeff
#         for sp_idx, coeff in products:
#             stoichiometry[sp_idx] += coeff
        
#         return reactants, products, stoichiometry
#     except Exception as e:
#         print(f"Error extracting reaction components for {reaction.name()}: {e}")
#         raise


# def simulate_metapopulation_dynamics(rn,
#                                      rate: str = 'mak',
#                                      num_patches: Optional[int] = None,
#                                      D_dict: Optional[Dict] = None,
#                                      x0_dict: Optional[Dict] = None,
#                                      spec_vector: Optional[List] = None,
#                                      t_span: Tuple[float, float] = (0, 20),
#                                      n_steps: int = 500,
#                                      connectivity_matrix: Optional[np.ndarray] = None,
#                                      additional_laws: Optional[Dict] = None,
#                                      method: str = 'LSODA',
#                                      rtol: float = 1e-8,
#                                      atol: float = 1e-10):
#     """
#     Simulate metapopulation dynamics with local reactions and dispersal.
    
#     Parameters:
#     -----------
#     rn : ReactionNetwork
#         pyCOT reaction network object
#     rate : str or list
#         Kinetic law name(s)
#     num_patches : int, optional
#         Number of discrete patches/populations. Default: 3
#     D_dict : dict, optional
#         Dispersal rates {species_name: D_value}
#         Default: random uniform [0.01, 0.2]
#     x0_dict : dict, optional
#         Initial patch concentrations {species_name: 1D_array(num_patches)}
#         Default: random uniform [0, 2.0]
#     spec_vector : list, optional
#         Kinetic parameters per reaction
#     t_span : tuple
#         Time interval (t_start, t_end)
#     n_steps : int
#         Number of time points
#     connectivity_matrix : array, optional
#         Connectivity matrix (num_patches × num_patches)
#         C[i,j] = probability/rate of dispersal from patch i to j
#         Rows must sum to 1 (conservation of mass)
#         Default: random normalized matrix
#     additional_laws : dict, optional
#         Custom kinetic functions
#     method : str
#         Integration method ('LSODA', 'RK45', etc.)
#     rtol, atol : float
#         Relative and absolute tolerances
        
#     Returns:
#     --------
#     t : array
#         Time points (n_steps,)
#     X_out : dict
#         Patch concentration time series {species: array(n_steps, num_patches)}
#     flux_out : dict
#         Patch flux time series {reaction: array(n_steps, num_patches)}
        
#     Notes:
#     ------
#     - Each patch evolves with local reaction dynamics
#     - Dispersal couples patches: species move according to connectivity
#     - Connectivity matrix normalized so rows sum to 1
#     - Dispersal flux: inflow from patch j to i minus outflow from i to j
    
#     Example:
#     --------
#     # 3 patches with asymmetric connectivity
#     connectivity = np.array([
#         [0.7, 0.2, 0.1],  # Patch 0: 70% stay, 20% to 1, 10% to 2
#         [0.1, 0.8, 0.1],  # Patch 1: 10% to 0, 80% stay, 10% to 2
#         [0.2, 0.2, 0.6]   # Patch 2: 20% to 0, 20% to 1, 60% stay
#     ])
    
#     ts, X, flux = simulate_metapopulation_dynamics(
#         rn, num_patches=3, connectivity_matrix=connectivity
#     )
#     """
#     np.random.seed(seed=42)
#     species = [specie.name for specie in rn.species()]
#     reactions = [reaction.name() for reaction in rn.reactions()]
    
#     rate = validate_rate_list(rate, len(reactions))
    
#     # Default number of patches
#     if num_patches is None:
#         num_patches = 3
    
#     # Default dispersal rates
#     if D_dict is None:
#         D_dict = {sp: np.round(np.random.uniform(0.01, 0.2), 3) for sp in species}
    
#     # Default initial conditions
#     if x0_dict is None:
#         x0_dict = {sp: np.round(np.random.uniform(0, 2.0, size=num_patches), 2) 
#                    for sp in species}
#     else:
#         # Validate x0_dict
#         for sp in species:
#             if sp not in x0_dict:
#                 raise ValueError(f"Missing initial condition for species '{sp}'")
#             if len(x0_dict[sp]) != num_patches:
#                 raise ValueError(
#                     f"Initial condition for '{sp}' must have length {num_patches}"
#                 )
    
#     # Default connectivity matrix (random, row-normalized)
#     if connectivity_matrix is None:
#         connectivity_matrix = np.random.uniform(0, 1, size=(num_patches, num_patches))
#         connectivity_matrix = np.round(connectivity_matrix, 2)
        
#         # Normalize rows to sum to 1
#         row_sums = np.sum(connectivity_matrix, axis=1, keepdims=True)
#         connectivity_matrix = connectivity_matrix / row_sums
#         connectivity_matrix = np.round(connectivity_matrix, 2)
        
#         # Fix rounding errors
#         for i in range(num_patches):
#             current_sum = np.sum(connectivity_matrix[i])
#             if current_sum != 1.0:
#                 diff = 1.0 - current_sum
#                 max_idx = np.argmax(connectivity_matrix[i])
#                 connectivity_matrix[i, max_idx] += diff
#                 connectivity_matrix[i, max_idx] = np.round(
#                     connectivity_matrix[i, max_idx], 2
#                 )
#     else:
#         # Validate connectivity matrix
#         if connectivity_matrix.shape != (num_patches, num_patches):
#             raise ValueError(
#                 f"Connectivity matrix must have shape ({num_patches}, {num_patches})"
#             )
#         if not np.allclose(np.sum(connectivity_matrix, axis=1), 1.0, atol=1e-6):
#             print("Warning: Normalizing connectivity matrix so rows sum to 1")
#             connectivity_matrix = connectivity_matrix / np.sum(
#                 connectivity_matrix, axis=1, keepdims=True
#             )
#             connectivity_matrix = np.round(connectivity_matrix, 2)
    
#     print("\nConnectivity matrix:\n", connectivity_matrix)
#     print(f"Row sums: {np.sum(connectivity_matrix, axis=1)}")
    
#     # Generate kinetic parameters
#     if spec_vector is None:
#         spec_vector = generate_default_parameters(rate, len(reactions), additional_laws)
    
#     # State flattening/reshaping utilities
#     def flatten_state(x_dict):
#         return np.concatenate([x_dict[sp] for sp in species])
    
#     def reshape_state(x):
#         return {sp: x[i*num_patches:(i+1)*num_patches] 
#                 for i, sp in enumerate(species)}
    
#     x0 = flatten_state(x0_dict)
    
#     # Build kinetic law registry
#     rate_laws = dict(KINETIC_REGISTRY)
#     if additional_laws:
#         rate_laws.update(additional_laws)
    
#     # Local reaction dynamics in each patch
#     def reaction_dynamics(Xdict, rate, spec_vector):
#         dxdt_dict = {sp: np.zeros(num_patches) for sp in species}
#         flux_dict = {r: np.zeros(num_patches) for r in reactions}
        
#         # Iterate over patches
#         for patch_idx in range(num_patches):
#             local_x = [Xdict[sp][patch_idx] for sp in species]
            
#             # Compute reaction rates
#             for r_idx, reaction in enumerate(rn.reactions()):
#                 reactants, products, stoichiometry = get_reaction_components(
#                     reaction, species
#                 )
                
#                 # Compute rate based on kinetic law
#                 kinetic = rate[r_idx]
#                 params = spec_vector[r_idx]
                
#                 if kinetic == 'mak':
#                     k = params[0]
#                     v_r = k
#                     if reactants:
#                         for sp_idx, stoich in reactants:
#                             v_r *= local_x[sp_idx] ** abs(stoich)
#                 elif kinetic == 'mmk':
#                     Vmax, Km = params
#                     sp_idx = reactants[0][0] if reactants else 0
#                     v_r = (Vmax * local_x[sp_idx] / (Km + local_x[sp_idx]) 
#                            if reactants else 0)
#                 elif kinetic == 'hill':
#                     Vmax, Kd, n = params
#                     sp_idx = reactants[0][0] if reactants else 0
#                     v_r = (Vmax * (local_x[sp_idx] ** n) / 
#                            (Kd ** n + local_x[sp_idx] ** n) if reactants else 0)
#                 elif kinetic in rate_laws:
#                     # Use registry for advanced kinetics
#                     species_idx = {sp: idx for idx, sp in enumerate(species)}
#                     v_r = rate_laws[kinetic](
#                         [(species[sp_idx], coef) for sp_idx, coef in reactants],
#                         local_x,
#                         species_idx,
#                         params
#                     )
#                 else:
#                     v_r = 0
                
#                 flux_dict[reaction.name()][patch_idx] = v_r
                
#                 # Apply stoichiometry
#                 for sp_idx, stoich in enumerate(stoichiometry):
#                     if stoich != 0:
#                         dxdt_dict[species[sp_idx]][patch_idx] += stoich * v_r
        
#         return dxdt_dict, flux_dict
    
#     # Dispersal term (connectivity-based)
#     def dispersal_term(Xdict):
#         dxdt_dict = {sp: np.zeros(num_patches) for sp in species}
        
#         for sp in species:
#             D = D_dict.get(sp, 0.0)  # Dispersal rate
#             X = Xdict[sp]  # Concentration vector
            
#             for i in range(num_patches):
#                 dxdt = 0
#                 for j in range(num_patches):
#                     # Inflow from patch j to patch i
#                     dxdt += D * connectivity_matrix[j, i] * X[j]
#                     # Outflow from patch i to patch j
#                     dxdt -= D * connectivity_matrix[i, j] * X[i]
#                 dxdt_dict[sp][i] = dxdt
        
#         return dxdt_dict
    
#     # Combined ODE system
#     def combined_ode(t, x):
#         x = np.maximum(x, 0)  # Non-negativity
#         Xdict = reshape_state(x)
#         dxdt_reac, _ = reaction_dynamics(Xdict, rate, spec_vector)
#         dxdt_disperse = dispersal_term(Xdict)
#         dxdt_total = {sp: dxdt_reac[sp] + dxdt_disperse[sp] for sp in species}
#         return flatten_state(dxdt_total)
    
#     # Integrate
#     t_eval = np.linspace(t_span[0], t_span[1], n_steps)
#     sol = solve_ivp(combined_ode, t_span, x0, t_eval=t_eval,
#                     method=method, rtol=rtol, atol=atol)
    
#     # Format output: patch time series
#     X_out = {sp: np.zeros((n_steps, num_patches)) for sp in species}
#     for i, xt in enumerate(sol.y.T):
#         xt_dict = reshape_state(xt)
#         for sp in species:
#             X_out[sp][i] = xt_dict[sp]
    
#     # Compute flux time series
#     flux_out = {r: np.zeros((n_steps, num_patches)) for r in reactions}
#     for time_idx, t_val in enumerate(sol.t):
#         x_val = sol.y[:, time_idx]
#         Xdict = reshape_state(x_val)
#         _, flux_dict = reaction_dynamics(Xdict, rate, spec_vector)
        
#         for r in reactions:
#             flux_out[r][time_idx] = flux_dict[r]
    
#     return sol.t, X_out, flux_out


# __all__ = ['simulate_metapopulation_dynamics']

