# Script 9: Metapopulation minima A <-> B en dos parches
#
# Red de reacciones (data/Examples_tests/RN_A_B.txt):
#
#     R1:   A => B;
#     R2:   B => A;
#
# Si la sintaxis de tus archivos de red difiere, ajusta el archivo a la
# convencion de autopoietic.txt. Solo hacen falta las dos reacciones: las
# constantes cineticas van en spec_vector, no en el archivo.
#
# Por que esta red y no otra. Ambas reacciones son unimoleculares, de modo
# que el campo vectorial de reaccion es LINEAL y el sistema completo
# (reaccion + transporte) es un sistema lineal de dimension
# n_especies * n_parches = 4, con solucion cerrada exp(G t) z0. Eso da una
# referencia analitica exacta contra la cual medir el error del splitting
# de Lie-Trotter que usa el modulo, en lugar de suponerlo.
#
# Los dos casos que se comparan:
#
#   Caso I   A y B migran con el mismo operador   L_A = L_B
#   Caso II  solo A migra                         L_B = 0
#
# El conmutador de los dos operadores vale
#
#     [A_op, B_op] = [[ 0             , k2 (L_B - L_A) ],
#                     [ k1 (L_A - L_B), 0              ]]
#
# es decir, es proporcional a (L_A - L_B). Se anula exactamente en el caso
# I, de modo que ahi el splitting es EXACTO con cualquier dt, y no se anula
# en el caso II, donde el error decae como O(dt). La lectura: con
# reacciones lineales el error de splitting no lo produce el transporte por
# si solo, sino la combinacion de transporte DIFERENCIAL entre especies con
# reacciones que las interconvierten.

# Import necessary libraries and modules
import os
import sys
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

import numpy as np
from scipy.linalg import expm

from pyCOT.io.functions import read_txt
from pyCOT.simulations.metapopulation import (
    simulate_metapopulation_dynamics,
    uniform_connectivity,
    total_mass,
)
from pyCOT.visualization.plot_dynamics import *

##########################################################################
# Red de reacciones
##########################################################################
file_path = 'data/Examples_tests/RN_A_B.txt'
rn = read_txt(file_path)

species = [s.name for s in rn.species()]
print("Species:", species)

reactions = [reaction.name() for reaction in rn.reactions()]
print("Reactions:", reactions)

##########################################################################
# Parametros
##########################################################################
rate = 'mak'

# Constantes cineticas, EN EL ORDEN de rn.reactions() impreso arriba.
# Verifica ese orden antes de correr: si read_txt devuelve [R2, R1] hay
# que intercambiarlas.
k1 = 1.0        # R1:  A -> B
k2 = 0.5        # R2:  B -> A
spec_vector = [[k1], [k2]]

num_patches = 2

# C es estocastica por filas y su diagonal es autorretencion, que no es un
# evento de transporte. Por eso la tasa de emigracion REALIZADA no es D
# sino D (1 - C[p,p]) = 0.6 * 0.5 = 0.3 en cada sentido.
connectivity_matrix = uniform_connectivity(num_patches, p_stay=0.5)
D = 0.6

# Condiciones iniciales deliberadamente antisimetricas: todo A en el
# parche 1, todo B en el parche 2, para que el transporte tenga algo que
# hacer desde el primer instante.
x0_dict = {'A': [10.0, 0.0],
           'B': [0.0, 10.0]}

t_span = (0, 10)
n_steps = 100
dt_couple = 0.01

# Tolerancias del integrador de la capa de reaccion, mas estrictas que las
# de fabrica (rtol=1e-8, atol=1e-10). Hacen falta aqui: en el caso I el
# error de splitting es exactamente cero, de modo que con las tolerancias
# por defecto lo que se mediria seria el error del integrador y no el del
# esquema. Con estas, el splitting domina y la comparacion tiene sentido.
rtol = 1e-11
atol = 1e-13

print("\nConnectivity matrix C:")
print(connectivity_matrix)
print("Tasa de emigracion realizada D(1 - C_pp):",
      D * (1.0 - np.diag(connectivity_matrix)))

##########################################################################
# CASO I: A y B migran con el mismo operador  ->  L_A = L_B
##########################################################################
print("\n" + "=" * 70)
print("CASO I: ambas especies migran con el mismo D")
print("=" * 70)

t_I, X_I, flux_I, info_I = simulate_metapopulation_dynamics(
    rn,
    rate=rate,
    modes='ode',
    num_patches=num_patches,
    spec_vector=spec_vector,
    x0_dict=x0_dict,
    connectivity_matrix=connectivity_matrix,
    D_dict={'A': D, 'B': D},
    t_span=t_span,
    n_steps=n_steps,
    dt_couple=dt_couple,
    rtol=rtol,
    atol=atol,
    return_info=True,
    check_invariants=True,
)

plot_species_dynamics_grid(t_I, X_I, species, num_patches=num_patches,
                           mode='per_patch',
                           filename='AB_casoI_per_patch.png')
plot_species_dynamics_grid(t_I, X_I, species, num_patches=num_patches,
                           mode='per_species',
                           filename='AB_casoI_per_species.png')

##########################################################################
# CASO II: solo A migra  ->  L_B = 0
##########################################################################
print("\n" + "=" * 70)
print("CASO II: solo A migra")
print("=" * 70)

t_II, X_II, flux_II, info_II = simulate_metapopulation_dynamics(
    rn,
    rate=rate,
    modes='ode',
    num_patches=num_patches,
    spec_vector=spec_vector,
    x0_dict=x0_dict,
    connectivity_matrix=connectivity_matrix,
    D_dict={'A': D},                     # B queda con D_default = 0
    t_span=t_span,
    n_steps=n_steps,
    dt_couple=dt_couple,
    rtol=rtol,
    atol=atol,
    return_info=True,
    check_invariants=True,
)

plot_species_dynamics_grid(t_II, X_II, species, num_patches=num_patches,
                           mode='per_patch',
                           filename='AB_casoII_per_patch.png')
plot_species_dynamics_grid(t_II, X_II, species, num_patches=num_patches,
                           mode='per_species',
                           filename='AB_casoII_per_species.png')

##########################################################################
# Conservacion: el transporte no mueve la masa total
##########################################################################
print("\n" + "=" * 70)
print("CONSERVACION DE MASA")
print("=" * 70)

masa_I = total_mass(X_I)
masa_II = total_mass(X_II)

print("\n  caso    especie   t=0        t=final    variacion")
for etiqueta, masa in (("I ", masa_I), ("II", masa_II)):
    for s in species:
        print(f"  {etiqueta}      {s:<8}{masa[s][0]:<11.6f}"
              f"{masa[s][-1]:<11.6f}{masa[s][-1] - masa[s][0]:+.6f}")

# A + B sumado sobre parches es invariante: las reacciones solo
# interconvierten A y B, y el transporte no crea ni destruye nada.
for etiqueta, masa in (("I ", masa_I), ("II", masa_II)):
    total = sum(masa[s] for s in species)
    print(f"\n  caso {etiqueta}: A+B total, min = {total.min():.10f}, "
          f"max = {total.max():.10f}   (debe ser constante = 20)")

##########################################################################
# Operadores de transporte y matrices de migracion
##########################################################################
print("\n" + "=" * 70)
print("OPERADORES DE TRANSPORTE")
print("=" * 70)

for etiqueta, info in (("I", info_I), ("II", info_II)):
    print(f"\n  --- Caso {etiqueta}  (dt efectivo = {info['dt']:.6g}, "
          f"n_sub = {info['n_sub']}) ---")
    for s in species:
        L = info['L_species'][s]
        P = info['P_mig'][s]
        print(f"\n    L_{s} =")
        print("     ", str(L).replace("\n", "\n      "))
        print(f"      suma de columnas: {L.sum(axis=0)}   (debe ser 0)")
        print(f"      tasa de salida -diag(L): {-np.diag(L)}")
        print(f"    P_{s} = expm(L^T dt):")
        print("     ", str(P).replace("\n", "\n      "))
        print(f"      suma de filas: {P.sum(axis=1)}   (debe ser 1)")

##########################################################################
# Referencia analitica: el sistema completo es lineal
##########################################################################
# Apilando z = (A1, A2, B1, B2) se tiene
#
#     A_op = M (x) I_2,  con  M = [[-k1,  k2],
#                                  [ k1, -k2]]
#     B_op = diag(L_A, L_B)
#     G    = A_op + B_op
#
# y la solucion exacta es z(t) = expm(G t) z0. El modulo resuelve el mismo
# sistema por splitting, de modo que la diferencia entre ambos ES el error
# de splitting, medido y no supuesto.
##########################################################################
print("\n" + "=" * 70)
print("REFERENCIA ANALITICA  exp(G t) z0")
print("=" * 70)

M = np.array([[-k1, k2],
              [k1, -k2]])
A_op = np.kron(M, np.eye(num_patches))

z0 = np.concatenate([np.asarray(x0_dict['A'], dtype=float),
                     np.asarray(x0_dict['B'], dtype=float)])

for etiqueta, info, t_mod, X_mod in (("I", info_I, t_I, X_I),
                                     ("II", info_II, t_II, X_II)):

    L_A = info['L_species']['A']
    L_B = info['L_species']['B']

    B_op = np.zeros((4, 4))
    B_op[:2, :2] = L_A
    B_op[2:, 2:] = L_B

    G = A_op + B_op

    # Conmutador, y la formula en bloques que predice su valor.
    Com = A_op @ B_op - B_op @ A_op
    predicho = np.zeros((4, 4))
    predicho[:2, 2:] = k2 * (L_B - L_A)
    predicho[2:, :2] = k1 * (L_A - L_B)

    # Trayectoria exacta sobre la misma malla temporal del modulo.
    Z = np.array([expm(G * ti) @ z0 for ti in t_mod])

    err_A = np.abs(X_mod['A'] - Z[:, 0:2]).max()
    err_B = np.abs(X_mod['B'] - Z[:, 2:4]).max()

    print(f"\n  --- Caso {etiqueta} ---")
    print(f"    ||[A_op, B_op]||_F          = {np.linalg.norm(Com):.6e}")
    print(f"    formula k*(L_A - L_B)       = {np.linalg.norm(predicho):.6e}"
          f"   (coincide: {np.allclose(Com, predicho)})")
    print(f"    L_A = L_B                   : {np.allclose(L_A, L_B)}")
    print(f"    error maximo del modulo, A  = {err_A:.3e}")
    print(f"    error maximo del modulo, B  = {err_B:.3e}")
    print(f"    dt efectivo                 = {info['dt']:.6g}")

##########################################################################
# Orden de convergencia del splitting
##########################################################################
# Se repite el caso II reduciendo dt_couple. El error debe decaer
# proporcionalmente al paso efectivo: pendiente 1 en log-log, que es lo
# que predice Lie-Trotter. En el caso I no hay nada que converger porque
# el error ya esta en el epsilon de maquina.
#
# Los valores de dt_couple deben quedar TODOS por debajo de
# dt_out = (t1 - t0)/(n_steps - 1) ~ 0.101. Por encima de ese valor el
# paso efectivo se satura en dt_out, porque n_sub = ceil(dt_out/dt_couple)
# nunca baja de 1, y dos dt_couple distintos darian el mismo error.
#
# Costo: el ultimo valor son (n_steps-1) * n_sub * num_patches llamadas a
# solve_ivp, unas 4000. Si quieres una corrida rapida, acorta la lista.
##########################################################################
print("\n" + "=" * 70)
print("ORDEN DE CONVERGENCIA  (caso II, solo A migra)")
print("=" * 70)

L_A_II = info_II['L_species']['A']
L_B_II = info_II['L_species']['B']

B_op_II = np.zeros((4, 4))
B_op_II[:2, :2] = L_A_II
B_op_II[2:, 2:] = L_B_II

G_II = A_op + B_op_II

dt_out = (t_span[1] - t_span[0]) / (n_steps - 1)
print(f"\n  dt_out = {dt_out:.6g}   (los dt_couple deben quedar por debajo)")
print(f"\n  {'dt_couple':>12}{'n_sub':>8}{'dt efectivo':>14}{'error max':>14}"
      f"{'orden':>10}")

err_previo = None
dt_previo = None

for dtc in (0.08, 0.04, 0.02, 0.01, 0.005):

    t_c, X_c, _, info_c = simulate_metapopulation_dynamics(
        rn,
        rate=rate,
        modes='ode',
        num_patches=num_patches,
        spec_vector=spec_vector,
        x0_dict=x0_dict,
        connectivity_matrix=connectivity_matrix,
        D_dict={'A': D},
        t_span=t_span,
        n_steps=n_steps,
        dt_couple=dtc,
        rtol=rtol,
        atol=atol,
        return_info=True,
        verbose=False,
    )

    Z_c = np.array([expm(G_II * ti) @ z0 for ti in t_c])
    err = max(np.abs(X_c['A'] - Z_c[:, 0:2]).max(),
              np.abs(X_c['B'] - Z_c[:, 2:4]).max())

    dt_efectivo = info_c['dt']

    # El orden solo es calculable si el paso efectivo cambio de verdad y
    # si ambos errores estan por encima del ruido de redondeo.
    calculable = (
        err_previo is not None
        and err > 1e-12 and err_previo > 1e-12
        and not np.isclose(dt_previo, dt_efectivo)
    )

    if calculable:
        orden = np.log(err_previo / err) / np.log(dt_previo / dt_efectivo)
        print(f"  {dtc:12.4f}{info_c['n_sub']:8d}{dt_efectivo:14.6g}"
              f"{err:14.3e}{orden:10.3f}")
    else:
        print(f"  {dtc:12.4f}{info_c['n_sub']:8d}{dt_efectivo:14.6g}"
              f"{err:14.3e}{'-':>10}")

    err_previo = err
    dt_previo = dt_efectivo

print("\n  El orden debe acercarse a 1: el splitting de Lie-Trotter es de")
print("  primer orden cuando los operadores no conmutan.")

##########################################################################
# Version hibrida: parche 1 determinista, parche 2 estocastico
##########################################################################
# Con modes = ['ode', 'ssa'] la masa continua que sale del parche 1 hacia
# el 2 pasa por el redondeo estocastico del acoplamiento hibrido, que
# preserva la esperanza pero no cada realizacion. La masa total deja de
# ser exactamente constante: fluctua alrededor de 20.
##########################################################################
print("\n" + "=" * 70)
print("VERSION HIBRIDA  modes = ['ode', 'ssa']")
print("=" * 70)

t_H, X_H, flux_H = simulate_metapopulation_dynamics(
    rn,
    rate=rate,
    modes=['ode', 'ssa'],
    num_patches=num_patches,
    spec_vector=spec_vector,
    x0_dict=x0_dict,
    connectivity_matrix=connectivity_matrix,
    D_dict={'A': D, 'B': D},
    t_span=t_span,
    n_steps=n_steps,
    dt_couple=dt_couple,
    rtol=rtol,
    atol=atol,
    seed=12345,
    verbose=False,
)

plot_species_dynamics_grid(t_H, X_H, species, num_patches=num_patches,
                           mode='per_patch',
                           filename='AB_hibrido_per_patch.png')

masa_H = total_mass(X_H)
total_H = sum(masa_H[s] for s in species)

print(f"\n  A+B total: inicial = {total_H[0]:.4f}, final = {total_H[-1]:.4f}")
print(f"             min = {total_H.min():.4f}, max = {total_H.max():.4f}, "
      f"media = {total_H.mean():.4f}")
print("\n  En los casos I y II, deterministas, ese total es exactamente 20")
print("  en todo instante. Aqui fluctua: el redondeo estocastico de la")
print("  frontera ODE/SSA conserva la masa solo en esperanza. Es el unico")
print("  punto del modulo donde la conservacion deja de ser exacta.")

print("\n" + "=" * 70)
print("LECTURAS")
print("=" * 70)
print("""
  Caso I   L_A = L_B, el conmutador se anula y el splitting es EXACTO:
           el error frente a exp(G t) z0 queda en el epsilon de maquina
           con cualquier dt_couple.

  Caso II  L_A != L_B, el conmutador es proporcional a (L_A - L_B) y el
           error decae como O(dt), con orden observado cercano a 1.

  Conclusion. Con reacciones lineales, el error de splitting no lo genera
  el transporte por si solo, sino el transporte DIFERENCIAL entre especies
  combinado con reacciones que las interconvierten. Con reacciones no
  lineales, como la red autopoietica, el error aparece ademas por la no
  linealidad, de modo que ahi dt_couple importa incluso si todas las
  especies comparten el mismo operador de transporte.
""")
print("=" * 70)