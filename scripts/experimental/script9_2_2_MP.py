# Script 9: Ejemplo trabajado de metapoblaciones.pdf, seccion 5
#
# Dos parches con REDES DISTINTAS, especies (a, b, c), y transporte
# DIRIGIDO de una sola especie. Es el caso que ejercita las dos
# caracteristicas centrales del planteamiento del PDF (seccion 3):
# transporte especie-dependiente y asimetrico.
#
# ---------------------------------------------------------------------
# Archivos de red que hay que crear
# ---------------------------------------------------------------------
#
#   data/Examples_tests/RN_abc_parche1.txt
#
#       R1:	a=>b;
#       R2:	b+c=>a;
#
#   data/Examples_tests/RN_abc_parche2.txt
#
#       R1:	=>c;
#       R2:	a=>2a;
#
# Ajusta la sintaxis a la de autopoietic.txt si difiere. Las constantes
# cineticas van en spec_vector, no en el archivo.
#
# ---------------------------------------------------------------------
# El sistema, tal como lo escribe el PDF
# ---------------------------------------------------------------------
#
# Parche 1:  a -> b  (k1 xa)    y    b + c -> a  (k2 xb xc)
#
#                  (-1  +1)                  (  k1 xa  )
#     Gamma^(1) =  (+1  -1)  ,   v^(1)   =   (k2 xb xc )
#                  ( 0  -1)
#
#     a' = -k1 a + k2 bc ,   b' = k1 a - k2 bc ,   c' = -k2 bc
#
# Parche 2:  inflow  vacio -> c  (orden cero, k_in)  y  a -> 2a  (k3 xa)
#
#                  ( 0  +1)                  ( k_in  )
#     Gamma^(2) =  ( 0   0)  ,   v^(2)   =   ( k3 xa )
#                  (+1   0)
#
#     a' = k3 a ,   b' = 0 ,   c' = k_in
#
# Acoplamiento: solo la especie a migra, con tasas dirigidas
#
#                  ( -T12   T21 )
#     L_a       =  (  T12  -T21 ) ,        L_b = L_c = 0
#
# y la ecuacion de a en el parche 2 queda
#
#     a'^(2) = k3 a^(2) + T12 a^(1) - T21 a^(2)
#
# ---------------------------------------------------------------------
# Sobre la especie b en el parche 2
# ---------------------------------------------------------------------
#
# El PDF declara b como especie del parche 2 con b' = 0 (fila nula en
# Gamma^(2)). En un archivo de red no se puede declarar una especie que
# no participa en ninguna reaccion, de modo que aqui b simplemente no
# forma parte de la red del parche 2 y el modulo la marca como ausente:
# present[b, 2] = False y la salida reporta NaN, no cero.
#
# Las dos versiones son dinamicamente EQUIVALENTES. En el PDF b^(2) tiene
# derivada nula y no migra (L_b = 0), luego es constante y no influye
# sobre nada. Aqui simplemente no se sigue. La diferencia es de
# contabilidad, no de dinamica.

# Import necessary libraries and modules
import os
import sys
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

import numpy as np
from scipy.linalg import null_space

from pyCOT.io.functions import read_txt
from pyCOT.simulations.metapopulation import (
    simulate_metapopulation_dynamics,
    total_mass,
)
from pyCOT.visualization.plot_dynamics import *

##########################################################################
# Redes de reaccion, una por parche
##########################################################################
file_parche1 = 'data/Examples_tests/RN_abc_parche1.txt'
file_parche2 = 'data/Examples_tests/RN_abc_parche2.txt'

rn1 = read_txt(file_parche1)
rn2 = read_txt(file_parche2)

rn = [rn1, rn2]                      # una red POR PARCHE
num_patches = 2

for p, red in enumerate(rn):
    print(f"\nParche {p + 1}:")
    print("  Species  :", [s.name for s in red.species()])
    print("  Reactions:", [r.name() for r in red.reactions()])

##########################################################################
# Matrices estequiometricas y leyes de conservacion locales
##########################################################################
# Seccion 7 del PDF: un vector en el nucleo POR LA IZQUIERDA de Gamma da
# una cantidad conservada por las reacciones de ese parche. Como el
# transporte conserva el total de cada especie por separado, la cantidad
#
#     sum_p  l^T x^(p)
#
# es conservada globalmente si y solo si l esta en el nucleo izquierdo de
# TODOS los Gamma^(p). Con redes heterogeneas esa interseccion suele ser
# trivial, y entonces una ley valida en un parche no sobrevive al
# sistema completo.
##########################################################################
print("\n" + "=" * 70)
print("MATRICES ESTEQUIOMETRICAS Y NUCLEOS IZQUIERDOS")
print("=" * 70)

Gammas = []

for p, red in enumerate(rn):
    sm = red.stoichiometry_matrix()
    G = np.asarray(sm, dtype=float)
    Gammas.append((list(sm.species), G))

    N = null_space(G.T)

    print(f"\n  Parche {p + 1}   (filas = {list(sm.species)}, "
          f"columnas = {list(sm.reactions)})")
    print("    Gamma =")
    print("     ", str(G).replace("\n", "\n      "))
    print(f"    nucleo izquierdo, dimension {N.shape[1]}:")
    if N.shape[1] == 0:
        print("      (trivial: las reacciones no conservan ninguna "
              "combinacion lineal)")
    else:
        for j in range(N.shape[1]):
            v = N[:, j]
            v = v / np.abs(v[np.abs(v) > 1e-12]).min()      # escala legible
            print(f"      l_{j + 1} = {np.round(v, 6)}")

# Interseccion de los nucleos: conservacion GLOBAL.
# Solo tiene sentido si ambos parches comparten el mismo orden de especies;
# se alinean por nombre antes de concatenar.
nombres_union = []
for nombres, _ in Gammas:
    for s in nombres:
        if s not in nombres_union:
            nombres_union.append(s)

bloques = []
for nombres, G in Gammas:
    G_ext = np.zeros((len(nombres_union), G.shape[1]))
    for i, s in enumerate(nombres):
        G_ext[nombres_union.index(s), :] = G[i, :]
    bloques.append(G_ext)

G_todos = np.hstack(bloques)
N_global = null_space(G_todos.T)

print(f"\n  Nucleo izquierdo COMUN a los dos parches, sobre "
      f"{nombres_union}:")
print(f"    dimension {N_global.shape[1]}")
if N_global.shape[1] == 0:
    print("    NO hay ninguna cantidad conservada globalmente. Una ley de")
    print("    conservacion valida en un parche no sobrevive cuando el otro")
    print("    parche lleva otra red.")

##########################################################################
# Parametros
##########################################################################
rate = 'mak'

# Constantes cineticas, en el orden de reacciones de CADA parche impreso
# arriba. Verifica ese orden antes de correr.
k1 = 1.0        # parche 1, R1:  a -> b
k2 = 0.5        # parche 1, R2:  b + c -> a
k_in = 0.2      # parche 2, R1:  vacio -> c      (orden cero)
k3 = 0.3        # parche 2, R2:  a -> 2a

# spec_vector POR PARCHE: las redes son distintas, de modo que no puede
# compartirse uno solo. El modulo distingue este caso porque la lista
# tiene num_patches entradas y cada una calza con el numero de
# reacciones de su propio parche.
spec_vector = [[[k1], [k2]],          # parche 1
               [[k_in], [k3]]]        # parche 2

# Transporte dirigido de la especie a. La entrada [p, q] es T[p -> q].
T12 = 0.5       # a: parche 1 -> parche 2
T21 = 0.1       # a: parche 2 -> parche 1

transport_rates = {'a': np.array([[0.0, T12],
                                  [T21, 0.0]])}
# b y c no aparecen en transport_rates y D_default = 0, de modo que
# L_b = L_c = 0, como pide el PDF.

x0_dict = {'a': [1.0, 0.5],
           'b': [1.0, 0.0],      # la entrada del parche 2 se ignora: b ausente
           'c': [2.0, 0.0]}

t_span = (0, 10)
n_steps = 101
dt_couple = 0.02

rtol = 1e-11
atol = 1e-13

##########################################################################
# Simulacion determinista (seccion 4 del PDF)
##########################################################################
print("\n" + "=" * 70)
print("SIMULACION DETERMINISTA")
print("=" * 70)

t, X, flux, info = simulate_metapopulation_dynamics(
    rn,
    rate=rate,
    modes='ode',
    num_patches=num_patches,
    spec_vector=spec_vector,
    x0_dict=x0_dict,
    transport_rates=transport_rates,
    t_span=t_span,
    n_steps=n_steps,
    dt_couple=dt_couple,
    rtol=rtol,
    atol=atol,
    # verbose=False,
    return_info=True,
    check_invariants=True,
)

species = info['species']
print("\nEspecies (union):", species)
print("Matriz de presencia (filas = especies, columnas = parches):")
print(info['present'].astype(int))
print("La fila de b tiene un 0 en el parche 2: ahi la salida es NaN, no cero.")

plot_species_dynamics_grid(t, X, species, num_patches=num_patches,
                           mode='per_patch',
                           filename='PDF5_1_per_patch.png')
plot_species_dynamics_grid(t, X, species, num_patches=num_patches,
                           mode='per_species', 
                           filename='PDF5_2_per_species.png')

##########################################################################
# Estructura del operador de transporte (secciones 3 y 7)
##########################################################################
print("\n" + "=" * 70)
print("OPERADOR DE TRANSPORTE  L_a")
print("=" * 70)

L_a = info['L_species']['a']
P_a = info['P_mig']['a']

S = 0.5 * (L_a + L_a.T)          # difusion
A = 0.5 * (L_a - L_a.T)          # deriva

print("\n  L_a =")
print("   ", str(L_a).replace("\n", "\n    "))
print(f"  suma de columnas: {L_a.sum(axis=0)}   (debe ser 0: conservacion)")
print(f"  tasa de salida -diag(L_a): {-np.diag(L_a)}   "
      f"(= [T12, T21] = [{T12}, {T21}])")
print(f"  simetrica: {np.allclose(L_a, L_a.T)}")

print("\n  parte simetrica (L+L^T)/2, difusion:")
print("   ", str(S).replace("\n", "\n    "))
print("  parte antisimetrica (L-L^T)/2, deriva:")
print("   ", str(A).replace("\n", "\n    "))
print(f"  |A|max = {np.abs(A).max():.6g}   "
      f"(= |T12 - T21|/2 = {abs(T12 - T21) / 2:.6g})")

ev = np.linalg.eigvals(L_a)
print(f"\n  autovalores de L_a: {np.round(ev, 6)}")
print(f"  det = {np.linalg.det(L_a):.3e},  traza = {np.trace(L_a):.6g}")
print("  Con m = 2 el determinante es identicamente cero y el espectro es")
print("  SIEMPRE real, por asimetrico que sea el transporte. La parte")
print("  imaginaria, que es circulacion de masa, necesita un ciclo de")
print("  longitud 3 o mas. Con dos parches la asimetria no produce")
print("  circulacion: produce una distribucion estacionaria sesgada.")

# Autovector nulo de L_a: como reparte el transporte la masa si se lo
# deja solo, sin reacciones.
w, V = np.linalg.eig(L_a)
v = np.real(V[:, np.argmin(np.abs(w))])
v = v / v.sum()
print(f"\n  distribucion estacionaria del transporte solo: {np.round(v, 6)}")
print(f"  razon a2/a1 = {v[1] / v[0]:.6g}   (= T12/T21 = {T12 / T21:.6g})")
print("  El transporte simetrico daria reparto uniforme. La asimetria")
print("  acumula masa en el parche hacia el que apunta la tasa mayor.")

print("\n  P_a = expm(L_a^T dt) con dt = %.6g:" % info['dt'])
print("   ", str(P_a).replace("\n", "\n    "))
print(f"  suma de filas: {P_a.sum(axis=1)}   (debe ser 1)")

##########################################################################
# Comprobaciones exactas sobre la solucion
##########################################################################
print("\n" + "=" * 70)
print("COMPROBACIONES EXACTAS")
print("=" * 70)

# (1) En el parche 2, c solo recibe el inflow de orden cero y no migra,
#     de modo que c^(2)(t) = c^(2)(0) + k_in t exactamente.
c2_exacto = x0_dict['c'][1] + k_in * t
err_c2 = np.abs(X['c'][:, 1] - c2_exacto).max()
print("\n  (1) c^(2)(t) = c^(2)(0) + k_in t")
print(f"      error maximo = {err_c2:.3e}   "
      f"(c^(2)(T) = {X['c'][-1, 1]:.6f}, exacto = {c2_exacto[-1]:.6f})")

# (2) El transporte conserva el total de a. Todo cambio del total viene de
#     las reacciones: -k1 a1 + k2 b1 c1 en el parche 1, y k3 a2 en el 2.
a1, a2 = X['a'][:, 0], X['a'][:, 1]
b1, c1 = X['b'][:, 0], X['c'][:, 0]

dtot_reacciones = -k1 * a1 + k2 * b1 * c1 + k3 * a2
total_a = total_mass(X)['a']

integral = np.concatenate([[0.0], np.cumsum(
    0.5 * (dtot_reacciones[1:] + dtot_reacciones[:-1]) * np.diff(t))])
err_a = np.abs((total_a - total_a[0]) - integral).max()

print("\n  (2) d/dt sum_p a^(p) = solo reacciones  (el transporte no aporta)")
print(f"      total a: inicial = {total_a[0]:.6f}, final = {total_a[-1]:.6f}")
print("      error maximo frente a la integral de las reacciones = "
      f"{err_a:.3e}")
print("      (limitado por el trapecio sobre la malla de salida y por el")
print("       sesgo O(dt) del splitting, no por el modulo)")

# (3) En el parche 1, a + b es conservada por las reacciones (nucleo
#     izquierdo l = (1,1,0)), pero el transporte de a la rompe: su
#     variacion es exactamente el balance de migracion.
ab1 = a1 + b1
balance = -T12 * a1 + T21 * a2
integral_ab = np.concatenate([[0.0], np.cumsum(
    0.5 * (balance[1:] + balance[:-1]) * np.diff(t))])
err_ab = np.abs((ab1 - ab1[0]) - integral_ab).max()

print("\n  (3) d/dt (a+b)^(1) = -T12 a^(1) + T21 a^(2)")
print("      una ley de conservacion LOCAL rota exactamente por el "
      "transporte")
print(f"      error maximo = {err_ab:.3e}")

##########################################################################
# Difusion frente a adveccion (seccion 7)
##########################################################################
# Se repite el escenario con transporte SIMETRICO de la misma intensidad
# total, T12 = T21 = (0.5 + 0.1)/2 = 0.3. Solo cambia la direccion, no la
# magnitud, de modo que cualquier diferencia proviene de la asimetria.
##########################################################################
print("\n" + "=" * 70)
print("DIFUSION FRENTE A ADVECCION")
print("=" * 70)

T_sim = 0.5 * (T12 + T21)

t_s, X_s, _, info_s = simulate_metapopulation_dynamics(
    rn,
    rate=rate,
    modes='ode',
    num_patches=num_patches,
    spec_vector=spec_vector,
    x0_dict=x0_dict,
    transport_rates={'a': np.array([[0.0, T_sim],
                                    [T_sim, 0.0]])},
    t_span=t_span,
    n_steps=n_steps,
    dt_couple=dt_couple,
    rtol=rtol,
    atol=atol,
    return_info=True,
    verbose=False,
)

plot_species_dynamics_grid(t_s, X_s, species, num_patches=num_patches,
                           mode='per_patch',
                           filename='PDF5_3_simetrico_per_patch.png')

L_s = info_s['L_species']['a']
print(f"\n  simetrico T12 = T21 = {T_sim}")
print(f"    autovalores: {np.round(np.linalg.eigvals(L_s), 6)}")
print(f"    simetrica: {np.allclose(L_s, L_s.T)}")
print(f"\n  dirigido   T12 = {T12}, T21 = {T21}")
print(f"    autovalores: {np.round(ev, 6)}")
print("    Mismo espectro: traza = -(T12+T21) = -0.6 en ambos casos, y")
print("    det = 0 siempre. La asimetria no cambia la velocidad de")
print("    relajacion, cambia hacia donde se relaja.")

print(f"\n  {'':>8}{'a^(1)':>12}{'a^(2)':>12}{'a2/a1':>12}")
print(f"  {'dirigido':>8}{X['a'][-1, 0]:12.5f}{X['a'][-1, 1]:12.5f}"
      f"{X['a'][-1, 1] / X['a'][-1, 0]:12.5f}")
print(f"  {'simetrico':>8}{X_s['a'][-1, 0]:12.5f}{X_s['a'][-1, 1]:12.5f}"
      f"{X_s['a'][-1, 1] / X_s['a'][-1, 0]:12.5f}")
print("\n  El transporte dirigido bombea a hacia el parche 2, donde la")
print("  autocatalisis lo amplifica. Es la co-dependencia que senala el")
print("  PDF: la produccion autocatalitica del parche 2 se alimenta del")
print("  parche 1 segun el balance de las tasas dirigidas.")

##########################################################################
# Formulacion mesoscopica (seccion 6)
##########################################################################
# En regimen de pocas copias se pasa a numeros enteros de moleculas. Dos
# clases de eventos: reacciones dentro de un parche, con propensidad
# a_j = k_j prod C(N_i, alpha_ij), y saltos de transporte, que son
# reacciones unimoleculares de propensidad T[q->p] N_i^(q).
#
# Con modes = ['ode', 'ssa'] el parche 1 es determinista y el 2 discreto.
# La masa continua que cruza la frontera pasa por el redondeo estocastico
# del acoplamiento, que preserva la esperanza pero no cada realizacion.
##########################################################################
print("\n" + "=" * 70)
print("FORMULACION MESOSCOPICA")
print("=" * 70)

OMEGA = 10.0
x0_meso = {s: [float(np.round(v * OMEGA)) for v in vals]
           for s, vals in x0_dict.items()}
print(f"\n  Condiciones iniciales escaladas por OMEGA = {OMEGA:g}: {x0_meso}")
print("  Las reacciones de orden cero y bimoleculares cambian de escala")
print("  con k -> k OMEGA^(1-orden); aqui solo se escalan las condiciones")
print("  iniciales, de modo que el regimen NO es el mismo que el")
print("  determinista de arriba. Es una ilustracion del algoritmo, no una")
print("  comparacion cuantitativa.")

for etiqueta, modos, nombre in (
        ("hibrido  ['ode','ssa']", ['ode', 'ssa'], 'PDF5_4_hibrido.png'),
        ("estocastico  'ssa'", 'ssa', 'PDF5_5_ssa.png')):

    t_m, X_m, _ = simulate_metapopulation_dynamics(
        rn,
        rate=rate,
        modes=modos,
        num_patches=num_patches,
        spec_vector=spec_vector,
        x0_dict=x0_meso,
        transport_rates=transport_rates,
        t_span=t_span,
        n_steps=n_steps,
        dt_couple=dt_couple,
        seed=12345,
        verbose=False,
    )

    plot_species_dynamics_grid(t_m, X_m, species, num_patches=num_patches,
                               mode='per_patch', filename=nombre)

    masa = total_mass(X_m)
    print(f"\n  {etiqueta}")
    for s in species:
        print(f"    {s}: inicial = {masa[s][0]:8.2f}, "
              f"final = {masa[s][-1]:8.2f}")

print("\n  En el caso hibrido el total de a puede fluctuar por el redondeo")
print("  de la frontera ODE/SSA, que conserva la masa solo en esperanza.")
print("  En el caso 'ssa' puro el transporte es un muestreo multinomial y")
print("  conserva la masa exactamente en cada realizacion.")

##########################################################################
print("\n" + "=" * 70)
print("LECTURAS")
print("=" * 70)
print("""
  Redes heterogeneas. Cada parche lleva su propio Gamma y su propio
  vector de velocidades. La ley de conservacion a+b del parche 1 no
  sobrevive globalmente, porque el parche 2 no la comparte: el nucleo
  izquierdo comun es trivial.

  Transporte especie-dependiente. Solo a migra. b y c quedan con L = 0,
  de modo que su unica dinamica es local. Aun asi, perturbar a cambia b y
  c del parche 1, porque las reacciones locales las acoplan.

  Transporte asimetrico. Con dos parches la asimetria NO produce
  circulacion: det(L_a) = 0 y el espectro es real sea cual sea la
  asimetria. Lo que produce es una distribucion estacionaria sesgada, con
  razon a2/a1 = T12/T21. Para ver autovalores complejos hacen falta al
  menos tres parches, porque la circulacion necesita un ciclo.

  Co-dependencia. El bombeo dirigido 1 -> 2 alimenta la autocatalisis del
  parche 2, que amplifica lo que recibe y devuelve parte al parche 1. Es
  justamente el mecanismo que el PDF senala al final de la seccion 5.
""")
print("=" * 70)