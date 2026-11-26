# Script 9: Metapoblacion de TRES parches con migracion de a y c
#
# Extension del ejemplo de dos parches (script9_2_2_MP.py). Tres parches
# con REDES DISTINTAS sobre las especies (a, b, c) y DOS operadores de
# transporte simultaneos y de naturaleza opuesta:
#
#     a  ->  anillo DIRIGIDO 1 -> 2 -> 3 -> 1, con tasas desiguales
#     c  ->  cadena SIMETRICA 1 <-> 2 <-> 3
#     b  ->  no migra  (L_b = 0)
#
# Tres especies, tres regimenes de transporte distintos en el mismo
# sistema. Es el planteamiento de transporte especie-dependiente de la
# seccion 3 llevado al caso donde de verdad se nota.
#
# ---------------------------------------------------------------------
# Por que TRES parches y no dos
# ---------------------------------------------------------------------
#
# El script de dos parches termina con esta observacion: con m = 2 el
# determinante de L es identicamente cero y el espectro es SIEMPRE real,
# por asimetrico que sea el transporte. La asimetria produce una
# distribucion estacionaria sesgada, no circulacion. Para ver autovalores
# complejos hacen falta al menos tres parches, porque la circulacion
# necesita un ciclo.
#
# Este script cierra ese ciclo y verifica la prediccion. Con el anillo
# dirigido 1 -> 2 -> 3 -> 1 de tasas T1, T2, T3:
#
#     traza L_a = -(T1 + T2 + T3),   det L_a = 0
#
# y los dos autovalores no nulos son raices de
#
#     lambda^2 + (T1+T2+T3) lambda + (T1 T2 + T1 T3 + T2 T3) = 0
#
# cuyo discriminante (T1+T2+T3)^2 - 4(T1T2+T1T3+T2T3) es NEGATIVO para
# las tasas elegidas aqui. Resultado: par complejo conjugado, es decir
# relajacion oscilatoria, que es la firma de la circulacion de masa.
#
# La medida limpia de esa circulacion es la corriente estacionaria. Si
# pi es el autovector nulo de L_a normalizado a suma 1, la corriente neta
# en la arista p -> q vale
#
#     J[p,q] = T[p,q] pi_p - T[q,p] pi_q
#
# Para CUALQUIER operador simetrico J es identicamente cero: eso es
# balance detallado. Para el anillo dirigido J es la MISMA constante en
# las tres aristas, es decir hay un flujo neto girando alrededor del
# ciclo sin que la distribucion cambie. En el anillo, pi_p es
# proporcional a 1/T_p, de modo que
#
#     J = 1 / (1/T1 + 1/T2 + 1/T3)
#
# Se comprueba numericamente mas abajo.
#
# ---------------------------------------------------------------------
# Archivos de red que hay que crear
# ---------------------------------------------------------------------
#
#   data/Examples_tests/RN_abc_parche1.txt      (ya existe)
#
#       R1:	a=>b;
#       R2:	b+c=>a;
#
#   data/Examples_tests/RN_abc_parche2.txt      (ya existe)
#
#       R1:	=>c;
#       R2:	a=>2a;
#
#   data/Examples_tests/RN_abc_parche3.txt      (NUEVO)
#
#       R1:	a+b=>c;
#       R2:	c=>b;
#
# Las constantes cineticas van en spec_vector, no en los archivos.
#
# ---------------------------------------------------------------------
# El sistema
# ---------------------------------------------------------------------
#
# Parche 1, CONVERSOR:   a -> b  (k1 a)   y   b + c -> a  (k2 b c)
#
#     a' = -k1 a + k2 bc ,  b' = k1 a - k2 bc ,  c' = -k2 bc
#
#     nucleo izquierdo l = (1,1,0):  a + b conservada por las reacciones
#
# Parche 2, FUENTE AUTOCATALITICA:   vacio -> c  (k_in)   y   a -> 2a  (k3 a)
#
#     a' = k3 a ,  c' = k_in                 (b no existe en este parche)
#
#     nucleo izquierdo trivial
#
# Parche 3, RECICLADOR:   a + b -> c  (k4 a b)   y   c -> b  (k5 c)
#
#     a' = -k4 ab ,  b' = -k4 ab + k5 c ,  c' = k4 ab - k5 c
#
#     nucleo izquierdo l = (0,1,1):  b + c conservada por las reacciones
#
# Cada parche conserva algo distinto, o nada. La interseccion de los tres
# nucleos izquierdos es trivial, de modo que NINGUNA de las dos leyes
# locales sobrevive al sistema completo. Se verifica abajo.
#
# ---------------------------------------------------------------------
# El ciclo espacial
# ---------------------------------------------------------------------
#
# Ningun parche cierra un ciclo por si solo. El sistema lo cierra en el
# espacio:
#
#     parche 2 produce a por autocatalisis
#          -> el anillo dirigido lleva a hacia el parche 3
#               -> el parche 3 convierte a + b en c
#                    -> la cadena simetrica difunde c hacia el parche 1
#                         -> el parche 1 usa c para regenerar a desde b
#
# Es co-dependencia espacial: la clausura del ciclo no esta en ninguna
# red local, esta en la combinacion de las tres con el transporte.
#
# ---------------------------------------------------------------------
# Una identidad exacta que sobrevive al caso de tres parches
# ---------------------------------------------------------------------
#
# En el script de dos parches la identidad exacta era
# c^(2)(t) = c^(2)(0) + k_in t, y se apoyaba en que c no migraba. Aqui c
# SI migra, de modo que esa identidad ya no vale. En su lugar se elige
#
#     k3 = T2         (autocatalisis del parche 2 = su tasa de salida)
#
# con lo cual la ecuacion de a en el parche 2 queda
#
#     a'^(2) = k3 a^(2) + T1 a^(1) - T2 a^(2) = T1 a^(1)
#
# porque los dos terminos en a^(2) se cancelan EXACTAMENTE. La
# autocatalisis del parche 2 solo alcanza para reponer lo que emigra: todo
# su crecimiento neto es lo que le llega del parche 1. Queda la identidad
#
#     a^(2)(t) = a^(2)(0) + T1 int_0^t a^(1)(s) ds
#
# que se comprueba por trapecio mas abajo.
#
# ---------------------------------------------------------------------
# Sobre la especie b en el parche 2
# ---------------------------------------------------------------------
#
# Igual que en el caso de dos parches, b no participa en ninguna reaccion
# del parche 2, de modo que no forma parte de esa red. El modulo la marca
# ausente, present[b, 2] = False, y la salida reporta NaN, no cero.

# Import necessary libraries and modules
import os
import sys
project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src_path = os.path.join(project_root, 'src')
sys.path.insert(0, src_path)

import numpy as np
from scipy.linalg import null_space, expm

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
file_parche3 = 'data/Examples_tests/RN_abc_parche3.txt'

rn1 = read_txt(file_parche1)
rn2 = read_txt(file_parche2)
rn3 = read_txt(file_parche3)

rn = [rn1, rn2, rn3]                 # una red POR PARCHE
num_patches = 3

etiquetas_parche = ['CONVERSOR', 'FUENTE AUTOCATALITICA', 'RECICLADOR']

for p, red in enumerate(rn):
    print(f"\nParche {p + 1}  ({etiquetas_parche[p]}):")
    print("  Species  :", [s.name for s in red.species()])
    print("  Reactions:", [r.name() for r in red.reactions()])

##########################################################################
# Matrices estequiometricas y leyes de conservacion locales
##########################################################################
# Un vector en el nucleo POR LA IZQUIERDA de Gamma da una cantidad
# conservada por las reacciones de ese parche. Como el transporte conserva
# el total de cada especie por separado, la cantidad
#
#     sum_p  l^T x^(p)
#
# es conservada globalmente si y solo si l esta en el nucleo izquierdo de
# TODOS los Gamma^(p). Con tres redes distintas esa interseccion es aun
# mas facil de romper que con dos.
##########################################################################
print("\n" + "=" * 70)
print("MATRICES ESTEQUIOMETRICAS Y NUCLEOS IZQUIERDOS")
print("=" * 70)

Gammas = []

for p, red in enumerate(rn):
    sm = red.stoichiometry_matrix()          # Matriz estequiometrica Gamma^(p)
    G = np.asarray(sm, dtype=float)          # Convertir a array para null_space
    Gammas.append((list(sm.species), G))

    N = null_space(G.T)                      # nucleo izquierdo de Gamma^(p)

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

# Interseccion de los nucleos: conservacion GLOBAL. Los parches no
# comparten ni el mismo conjunto ni el mismo orden de especies, de modo
# que se alinean por nombre antes de concatenar.
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

print(f"\n  Nucleo izquierdo COMUN a los tres parches, sobre "
      f"{nombres_union}:")
print(f"    dimension {N_global.shape[1]}")
if N_global.shape[1] == 0:
    print("    NO hay ninguna cantidad conservada globalmente. El parche 1")
    print("    conserva a+b y el parche 3 conserva b+c, pero son leyes")
    print("    distintas y el parche 2 no comparte ninguna: la interseccion")
    print("    de los tres nucleos es trivial.")

##########################################################################
# Parametros
##########################################################################
rate = 'mak'

# Tasas de transporte. La entrada [p, q] es T[p -> q].
#
#   a: anillo DIRIGIDO 1 -> 2 -> 3 -> 1, con tasas DESIGUALES. La
#      desigualdad sesga la distribucion estacionaria; el hecho de que sea
#      un ciclo es lo que produce la circulacion.
T1 = 0.5        # a: parche 1 -> parche 2
T2 = 0.3        # a: parche 2 -> parche 3
T3 = 0.2        # a: parche 3 -> parche 1

T_a = np.zeros((num_patches, num_patches))
T_a[0, 1] = T1
T_a[1, 2] = T2
T_a[2, 0] = T3

#   c: cadena SIMETRICA 1 <-> 2 <-> 3 (sin arista 1 <-> 3). Simetrica, de
#      modo que cumple balance detallado: corriente estacionaria nula.
Tc = 0.25

T_c = np.zeros((num_patches, num_patches))
T_c[0, 1] = T_c[1, 0] = Tc
T_c[1, 2] = T_c[2, 1] = Tc

transport_rates = {'a': T_a, 'c': T_c}
# b no aparece en transport_rates y D_default = 0, de modo que L_b = 0.

# Constantes cineticas, en el orden de reacciones de CADA parche impreso
# arriba. Verifica ese orden antes de correr.
k1 = 1.0        # parche 1, R1:  a -> b
k2 = 0.5        # parche 1, R2:  b + c -> a
k_in = 0.2      # parche 2, R1:  vacio -> c        (orden cero)
k3 = T2         # parche 2, R2:  a -> 2a           (igualada a la tasa de
                #                                   salida, ver cabecera)
k4 = 0.4        # parche 3, R1:  a + b -> c
k5 = 0.25       # parche 3, R2:  c -> b

# spec_vector POR PARCHE: las redes son distintas, de modo que no puede
# compartirse uno solo. El modulo distingue este caso porque la lista
# tiene num_patches entradas y cada una calza con el numero de reacciones
# de su propio parche.
spec_vector = [[[k1], [k2]],          # parche 1
               [[k_in], [k3]],        # parche 2
               [[k4], [k5]]]          # parche 3

# Condiciones iniciales. Toda la a arranca en el parche 1 y toda la c en
# el parche 3: el parche 2 no tiene con que autocatalizar y el parche 1 no
# tiene con que regenerar a. Ambos dependen del transporte para arrancar.
x0_dict = {'a': [3.0, 0.0, 0.0],
           'b': [1.0, 0.0, 1.0],      # la entrada del parche 2 se ignora
           'c': [0.0, 0.0, 2.0]}

t_span = (0, 20)
n_steps = 201
dt_couple = 0.02

rtol = 1e-11
atol = 1e-13

##########################################################################
# Simulacion determinista
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
                           filename='MP3_1_per_patch.png')
plot_species_dynamics_grid(t, X, species, num_patches=num_patches,
                           mode='per_species',
                           filename='MP3_2_per_species.png')

##########################################################################
# Estructura de los operadores de transporte
##########################################################################
print("\n" + "=" * 70)
print("OPERADORES DE TRANSPORTE  L_a  (anillo dirigido)  y  L_c  (cadena)")
print("=" * 70)


def reporte_operador(nombre, L, T):
    """Descomposicion, espectro, estacionario y corriente de un L."""
    S = 0.5 * (L + L.T)              # difusion
    A = 0.5 * (L - L.T)              # deriva

    print(f"\n  --- especie {nombre} ---")
    print("    L =")
    print("     ", str(np.round(L, 6)).replace("\n", "\n      "))
    print(f"    suma de columnas: {np.round(L.sum(axis=0), 12)}   "
          f"(debe ser 0: conservacion)")
    print(f"    tasa de salida -diag(L): {-np.diag(L)}")
    print(f"    simetrica: {np.allclose(L, L.T)}")
    print(f"    |parte antisimetrica|max = {np.abs(A).max():.6g}   "
          f"|parte simetrica|max = {np.abs(S).max():.6g}")

    ev = np.linalg.eigvals(L)
    orden = np.argsort(-ev.real)
    ev = ev[orden]
    print(f"    autovalores: {np.round(ev, 6)}")
    print(f"    traza = {np.trace(L):.6g},  det = {np.linalg.det(L):.3e}")
    print(f"    |Im| maxima = {np.abs(ev.imag).max():.6g}   -> "
          f"{'HAY circulacion' if np.abs(ev.imag).max() > 1e-12 else 'sin circulacion'}")

    # Distribucion estacionaria: autovector nulo normalizado a suma 1.
    w, V = np.linalg.eig(L)
    pi = np.real(V[:, np.argmin(np.abs(w))])
    pi = pi / pi.sum()
    print(f"    distribucion estacionaria pi = {np.round(pi, 6)}")

    # Corriente neta por arista en el estacionario. Cero para todo
    # operador simetrico: eso es balance detallado.
    J = np.zeros_like(T)
    for p in range(T.shape[0]):
        for q in range(T.shape[0]):
            if p != q:
                J[p, q] = T[p, q] * pi[p] - T[q, p] * pi[q]
    print("    corriente neta J[p,q] = T[p,q] pi_p - T[q,p] pi_q:")
    print("     ", str(np.round(J, 6)).replace("\n", "\n      "))
    print(f"    |J|max = {np.abs(J).max():.6g}")

    return pi, J


pi_a, J_a = reporte_operador('a', info['L_species']['a'], T_a)
pi_c, J_c = reporte_operador('c', info['L_species']['c'], T_c)

print(f"\n  --- especie b ---")
print(f"    L_b = 0, no migra. |L_b|max = "
      f"{np.abs(info['L_species']['b']).max():.6g}")

##########################################################################
# Las dos predicciones cerradas sobre el anillo dirigido
##########################################################################
print("\n" + "=" * 70)
print("PREDICCIONES CERRADAS SOBRE EL ANILLO DIRIGIDO")
print("=" * 70)

# (i) pi_p proporcional a 1/T_p
pi_teorico = (1.0 / np.array([T1, T2, T3]))
pi_teorico = pi_teorico / pi_teorico.sum()
print(f"\n  (i) pi proporcional a 1/T_p")
print(f"      teorico  = {np.round(pi_teorico, 6)}")
print(f"      numerico = {np.round(pi_a, 6)}")
print(f"      error maximo = {np.abs(pi_a - pi_teorico).max():.3e}")

# (ii) corriente identica en las tres aristas, J = 1/sum(1/T_p)
J_teorico = 1.0 / (1.0 / T1 + 1.0 / T2 + 1.0 / T3)
J_aristas = np.array([J_a[0, 1], J_a[1, 2], J_a[2, 0]])
print(f"\n  (ii) corriente constante alrededor del ciclo")
print(f"      teorica  = {J_teorico:.10f}")
print(f"      numerica = {np.round(J_aristas, 10)}")
print(f"      error maximo = {np.abs(J_aristas - J_teorico).max():.3e}")
print("      La misma corriente en las tres aristas y una distribucion que")
print("      no cambia: eso es circulacion pura. Imposible con m = 2.")

# (iii) el par complejo predicho por el polinomio caracteristico
suma = T1 + T2 + T3
producto = T1 * T2 + T1 * T3 + T2 * T3
disc = suma ** 2 - 4.0 * producto
print(f"\n  (iii) autovalores no nulos: raices de "
      f"lambda^2 + {suma:g} lambda + {producto:g}")
print(f"      discriminante = {disc:.6g}   "
      f"({'NEGATIVO -> par complejo' if disc < 0 else 'positivo -> reales'})")
lam = np.roots([1.0, suma, producto])
print(f"      teoricos  = {np.round(lam, 6)}")
ev_a = np.linalg.eigvals(info['L_species']['a'])
ev_a = ev_a[np.argsort(-np.abs(ev_a))][:2]
print(f"      numericos = {np.round(ev_a, 6)}")

##########################################################################
# Transporte puro: referencia analitica exp(L t) x0
##########################################################################
# Si todas las constantes cineticas son cero, la capa de reaccion es la
# identidad y el modulo solo aplica migracion. Como la matriz de migracion
# es P = expm(L^T dt) y los subpasos se componen, el resultado debe
# coincidir con expm(L t) x0 a precision de maquina. Es el analogo de la
# referencia analitica del script de dos parches, y valida el operador de
# transporte por separado, incluida la circulacion.
##########################################################################
print("\n" + "=" * 70)
print("TRANSPORTE PURO: referencia analitica  exp(L t) x0")
print("=" * 70)

spec_cero = [[[0.0], [0.0]] for _ in range(num_patches)]

t_T, X_T, _, info_T = simulate_metapopulation_dynamics(
    rn,
    rate=rate,
    modes='ode',
    num_patches=num_patches,
    spec_vector=spec_cero,
    x0_dict=x0_dict,
    transport_rates=transport_rates,
    t_span=t_span,
    n_steps=n_steps,
    dt_couple=dt_couple,
    rtol=rtol,
    atol=atol,
    return_info=True,
    verbose=False,
)

plot_species_dynamics_grid(t_T, X_T, species, num_patches=num_patches,
                           mode='per_species',
                           filename='MP3_3_transporte_puro_per_species.png')

print("\n  especie   error maximo frente a exp(L t) x0")
for s in ('a', 'c'):
    L = info_T['L_species'][s]
    x0_s = np.asarray(x0_dict[s], dtype=float)
    ref = np.array([expm(L * ti) @ x0_s for ti in t_T])
    err = np.abs(X_T[s] - ref).max()
    print(f"  {s:<9} {err:.3e}")

print("\n  Sin reacciones el splitting es exacto: la capa de reaccion es la")
print("  identidad y las matrices de migracion se componen exactamente,")
print("  expm(L^T dt)^n = expm(L^T n dt). El residuo que queda es del")
print("  integrador, no del esquema.")

# La firma temporal de la circulacion: bajo transporte puro, a en un
# parche NO se acerca monotonamente a su valor estacionario, sino que lo
# sobrepasa. c, simetrica, si es monotona.
print("\n  Firma temporal de la circulacion (transporte puro):")
for s in ('a', 'c'):
    total_s = np.nansum(X_T[s], axis=1)
    pi_s = pi_a if s == 'a' else pi_c
    for p in range(num_patches):
        serie = X_T[s][:, p]
        objetivo = pi_s[p] * total_s[0]
        # Un sobrepaso existe si la serie cruza su valor final y vuelve.
        exceso = serie - objetivo
        cambios = np.sum(np.abs(np.diff(np.sign(exceso[np.abs(exceso) > 1e-9]))) > 0)
        print(f"    {s}^({p + 1}): estacionario = {objetivo:8.5f}, "
              f"cruces del estacionario = {cambios}")
print("    Mas de un cruce significa oscilacion amortiguada alrededor del")
print("    estacionario, que es lo que produce el par complejo. Con un")
print("    operador simetrico cada serie cruza a lo sumo una vez.")

##########################################################################
# Comprobaciones exactas sobre la solucion con reacciones
##########################################################################
print("\n" + "=" * 70)
print("COMPROBACIONES EXACTAS  (sistema completo)")
print("=" * 70)


def integral_trapecio(f, t):
    """Integral acumulada de f sobre la malla t, por trapecio."""
    return np.concatenate([[0.0], np.cumsum(
        0.5 * (f[1:] + f[:-1]) * np.diff(t))])


a1, a2, a3 = X['a'][:, 0], X['a'][:, 1], X['a'][:, 2]
b1, b3 = X['b'][:, 0], X['b'][:, 2]
c1, c2, c3 = X['c'][:, 0], X['c'][:, 1], X['c'][:, 2]

# (1) Identidad exacta del parche 2: como k3 = T2, los terminos en a^(2)
#     se cancelan y queda a'^(2) = T1 a^(1).
a2_pred = x0_dict['a'][1] + T1 * integral_trapecio(a1, t)
err_a2 = np.abs(a2 - a2_pred).max()
print(f"\n  (1) a^(2)(t) = a^(2)(0) + T1 int a^(1)     [porque k3 = T2 = {T2}]")
print(f"      error maximo = {err_a2:.3e}   "
      f"(a^(2)(T) = {a2[-1]:.6f}, predicho = {a2_pred[-1]:.6f})")
print("      La autocatalisis del parche 2 solo repone lo que emigra: todo")
print("      su crecimiento neto proviene del parche 1.")

# (2) El transporte conserva el total de cada especie. Todo cambio del
#     total viene de las reacciones.
dtot = {
    'a': -k1 * a1 + k2 * b1 * c1 + k3 * a2 - k4 * a3 * b3,
    'b': k1 * a1 - k2 * b1 * c1 - k4 * a3 * b3 + k5 * c3,
    'c': -k2 * b1 * c1 + k_in + k4 * a3 * b3 - k5 * c3,
}

print("\n  (2) d/dt sum_p s^(p) = solo reacciones  (el transporte no aporta)")
print(f"      {'especie':<9}{'inicial':>11}{'final':>12}{'error max':>13}")
for s in species:
    total_s = total_mass(X)[s]
    err = np.abs((total_s - total_s[0]) - integral_trapecio(dtot[s], t)).max()
    print(f"      {s:<9}{total_s[0]:11.6f}{total_s[-1]:12.6f}{err:13.3e}")
print("      (limitado por el trapecio sobre la malla de salida y por el")
print("       sesgo O(dt) del splitting, no por el modulo)")

# (3) Las dos leyes locales, rotas exactamente por el transporte de la
#     especie que migra. En el parche 1, b no migra, de modo que toda la
#     variacion de a+b viene del balance de a. En el parche 3, a y c
#     migran ambas, pero b no: la variacion de b+c viene del balance de c.
L_a_num = info['L_species']['a']
L_c_num = info['L_species']['c']

ab1 = a1 + b1
balance_a1 = L_a_num[0, 0] * a1 + L_a_num[0, 1] * a2 + L_a_num[0, 2] * a3
err_ab1 = np.abs((ab1 - ab1[0]) - integral_trapecio(balance_a1, t)).max()

bc3 = b3 + c3
balance_c3 = L_c_num[2, 0] * c1 + L_c_num[2, 1] * c2 + L_c_num[2, 2] * c3
err_bc3 = np.abs((bc3 - bc3[0]) - integral_trapecio(balance_c3, t)).max()

print("\n  (3) leyes de conservacion LOCALES rotas exactamente por el "
      "transporte")
print(f"      d/dt (a+b)^(1) = sum_q L_a[1,q] a^(q)   "
      f"error maximo = {err_ab1:.3e}")
print(f"      d/dt (b+c)^(3) = sum_q L_c[3,q] c^(q)   "
      f"error maximo = {err_bc3:.3e}")
print("      Cada parche conserva algo distinto y cada ley la rompe la")
print("      especie que migra en ese parche.")

##########################################################################
# Anillo dirigido frente a anillo simetrico de la misma intensidad
##########################################################################
# Se repite el escenario con un anillo SIMETRICO cuya traza es la misma,
# -(T1+T2+T3). Con tres aristas no dirigidas de tasa s cada patch tiene
# dos vecinos, de modo que la traza vale -6s y hay que tomar
# s = (T1+T2+T3)/6. Solo cambia la direccionalidad, no la intensidad
# total, de modo que cualquier diferencia proviene de la circulacion.
##########################################################################
print("\n" + "=" * 70)
print("ANILLO DIRIGIDO FRENTE A ANILLO SIMETRICO")
print("=" * 70)

s_sim = (T1 + T2 + T3) / 6.0

T_a_sim = np.zeros((num_patches, num_patches))
for p in range(num_patches):
    T_a_sim[p, (p + 1) % num_patches] = s_sim
    T_a_sim[p, (p - 1) % num_patches] = s_sim

t_s, X_s, _, info_s = simulate_metapopulation_dynamics(
    rn,
    rate=rate,
    modes='ode',
    num_patches=num_patches,
    spec_vector=spec_vector,
    x0_dict=x0_dict,
    transport_rates={'a': T_a_sim, 'c': T_c},
    t_span=t_span,
    n_steps=n_steps,
    dt_couple=dt_couple,
    rtol=rtol,
    atol=atol,
    return_info=True,
    verbose=False,
)

plot_species_dynamics_grid(t_s, X_s, species, num_patches=num_patches,
                           mode='per_species',
                           filename='MP3_4_anillo_simetrico_per_species.png')

L_a_sim = info_s['L_species']['a']
ev_sim = np.linalg.eigvals(L_a_sim)

print(f"\n  simetrico  s = {s_sim:.6g} en cada arista")
print(f"    autovalores: {np.round(ev_sim, 6)}")
print(f"    simetrica: {np.allclose(L_a_sim, L_a_sim.T)}   "
      f"|Im|max = {np.abs(ev_sim.imag).max():.3e}")
print(f"\n  dirigido   T1 = {T1}, T2 = {T2}, T3 = {T3}")
print(f"    autovalores: {np.round(np.linalg.eigvals(L_a_num), 6)}")
print(f"    simetrica: {np.allclose(L_a_num, L_a_num.T)}   "
      f"|Im|max = {np.abs(np.linalg.eigvals(L_a_num).imag).max():.6g}")
print(f"\n    traza identica en ambos casos: {np.trace(L_a_num):.6g} y "
      f"{np.trace(L_a_sim):.6g}")
print("    Misma intensidad total de transporte, misma velocidad global de")
print("    relajacion. Lo unico que cambia es que uno circula y el otro no.")

print(f"\n  {'':>10}{'a^(1)':>12}{'a^(2)':>12}{'a^(3)':>12}{'total a':>12}")
print(f"  {'dirigido':>10}{X['a'][-1, 0]:12.5f}{X['a'][-1, 1]:12.5f}"
      f"{X['a'][-1, 2]:12.5f}{np.nansum(X['a'][-1]):12.5f}")
print(f"  {'simetrico':>10}{X_s['a'][-1, 0]:12.5f}{X_s['a'][-1, 1]:12.5f}"
      f"{X_s['a'][-1, 2]:12.5f}{np.nansum(X_s['a'][-1]):12.5f}")

print(f"\n  {'':>10}{'b total':>12}{'c total':>12}")
print(f"  {'dirigido':>10}{np.nansum(X['b'][-1]):12.5f}"
      f"{np.nansum(X['c'][-1]):12.5f}")
print(f"  {'simetrico':>10}{np.nansum(X_s['b'][-1]):12.5f}"
      f"{np.nansum(X_s['c'][-1]):12.5f}")

print("\n  El anillo dirigido impone el recorrido 2 -> 3 de la a producida")
print("  por autocatalisis, que es justo donde se la consume. El anillo")
print("  simetrico la reparte en ambos sentidos y una parte vuelve al")
print("  parche 1 sin pasar por el reciclador.")

##########################################################################
print("\n" + "=" * 70)
print("LECTURAS")
print("=" * 70)
print("""
  Circulacion. Con tres parches y un anillo dirigido el espectro de L_a
  tiene un par complejo conjugado y la corriente estacionaria es la misma
  constante en las tres aristas. La distribucion no cambia y sin embargo
  hay flujo neto girando: es circulacion pura, y es exactamente lo que el
  caso de dos parches no puede producir, porque det(L) = 0 obliga a que el
  espectro sea real y la unica firma de la asimetria sea un estacionario
  sesgado.

  Dos transportes a la vez. a circula, c difunde y b no se mueve. Los tres
  regimenes conviven en el mismo sistema sin interferirse en el operador,
  porque el transporte es diagonal por especie. Se acoplan solo a traves
  de las reacciones locales.

  Dos leyes locales, ninguna global. El parche 1 conserva a+b y el parche
  3 conserva b+c. Son distintas, el parche 2 no comparte ninguna y la
  interseccion de los tres nucleos izquierdos es trivial. Cada ley local se
  rompe por la especie que migra en su parche, y el balance de migracion da
  la magnitud exacta de la rotura.

  Clausura espacial. Ningun parche cierra un ciclo. El parche 2 produce a,
  el anillo la lleva al parche 3, que la convierte en c, y la cadena
  difunde esa c hasta el parche 1, que regenera a desde b. La organizacion
  que cierra no esta en ninguna red local: esta en las tres mas el
  transporte. Es el punto que interesa desde COT, porque una clausura que
  solo existe a nivel espacial no se detecta analizando los parches por
  separado.
""")
print("=" * 70)
