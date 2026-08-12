# Inflow-regime analysis: e_coli_core

Source: `data/biochemical_databases/BiGG/bigg_e_coli_core.txt`

Scenarios: 4 (3 named regime(s) + the isolated baseline)

Organizations reported below are **fundamental** organizations only (ERCs combined via fundamental synergy/complementarity, per Veloz & Bassi's productive-novelty theory) -- not the full naive enumeration. `spurious lower bound` is a cheap, non-enumerating count of additional (uncomputed) organizations implied by free-species combinatorics alone; see organizations.py's module docstring.

| scenario | food species | ERCs | P-ERCs | req mean/max | prod mean/max | synergies | complementarities | fundamental orgs (sizes) | spurious lower bound | Hasse diagram |
|---|---|---|---|---|---|---|---|------|
| isolated (no inflow) | (none) | 77 | 38 | 0.71 / 3 | 5.79 / 63 | 499 | 124 | 12 [2, 3, 3, 4, 4, 5, 6, 7, 8, 10, 10, 14] | 21,069,022,323 | [hasse_isolated_no_inflow.png](hasse_isolated_no_inflow.png) (12/12 agree) |
| aerobic_core | co2_e, glc__D_e, h2o_e, h_e, nh4_e, o2_e, pi_e | 49 | 34 | 0.33 / 2 | 6.80 / 53 | 110 | 22 | 12 [13, 15, 15, 16, 16, 17, 26, 43, 57, 59, 66, 68] | 241,447 | [hasse_aerobic_core.png](hasse_aerobic_core.png) (12/12 agree) |
| anaerobic_core | co2_e, glc__D_e, h2o_e, h_e, nh4_e, pi_e | 50 | 35 | 0.34 / 2 | 6.68 / 53 | 119 | 24 | 12 [11, 13, 13, 14, 14, 15, 24, 41, 55, 57, 64, 66] | 505,093 | [hasse_anaerobic_core.png](hasse_anaerobic_core.png) (12/12 agree) |
| carbon_starvation_core | co2_e, h2o_e, h_e, nh4_e, o2_e, pi_e | 51 | 34 | 0.35 / 2 | 6.39 / 53 | 117 | 22 | 7 [12, 14, 14, 15, 15, 16, 25] | 249,815 | [hasse_carbon_starvation_core.png](hasse_carbon_starvation_core.png) (7/7 agree) |

(No EPM/ESPM/organization computation performed unless --organizations was passed; Phases A/B are req/prod structural statistics only, matching the network's own wiring under each regime.)

