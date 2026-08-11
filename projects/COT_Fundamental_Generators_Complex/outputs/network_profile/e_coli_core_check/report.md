# Network structure profile: e_coli_core_check

Source file: `data/biochemical_databases/BiGG/bigg_e_coli_core.txt`

## Step 1: ERCs and req/prod structure

- Species: 72   Reactions: 141   ERCs: 49
- P-ERCs (persistent, req=0): 34
- req size per ERC: mean=0.33, max=2
- prod size per ERC: mean=6.80, max=53
- Explicit inflow species (9): co2_e, glc__D_e, h2o_c, h2o_e, h_e, nh4_e, o2_e, pep_c, pi_e
- Explicit outflow species (22): ac_e, acald_e, akg_e, co2_e, etoh_e, for_e, fru_e, fum_e, glc__D_e, gln__L_e, glu__L_e, h2o_c, h2o_e, h_e, lac__D_e, mal__L_e, nh4_e, o2_e, pep_c, pi_e, pyr_e, succ_e
- Structural-only sources (produced, never consumed, no explicit inflow reaction): (none)
- Structural-only sinks (consumed, never produced, no explicit outflow reaction): (none)

## Step 2: ERC hierarchy and fundamental relations

- Hasse (containment) edges: 67
- Comparable / incomparable ERC pairs: 124 / 1052
- Fundamental synergies: 110  (density 2.24 per ERC)
- Fundamental complementarities: 22  (density 0.45 per ERC)
- Max synergy out-degree: 16

Files in this directory: `erc_table.csv`, `req_prod_histogram.png`, `synergy_complementarity_degree_histogram.png`.

(No organization / EPM / ESPM computation performed -- structural profile only.)
