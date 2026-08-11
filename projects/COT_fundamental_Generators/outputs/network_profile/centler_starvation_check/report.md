# Network structure profile: centler_starvation_check

Source file: `data/Examples_tests/Centler2006_EcoliSugar/centler_starvation.txt`

## Step 1: ERCs and req/prod structure

- Species: 92   Reactions: 168   ERCs: 23
- P-ERCs (persistent, req=0): 10
- req size per ERC: mean=0.57, max=1
- prod size per ERC: mean=8.13, max=63
- Explicit inflow species (25): ADP, AMP, ATP, PromCrp, PromCya, PromEI, PromEIIA, PromEIIBC, PromFbp, PromFda, PromGap, PromGlcT, PromGlk, PromGlpD, PromGlpFK, PromGlpR, PromGpm, PromHPr, PromLacI, PromLacZY, PromPfk, PromPgi, PromPyk, PromTpi, RNAP
- Explicit outflow species (65): ADP, AMP, ATP, Allo, Crp, CrpmRNA, Cya, CyamRNA, DHAP, EI, EIIA, EIIAP, EIIAmRNA, EIIBC, EIIBCmRNA, EImRNA, FBP, Fbp, FbpmRNA, Fda, FdamRNA, Fru6P, G3P, Gap, GapmRNA, Glc, Glc6P, GlcT, GlcTmRNA, Glk, GlkmRNA, GlpD, GlpDmRNA, GlpF, GlpFKmRNA, GlpFKmRNA1, GlpK, GlpR, GlpRmRNA, Gly ... (truncated)
- Structural-only sources (produced, never consumed, no explicit inflow reaction): (none)
- Structural-only sinks (consumed, never produced, no explicit outflow reaction): Glcex, Lacex

## Step 2: ERC hierarchy and fundamental relations

- Hasse (containment) edges: 23
- Comparable / incomparable ERC pairs: 61 / 192
- Fundamental synergies: 4  (density 0.17 per ERC)
- Fundamental complementarities: 13  (density 0.57 per ERC)
- Max synergy out-degree: 2

Files in this directory: `erc_table.csv`, `req_prod_histogram.png`, `synergy_complementarity_degree_histogram.png`.

(No organization / EPM / ESPM computation performed -- structural profile only.)
