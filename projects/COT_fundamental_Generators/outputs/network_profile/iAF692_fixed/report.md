# Network structure profile: iAF692_fixed

Source file: `data/biochemical_databases/BiGG/bigg_iAF692.txt`

## Step 1: ERCs and req/prod structure

- Species: 628   Reactions: 900   ERCs: 517
- P-ERCs (persistent, req=0): 169
- req size per ERC: mean=0.84, max=3
- prod size per ERC: mean=11.44, max=370
- Explicit inflow species (29): aicar_c, chor_c, co2_c, co2_e, cobalt2_e, cys__L_e, dhpt_c, dkdofp_c, fum_c, g3p_c, h2o_c, h2o_e, h2s_e, h_c, h_e, indole_c, lys__L_c, meoh_e, n2_e, na1_e, nac_e, nh4_e, ni2_e, pep_c, pi_c, pi_e, ppbng_c, ppi_c, so3_e
- Explicit outflow species (73): ac_e, actn__R_e, ala__L_e, alac__S_e, btn_e, ca2_e, cbi_e, cbl1_e, cbl1hbi_e, cd2_e, ch4_e, ch4s_e, cit_e, cl_e, co2_e, co_e, cobalt2_e, cu2_e, cys__L_e, dma_e, dms_e, dohdu_c, etha_e, fe2_e, fe3_e, fol_e, gcald_e, glcn_e, glu1sa_c, glu__L_e, gly_e, glyald_e, glyb_e, glyc_e, h2_e, h2o_c, h2o_e, h_e, ile__L_e, ind3ac_e ... (truncated)
- Structural-only sources (produced, never consumed, no explicit inflow reaction): 2ins_c, 3uib_c, 4hba_c, 5mdr1p_c, 5odhf2a_c, 6ax_c, 6pgl_c, 6pthp_c, aconm_c, alatrna_c, argtrna_c, asptrna_c, btamp_c, cala_c, camp_c, cmaphis_c, cystrna_c, dha_c, dtdprmn_c, fe3_c ... (truncated)
- Structural-only sinks (consumed, never produced, no explicit outflow reaction): 2pglyc_c, 4cml_c, 4mhetz_c, 56dthm_c, 56dura_c, 5mta_c, 6ax6ax_c, acon_T_c, ahdt_c, atrz_c, caphis_c, csn_c, dxyl5p_c, fald_c, fmn_c, fru_c, gua_c, hgbam_c, iasp_c, nmn_c ... (truncated)

## Step 2: ERC hierarchy and fundamental relations

- Hasse (containment) edges: 733
- Comparable / incomparable ERC pairs: 2677 / 130709
- Fundamental synergies: 882  (density 1.71 per ERC)
- Fundamental complementarities: 425  (density 0.82 per ERC)
- Max synergy out-degree: 91

Files in this directory: `erc_table.csv`, `req_prod_histogram.png`, `synergy_complementarity_degree_histogram.png`.

(No organization / EPM / ESPM computation performed -- structural profile only.)
