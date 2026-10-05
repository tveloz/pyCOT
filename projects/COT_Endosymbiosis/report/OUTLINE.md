# Report outline (working)

Title (working): "Endosymbiotic fusion of elementary persistent modules: a
computable, structural model of major evolutionary transitions"

1. **Motivation** — the E. coli trace-metal fusion finding, generalized into
   two questions (fusigenic inflow search; endosymbiosis as extreme fusion).
2. **Framework** — from THEORY.md: fusigenic score definition, compartmentalized
   union + interface + synergy construction, complexification metrics
   (hybrid EPMs, ΔE0, Δmax size, count fusion vs proliferation).
3. **Toy model** — glycolysis host / TCA+ETC+ISC symbiont, hand-verified.
   Table: host alone / symbiont alone / merged. Figure: network + hybrid EPM
   highlight.
4. **Fusigenic inflow search** — toy validation (fescluster_h rediscovered
   independently of the merge machinery — nice cross-check). Real e_coli_core
   exhaustive sweep: pathway-intermediate shortcuts (2pg/3pg/pep, TCA
   intermediates) beat trace metals for THIS reconstruction — a real,
   reconstruction-detail-dependent finding, contrasted with the companion
   paper's iAF1260/iJO1366/iML1515 native-medium trace-metal result. Targeted
   iAF1260 trace-metal check for direct cross-validation.
5. **Real endosymbiosis mergers**
   5a. Mitochondrial-type (e_coli_core self-derived: fermentative host vs
       full aerobic symbiont). 100% hybrid EPMs, mean size 15.1→28.2,
       E0 unchanged — "re-integration" mechanism.
   5b. Chloroplast-type (S. cerevisiae host / Synechocystis symbiont,
       genuinely cross-organism). Dark vs light symbiont variants.
   Table comparing both real cases + toy model on the same metrics.
6. **Two mechanisms of endosymbiotic complexification** (the theoretical
   payoff): capability-gain (E0 grows — toy model) vs. re-integration
   (E0 constant, EPM hybridization — mitochondrial-type). Discuss whether
   chloroplast-type falls into one or the other or a third pattern.
7. **Discussion**
   - Falsifiable predictions: gene-retention in reduced organelle/symbiont
     genomes should track "interface load-bearing-ness" (which we can compute
     via interface-reaction knockout, per THEORY.md §2.3).
   - Limits: steady-state only, no dynamics/fitness; interface and synergy
     reactions are curated hypotheses, not derived; BiGG reconstructions are
     modern free-living relatives, not literal ancestral lineages.
   - Relationship to Question 1: the search independently rediscovering the
     toy model's designed dependency is itself a small piece of evidence the
     two frameworks are capturing the same real structural phenomenon from
     two different angles.
   - Future work: interface-richness sweep (single-metabolite to fully-open),
     interface-fragility/knockout analysis, extending to the 108-organism
     catalogue for many host/symbiont pairs, ESPM-level (not just EPM)
     merger analysis once genome-scale ESPM is tractable.

## Figures planned
- Fig 1: toy model network + hybrid EPM (reuse style from companion paper's
  fig0_illustrative_example).
- Fig 2: fusigenic-search ranked bar chart (toy + e_coli_core).
- Fig 3: before/after EPM count+size comparison across all cases (toy,
  mitochondrial-type, chloroplast-type dark, chloroplast-type light) —
  grouped bar/dot plot, same visual language as companion paper's fig1/fig4.
- Fig 4 (maybe): interface-fragility / knockout result if time permits.

## Status: drafting now, will fill in numbers as remaining runs complete.
See PROGRESS.md for the authoritative up-to-date numbers; this file is the
structural skeleton only.
