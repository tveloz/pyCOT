# COT Fundamental Generators — Complex

This project holds batch/sweep statistics, catalogue-building, and visualization work analyzing complex (large, genome-scale) reaction networks with the fundamental-generators (ERC/EPM/ESPM) approach — for example the BiGG and BioModels genome-scale collections, with networks ranging into the thousands of reactions. It includes batch sweep runners, resumable ESPM computation with checkpointing, risk/economics analyses of dead-end species, and network-structure profiling and visualization tools.

The core computational engine (types, closure, ERC discovery, hierarchy, synergy, complementarity, fundamental graph, EPM/ESPM search, maxSemiOrganization) now lives in `src/pyCOT/analysis/organizations/`. This project consumes that library — it does not reimplement it. `benchmarks/bigg_synergy_explore.py` additionally imports the `cot_gen` reporting package from the sibling `projects/COT_Fundamental_Generators_Exploration/` project.
