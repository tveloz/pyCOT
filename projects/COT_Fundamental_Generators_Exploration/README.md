# COT Fundamental Generators — Exploration

This project holds exploratory, pedagogical, and correctness-validation work built on the fundamental-generators (ERC/EPM/ESPM) approach to computing chemical organizations, targeted at small networks (on the order of 100 reactions at most). It includes oracle-based validation scripts, interactive single-network deep-dives, hierarchy/degeneracy exploration tools, and a reproduction of the Centler et al. (2006) worked example.

The core computational engine (types, closure, ERC discovery, hierarchy, synergy, complementarity, fundamental graph, EPM/ESPM search, maxSemiOrganization) now lives in `src/pyCOT/analysis/organizations/`. This project consumes that library — it does not reimplement it. The `cot_gen/` package here is project-specific reporting/exploration tooling built on top of the core library, and `oracles/` provides independent brute-force reference implementations used to validate the core engine against small, hand-checkable networks.

Note: `cot_gen/` is also imported cross-project by `projects/Decomposition_Theorem/` and `projects/RAF_Comparison/` (via `sys.path` insertion) and by `projects/COT_Fundamental_Generators_Complex/benchmarks/bigg_synergy_explore.py` — moving or renaming this package requires updating those dependents.
