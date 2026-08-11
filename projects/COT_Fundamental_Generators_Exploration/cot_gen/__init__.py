"""
cot_gen — COT_Fundamental_Generators_Exploration-specific reporting/exploration tools.

The core generative-organization engine (types, closure, ERC discovery,
hierarchy, synergy, complementarity, fundamental graph, EPM/ESPM search,
maxSemiOrganization, generators, metrics, io) has moved to the pyCOT core
library at `pyCOT.analysis.organizations` — see that package's docstring
for the full pipeline and a quick-start example. This package now only
holds project-specific tooling built on top of that engine:

  deep_report / deep_report_viz  : instrumented EPM/ESPM search + degeneracy
                                    statistics, HTML/PNG report generation
  explorer                        : interactive GraphExplorer for manual
                                    ERC-by-ERC generator exploration
  metanetwork                     : ERC/hierarchy/synergy/complementarity/
                                    generators summary view (MetaNetwork)
  results_io                      : CSV row builder / writer for batch runs

Implements the theory from Veloz & Bassi (2025) "Synergy and Complementarity:
The Generative Basis of Chemical Organizations".
"""
__version__ = "0.1.0"
