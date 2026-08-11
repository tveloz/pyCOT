"""
pyCOT.analysis.organizations — Genome-scale organization computation.

Implements the theory from Veloz & Bassi (2025) "Synergy and Complementarity:
The Generative Basis of Chemical Organizations", ported from the
projects/COT_Fundamental_Generators_Exploration/cot_gen research prototype into the
core library as its canonical, actively-maintained home.

Pipeline (see organizations.py for the full narrative):
  Stage 0: Preprocessing & E0            io_pyCOT, cot_types, closure_oracle
  Stage 1: Closure + ERCs + MinBas       closure, erc
  Stage 2: ERC Hierarchy                 hierarchy
  Stage 3: Fundamental synergies         synergy
  Stage 4: Fundamental complementarities complementarity
  Stage 5: EPM / ESPM exploration        fundamental_graph, epm
  Stage 6: LP self-maintenance           self_maintenance
  Orchestration                          organizations.compute_organizations

Quick start
-----------
    from pyCOT.io.functions import read_txt
    from pyCOT.analysis.organizations import compute_organizations

    rn = read_txt("data/Examples_tests/Centler2006_EcoliSugar/centler_glucose.txt")
    result = compute_organizations(rn, network_id="centler_glucose")
    print(f"{len(result.organizations)} verified organizations "
          f"out of {len(result.semiorganizations)} semi-organizations")
"""
from .cot_types import BitSet, RNData, ERCData
from .io_pyCOT import build_rndata, load_rndata
from .erc import compute_ercs
from .hierarchy import build_hierarchy, HierarchyData
from .synergy import compute_synergies_basis_first, SynergyResult
from .complementarity import compute_complementarities, CompResult
from .epm import compute_epms, compute_espm, EPMResult, ESPMResult
from .max_semiorg import compute_max_semiorganization, max_semiorganization_species
from .self_maintenance import minimize_sv, check_self_maintenance, diagnose_self_maintenance
from .organizations import compute_organizations, OrganizationsResult, SemiOrganization

__all__ = [
    "BitSet", "RNData", "ERCData",
    "build_rndata", "load_rndata",
    "compute_ercs",
    "build_hierarchy", "HierarchyData",
    "compute_synergies_basis_first", "SynergyResult",
    "compute_complementarities", "CompResult",
    "compute_epms", "compute_espm", "EPMResult", "ESPMResult",
    "compute_max_semiorganization", "max_semiorganization_species",
    "minimize_sv", "check_self_maintenance", "diagnose_self_maintenance",
    "compute_organizations", "OrganizationsResult", "SemiOrganization",
]
