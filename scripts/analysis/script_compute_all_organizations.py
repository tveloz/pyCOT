#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple test script for compute_all_organizations function
"""
import os
import sys

# Add the root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import pyCOT modules

from pyCOT.visualization.plot_dynamics import plot_series_ode
# Import required modules
from pyCOT.analysis.Persistent_Modules_Generator import compute_all_organizations  # Organization computation
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import time
import json
from pathlib import Path
from pyCOT.core.rn_rustworkx import ReactionNetwork
from pyCOT.io.functions import read_txt
from pyCOT.visualization.rn_visualize import hierarchy_visualize_html, rn_visualize_html

# Import the reaction network library and persistent modules
from pyCOT.analysis.ERC_Hierarchy import ERC, ERC_Hierarchy
from pyCOT.analysis.ERC_Synergy_Complementarity import build_erc_sorn
from pyCOT.analysis.Persistent_Modules_Generator import (
    ElementarySO,
    Organization,
    compute_elementary_sos
)

"""Simple test of the hierarchical organization computation."""


print("Loading reaction network...")

# Load a reaction network
#file_path = 'networks/testing/Farm.txt'
#file_path = 'networks/Navarino/RN_IN_05.txt'
file_path = 'networks/Conflict_Theory/Resource_Community_Insurgency_Loops_model3.txt'
RN = read_txt(file_path)

# Verify we got a proper ReactionNetwork object
if not isinstance(RN, ReactionNetwork):
    raise TypeError(f"read_txt returned {type(RN)}, expected ReactionNetwork")

print(f"✅ Loaded network: {len(RN.species())} species, {len(RN.reactions())} reactions")

print("\n" + "="*60)
print("COMPUTING ALL ORGANIZATIONS")
print("="*60)

# Run the main computation
try:
    results = compute_all_organizations(
        RN,
        max_generator_size=15,
        max_organization_size=20,
        verbose=True
    )

    # Extract results from dictionary
    elementary_sos = results['elementary_sos']
    elementary_organizations = results['elementary_organizations']
    all_organizations = results['all_organizations']
    statistics = results['statistics']

    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)

    print(f"Elementary semi-organizations: {len(elementary_sos)}")
    print(f"Elementary organizations: {len(elementary_organizations)}")
    print(f"Total organizations: {len(all_organizations)}")

    if all_organizations:
        print(f"\nOrganizations found:")
        for i, org in enumerate(all_organizations):
            # Handle both elementary (closure_species) and non-elementary (combined_closure)
            if hasattr(org, 'combined_closure'):
                species_names = [sp.name for sp in org.combined_closure]
                org_type = "non-elementary"
            else:
                species_names = [sp.name for sp in org.closure_species]
                org_type = "elementary"
            print(f"  Organization {i+1} [{org_type}]: {len(species_names)} species - {set(species_names)}")

    print("\n✅ Test completed successfully!")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

