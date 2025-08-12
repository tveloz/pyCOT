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

from pyCOT.plot_dynamics import plot_series_ode
# Import required modules
from pyCOT.Persistent_Modules import compute_all_organizations  # Our new function
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import time
import json
from pathlib import Path

# Import the reaction network library and persistent modules
from pyCOT.io.functions import *
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy
from pyCOT.ERC_Synergy_Complementarity import build_erc_sorn
from pyCOT.Persistent_Modules import (
    compute_persistent_modules,
    build_irreducible_generators,
    IrreducibleGenerator, 
    ElementarySO,
    identify_p_ercs,
    analyze_generator_statistics,
    compute_elementary_sos
)

"""Simple test of the hierarchical organization computation."""

print("Loading reaction network...")

# Load a reaction network (adjust path as needed)
# You can replace this with your own network file
# Example: Load from file (replace with your network file)
file_path='networks/testing/Farm.txt'
file_path = 'networks/RandomAlife/RN_Ns_20_Norg_4_id_12.txt'
file_path = 'networks/RandomAlife/RN_Ns_40_Norg_12_id_358.txt'
#file_path = 'networks/Navarino/RN_IN_05.txt'

RN = read_txt(file_path)
print(f"✓ Loaded network: {len(RN.species())} species, {len(RN.reactions())} reactions")

print("\n" + "="*60)
print("COMPUTING ALL ORGANIZATIONS")
print("="*60)

# Run the main computation
elementary_sos, elementary_organizations, all_organizations, all_semi_organizations, statistics, computation_data = compute_all_organizations(
    RN, 
    max_module_size=100,  # Small size for testing
    max_generator_size=100,  # Small size for testing
    verbose=True
)

print("\n" + "="*60)
print("RESULTS SUMMARY")
print("="*60)

print(f"Total semi-organizations found: {len(all_semi_organizations)}")
print(f"Total organizations found: {len(all_organizations)}")

if all_organizations:
    print(f"\nOrganizations by size:")
    org_by_size = {}
    for org in all_organizations:
        size = org.size()
        org_by_size[size] = org_by_size.get(size, 0) + 1
    
    for size in sorted(org_by_size.keys()):
        count = org_by_size[size]
        print(f"  Size {size}: {count} organizations")
    
    print(f"\nFirst few organizations:")
    for i, org in enumerate(all_organizations[:3]):
        print(f"  Organization {i+1}:")
        print(f"    - Size: {org.size()} elementary SOs")
        print(f"    - Species: {len(org.combined_closure)} total")
        print(f"    - ERCs: {len(org.constituent_ercs)} total")

if all_semi_organizations:
    print(f"\nSemi-organizations by size:")
    so_by_size = {}
    for so in all_semi_organizations:
        size = so.size()
        so_by_size[size] = so_by_size.get(size, 0) + 1
    
    for size in sorted(so_by_size.keys()):
        count = so_by_size[size]
        print(f"  Size {size}: {count} semi-organizations")

# Show final statistics
final_stats = statistics['final_results']
conversion_rate = final_stats['conversion_rate']
print(f"\nConversion rate: {conversion_rate:.1%} (SOs → Organizations)")

print("\n" + "="*60)
print("TEST COMPLETE")
print("="*60)

# Convert all semi-organizations to species name sets
all_semi_organizations_sets = [
    set(sp.name for sp in so.combined_closure) 
    for so in all_semi_organizations
]


# Convert to species name sets
elementary_orgs_sets = [
    set(sp.name for sp in org.closure_species) 
    for org in elementary_organizations
]

# Separate elementary from non-elementary organizations

non_elementary_orgs = [org for org in all_organizations if org not in elementary_organizations]

non_elementary_orgs_sets = [
    set(sp.name for sp in org.combined_closure) 
    for org in non_elementary_orgs
]

# Use the converted sets for visualization with three color groups
from pyCOT.rn_visualize import hierarchy_visualize_html
hierarchy_visualize_html(
    all_semi_organizations_sets,
    lst_color_subsets=[
        ("yellow", all_semi_organizations_sets), # All semi-organizations in yellow
        ("green", elementary_orgs_sets),  # Elementary organizations in green
        ("cyan", non_elementary_orgs_sets)       # Non-elementary organizations in cyan     
    ]
)

print("✓ Hierarchy visualization saved as: hierarchy_visualization.html")


