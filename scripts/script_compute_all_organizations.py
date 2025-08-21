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
from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.io.functions import read_txt

# Import the reaction network library and persistent modules
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

# Load a reaction network
file_path = 'networks/testing/Farm.txt'
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
    result = compute_all_organizations(
        RN, 
        max_module_size=10,
        max_generator_size=8,
        verbose=True
    )
    
    elementary_sos, elementary_organizations, all_organizations, all_semi_organizations, statistics, computation_data = result
    
    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    
    print(f"Elementary semi-organizations: {len(elementary_sos)}")
    print(f"Elementary organizations: {len(elementary_organizations)}")
    print(f"Total organizations: {len(all_organizations)}")
    
    if all_organizations:
        print(f"\nFirst few organizations:")
        for i, org in enumerate(all_organizations[:3]):
            species_count = len(org.closure_species)
            erc_count = len(org.constituent_ercs)
            print(f"  Organization {i+1}: {species_count} species, {erc_count} ERCs")
    
    print("\n✅ Test completed successfully!")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

