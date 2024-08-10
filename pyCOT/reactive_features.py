#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 19:06:00 2024

Functions that help to analyze basic features of the reaction network

@author: tveloz
"""
# import time

# import networkx as nx
# import matplotlib.pyplot as plt
import pandas as pd

from pyCOT.reaction_network import *
from pyCOT.closure_structure import *
from pyCOT.file_manipulation import *

def reac_analysis(RN):
    columns=RN.SpStr
    indexes=('inflow','outflow','trans','synth','decomp','single_rep','multi_rep')
    # Create a DataFrame with a single row of zeros and specified index
    # Set each value in the first row to zero individually
    df_sp = pd.DataFrame(columns=columns,index=indexes)
    df_sp.columns = columns
    df_sp.index = indexes
    for col in columns:
        for i in indexes:
            df_sp.loc[i, col] = 0.0
    # Display the DataFrame
    
    for r in RN.RnStr:
        supp=RN.get_supp_from_reactions(r)
        prod=RN.get_prod_from_reactions(r)
        if supp==[]:
            for s in prod:
                df_sp.loc["inflow",s]=df_sp.loc["inflow",s]+1.0
        elif prod==[]:
            for s in supp:
                df_sp.loc["outflow",s]=df_sp.loc["outflow",s]+1.0
        else:
            if (len(supp)==1 and len(prod)==1):
                for s in set(supp).union(set(prod)):
                    df_sp.loc["trans",s]=df_sp.loc["trans",s]+1.0
            elif len(supp)>len(prod):
                for s in set(supp).union(set(prod)):
                    df_sp.loc["synth",s]=df_sp.loc["synth",s]+1.0
            elif len(supp)<len(prod):
                for s in set(supp).union(set(prod)):
                    df_sp.loc["decomp",s]=df_sp.loc["synth",s]+1.0
            else: 
                diff=sum(1 for sp, pr in zip(supp, prod) if sp != pr)
                if diff==1:
                    for s in set(supp).union(set(prod)):
                        df_sp.loc["single_rep",s]=df_sp.loc["single_rep",s]+1.0
                else:
                    for s in set(supp).union(set(prod)):
                        df_sp.loc["multi_rep",s]=df_sp.loc["multi_rep",s]+1.0
    return df_sp
    