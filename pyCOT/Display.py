#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 15:31:23 2023

@author: tveloz
"""
import re

import numpy as np
from bitarray import bitarray
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt

from pyCOT_constructor import *

def parse_chemical_reactions(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    species_set = set(re.findall(r'\b[a-zA-Z]\b', ''.join(lines)))
    species_mapping = {species: index for index, species in enumerate(sorted(species_set))}
    
    sp_bt = bitarray([False] * len(species_mapping))
    sp_vec = np.zeros(len(species_mapping), dtype=int)
    sp_str = list(species_mapping.keys())

    rn_bt_s = bitarray()
    rn_bt_p = bitarray()
    rn_vec_s_list = []
    rn_vec_p_list = []
    rn_str_list = []

    for index, line in enumerate(lines):
        reaction_str = line.strip()
        reactants, products = map(str.strip, reaction_str.split('->'))

        rn_bt_s.append(True)  # Assume all reactions require support for simplicity
        rn_bt_p.append(True)  # Assume all reactions produce support for simplicity

        rn_vec_s = np.zeros(len(species_mapping), dtype=int)
        rn_vec_p = np.zeros(len(species_mapping), dtype=int)

        for reactant in re.findall(r'\b[a-zA-Z]\b', reactants):
            if reactant in species_mapping:
                rn_vec_s[species_mapping[reactant]] += 1

        for product in re.findall(r'\b[a-zA-Z]\b', products):
            if product in species_mapping:
                rn_vec_p[species_mapping[product]] += 1

        rn_vec_s_list.append(rn_vec_s)
        rn_vec_p_list.append(rn_vec_p)
        rn_str_list.append(f'r_{index+1}')  # Label each reaction as r_i where i is the index

    return ReactionNetwork(SpBt=sp_bt, SpVec=sp_vec, SpStr=sp_str, RnBtS=rn_bt_s, RnBtP=rn_bt_p,
                 RnVecS=np.array(rn_vec_s_list), RnVecP=np.array(rn_vec_p_list), RnStr=rn_str_list)

def build_graph(reaction_network):
    G = nx.DiGraph()
    i=0;
    for reaction in reaction_network: 
        if len(reaction) == 2:  # Transformation reaction
            reactant, product = reaction
            ri="r"+str(i)
            G.add_edge(reactant, ri)
            G.add_edge(ri, product)
            i=i+1
        elif len(reaction) == 3:  # Synthesis reaction
            reactant1, reactant2, product = reaction
            ri="r"+str(i)
            G.add_edge(reactant1, ri)
            G.add_edge(reactant2, ri)
            G.add_edge(ri,product)
            i=i+1
    return G

def visualize_graph(G):
    pos = nx.circular_layout(G)  # Position nodes in a circular arrangement
    pos = graphviz_layout(G, prog="dot")
    nx.draw(G, pos, with_labels=True, node_size=150, node_color="skyblue", font_size=6, font_weight="bold")
    plt.title("Reaction Network Graph")

    # Enable interactive mode
    plt.ion()

    def on_press(event):
        if event.inaxes is None:
            return
        if event.button == 1:  # Left mouse button
            node_labels = nx.get_node_attributes(G, 'label')
            clicked_node = None
            for node, coords in pos.items():
                if (event.xdata - coords[0])**2 + (event.ydata - coords[1])**2 < 0.02:  # Tolerance for clicking
                    clicked_node = node
                    break
            if clicked_node is not None:
                print(f"Clicked node: {clicked_node}")

    plt.gcf().canvas.mpl_connect('button_press_event', on_press)
    plt.show()
