#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 18:55:14 2023

@author: tveloz
"""
import re
import copy
import random as rm
from itertools import combinations
from numbers import Real

import numpy as np
import pandas as pd
from bs4 import BeautifulSoup as bs
from bitarray import bitarray
from bitarray import frozenbitarray as fbt
# from scipy.optimize import linprog
import matplotlib.pyplot as plt
# from pyvis.network import Network
import networkx as nx

# from pyCOT_relational_properties import *
from Display import *


########################################################################################################
#############Vector of positions to bt transformations##################################################
########################################################################################################

def get_id_from_bt(bt: bitarray) -> list[int]:
    """Function that returns a vector from bitarray representation."""
    vec = [i for i in range(len(bt)) if bt[i] == 1]
    return vec

def get_bt_from_id(vec: Sequence[int], size) -> bitarray:
    """Function that returns bitarray from vector representation."""
    bt_array = bitarray(size)
    bt_array.setall(0)
    for i in vec:
        bt_array[i] = 1
    return bt_array

def get_bt_abstraction_from_vector(vec: Sequence[Real], t: int):
    """Function that returns bitarray of species with value larger than t in a vector."""
    bt = bitarray()
    for i in range(len(vec)):
        if vec[i] > t:
            bt.append(True)
        else:
            bt.append(False)
    return(bt)

def get_bt_from_names(names: Sequence[str], names_catalog: Sequence[str]) -> bitarray:
    """Function that returns bitarray from a sequence of names."""
    bt = bitarray(len(names_catalog))
    bt.setall(0)
    for name in names:
        bt[names_catalog.index(name)] = 1
    return bt

def get_bt_from_species(species: str | Sequence[str], species_catalog: Sequence[str]) -> bitarray:
    """
    Function that returns bitarray from a sequence of species names.

    Parameters
    ----------
    species : str or Sequence[str]
        Sequence of species names to represent as bitarray.
    species_catalog : Sequence[str]
        Sequence of all species names.

    Returns
    -------
    bt : bitarray
        Bitarray representation of species.
    """
    if isinstance(species, str):
        species = [species]
    
    for specie in species:
        if specie not in species_catalog:
            raise ValueError("species is not a sequence of recognized species")
        
    return get_bt_from_names(species, species_catalog)

def get_bt_from_reactions(reactions: str | Sequence[str], reaction_catalog: Sequence[str]) -> bitarray:
    """
    Function that returns bitarray from a sequence of reaction names.

    Parameters
    ----------
    reactions : str or Sequence[str]
        Sequence of reaction names to represent as bitarray.
    reaction_catalog : Sequence[str]
        Sequence of all reaction names.

    Returns
    -------
    bt : bitarray
        Bitarray representation of reaction.
    """
    if isinstance(reactions, str):
        reactions = [reactions]
    
    for reaction in reactions:
        if reaction not in reaction_catalog:
            raise ValueError("reaction is not a sequence of recognized reactions")
        
    return get_bt_from_names(reactions, reaction_catalog)
    