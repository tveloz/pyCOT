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

def get_bt_abstraction_from_vector(vec: Sequence[int], t: int):
    """Function that returns bitarray of species with value larger than t in a vector."""
    bt = bitarray()
    for i in range(len(vec)):
        if vec[i] > t:
            bt.append(True)
        else:
            bt.append(False)
    return(bt)
