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
from bitarray import bitarray as bt
from bitarray import frozenbitarray as fbt
# from scipy.optimize import linprog
import matplotlib.pyplot as plt
# from pyvis.network import Network
import networkx as nx

# from pyCOT_relational_properties import *
from Display import *


# Function that returns bitarray from vector representation    

