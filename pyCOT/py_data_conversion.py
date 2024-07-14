#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 18:55:14 2023

@author: tveloz
"""
from collections.abc import Sequence
from numbers import Real

from bitarray import bitarray
from numpy import ndarray

########################################################################################################
#############Vector of positions to bt transformations##################################################
########################################################################################################

def get_bt_abstraction_from_vector(vec: Sequence[Real], t: Real) -> bitarray:
    """Function that returns bitarray of species with value larger than t in a vector."""
    bt = bitarray()
    for i in range(len(vec)):
        if vec[i] > t:
            bt.append(True)
        else:
            bt.append(False)
    return bt

# TODO: Reconsiderar "names" en la signature de la función
def get_bt_from_names(names: Sequence[str], names_catalog: Sequence[str]) -> bitarray:
    """Function that returns bitarray from a sequence of names."""
    bt = bitarray(len(names_catalog))
    bt.setall(0)
    for name in names:
        bt[names_catalog.index(name)] = 1
    return bt

# TODO: Simplificar get_bt_from_species and from_reactions por get_bt_from_names
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
    
def filter_by_presence(elements: Sequence, presence_array: bitarray) -> Sequence:
    if len(elements) != len(presence_array):
        raise ValueError("elements and presence must have the same length")
    return [elements[i] for i, presence in zip(elements, presence_array) if presence]

def get_species_from_bt(species_presence: bitarray, species_catalog: Sequence[str]) -> list[str]:
    """
    Obtain the corresponding species names from a bitarray.

    Parameters
    ----------
    species_presence : bitarray
        Bitarray of species presence.
    species_catalog : Sequence[str]
        Sequence of all species names.

    Returns
    -------
    list[str]
        List of species names.
    """
    if len(species_catalog) != len(species_presence):
        raise ValueError("species_catalog and species_presence must have the same length")
    return filter_by_presence(species_catalog, species_presence)

def get_reactions_from_bt(reactions_presence: bitarray, reaction_catalog: Sequence[str]) -> list[str]:
    """
    Obtain the corresponding reaction names from a bitarray.

    Parameters
    ----------
    reactions_presence : bitarray
        Bitarray of reaction presence.
    reaction_catalog : Sequence[str]
        Sequence of all reaction names.

    Returns
    -------
    list[str]
        List of reaction names.
    """
    if len(reaction_catalog) != len(reactions_presence):
        raise ValueError("reaction_catalog and reactions_presence must have the same length")
    return filter_by_presence(reaction_catalog, reactions_presence)

def get_species_bt_from_reaction(reaction_name: str, reaction_set: Sequence[str], matrix: ndarray, t: Real = 0) -> bitarray:
    """
    Function that returns a bitarray for the reaction' support or products given a threshold t.
    
    Parameters
    ----------
    reaction_name : str
        Name of the reaction.
    reaction_set : Sequence[str]
        Sequence of all reaction names.
    matrix : ndarray
        Matrix of the reactions' support or products.
    t : Real
        Threshold.

    Returns
    -------
    bt : bitarray
        Bitarray representation of reaction's support or products.
    """

    if reaction_name not in reaction_set:
        raise ValueError(f"Reaction '{reaction_name}' not found in the reactions set.")
    
    reaction_index = reaction_set.index(reaction_name)
    species = matrix[reaction_index]
    return get_bt_abstraction_from_vector(species, t)