#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 15:31:23 2023

@author: tveloz
"""
import os
import re

from bitarray import bitarray as bt
import numpy as np
from collections import OrderedDict
from bs4 import BeautifulSoup as bs
import requests

from pyCOT_constructor import *


def extract_species_and_reactions(file_path):
    species_set = OrderedDict()
    reactions_set = OrderedDict()

    with open(file_path, 'r') as file:
        reaction_lines = file.read().split(';')

    for line in reaction_lines:
        line = line.strip()
        if line and not line.startswith('#'):
            # Extract reaction name and equation
            parts = line.split(':')
            reaction_name = parts[0].strip()
            reaction_equation = parts[1].strip()

            # Extract species from reactants and products
            species = re.findall(r'(\d*)?([a-zA-Z_]\w*)', reaction_equation)
            for coefficient, species_name in species:
                species_set[species_name] = None  # Using None as a placeholder

            # Add reaction name to the set
            #reac_name = re.findall(r'[Rr](_*\d+)+', reaction_name)  
            reactions_set[reaction_name] = None  # Using None as a placeholder

    unique_species = list(species_set.keys())
    reactions = list(reactions_set.keys())
    return unique_species, reactions

def build_stoichiometric_vectors(file_path, species_set):
    reactants_vectors = []
    products_vectors = []

    with open(file_path, 'r') as file:
        reaction_lines = file.read().split(';')

    for line in reaction_lines:
        line = line.strip()
        if line and not line.startswith('#'):
            # Extract reaction equation
            parts = line.split(':')
            reaction_equation = parts[1].strip()
            reactants, products = reaction_equation.split('=>')
            reactants = reactants.strip()
            products = products.strip()

            # Initialize vectors for reactants and products
            reactants_vector = np.zeros(len(species_set), dtype=int)
            products_vector = np.zeros(len(species_set), dtype=int)

            # Extract species and coefficients from reactants and products
            species_and_coefficients_reactants = re.findall(r'(\d*)?([a-zA-Z_]\w*)', reactants)
            species_and_coefficients_products = re.findall(r'(\d*)?([a-zA-Z_]\w*)', products)
            
            for coefficient, species_name in species_and_coefficients_reactants:
                species_index = species_set.index(species_name)
                stoichiometric_coefficient = int(coefficient) if coefficient else 1
                reactants_vector[species_index] = stoichiometric_coefficient
            for coefficient, species_name in species_and_coefficients_products:
                species_index = species_set.index(species_name)
                stoichiometric_coefficient = int(coefficient) if coefficient else 1
                products_vector[species_index] = stoichiometric_coefficient    

            reactants_vectors.append(reactants_vector)
            products_vectors.append(products_vector)

    return reactants_vectors, products_vectors
def load_pyCOT_from_Txt(file_path):
    input_file_name = file_path    
    base_name, extension = os.path.splitext(input_file_name)
    if extension!='.txt':
        print("load_pyCOT_from_file: Wrong extension in input!")
        return 0
    species_set, reactions_list = extract_species_and_reactions(file_path)
    
    SpStr = list(species_set)
    SpBt = bt(len(SpStr))
    SpBt.setall(True)
    RnStr = reactions_list
    RnBt = bt(len(RnStr))
    RnBt.setall(True)
    RnMsupp, RnMprod = build_stoichiometric_vectors(file_path, species_set)
    
    return ReactionNetwork(SpStr, SpBt, RnStr, RnBt, RnMsupp, RnMprod)

def load_pyCOT_from_Sbml(file):
    
    # Input XML file name
    input_file_name = file
    base_name, extension = os.path.splitext(input_file_name)
    if extension!='.xml':
        print("load_pyCOT_from_Sbml: Wrong extension in input!")
        return 0
    # Open the XML file and read its contents
    with open(input_file_name, "r") as file:
        xml_content = file.read()
    
    # Parse the XML content using BeautifulSoup with the "xml" parser
    soup = bs(xml_content, "xml")
    
    # Extract file name without extension
    file_name_without_extension = os.path.splitext(input_file_name)[0]
    
    # Output file name with .txt extension
    output_file_name = file_name_without_extension + ".txt"
    
    # Extract all <reaction> tags
    reaction_tags = soup.find_all("reaction")
    
    # Create a new text file for writing
    with open(output_file_name, "w") as output_file:
        # Process each <reaction> tag
        for reaction_tag in reaction_tags:
            # Extract the reaction id and reversible attribute
            reaction_id = reaction_tag["id"]
            reversible = reaction_tag.get("reversible", "false") == "true"
    
            # Extract reactants, products, and modifiers for the reaction if <listOfReactants>, <listOfProducts>, and <listOfModifiers> tags exist
            reactants_tag = reaction_tag.find("listOfReactants")
            products_tag = reaction_tag.find("listOfProducts")
            modifiers_tag = reaction_tag.find("listOfModifiers")
            
            # Extract reactants, products, and modifiers with stoichiometry if the tags exist, otherwise set to empty lists with stoichiometry 1
            if reactants_tag:
                reactants = [(reactant["species"], int(reactant.get("stoichiometry", "1"))) for reactant in reactants_tag.find_all("speciesReference")]
            else:
                reactants = []
    
            if products_tag:
                products = [(product["species"], int(product.get("stoichiometry", "1"))) for product in products_tag.find_all("speciesReference")]
            else:
                products = []
    
            if modifiers_tag:
                modifiers = [(modifier["species"], int(modifier.get("stoichiometry", "1"))) for modifier in modifiers_tag.find_all("modifierSpeciesReference")]
            else:
                modifiers = []
    
            # Filter out empty categories
            reactants_str = " + ".join([f"{stoich}{species}" for species, stoich in reactants if stoich != 0])
            products_str = " + ".join([f"{stoich}{species}" for species, stoich in products if stoich != 0])
            modifiers_str = " + ".join([f"{stoich}{species}" for species, stoich in modifiers if stoich != 0])
            # Write the reaction id and the reaction as a sum of stoichiometry*species => sum of stoichiometry*species with modifiers on both sides of the reaction, ending with a semicolon
            reaction_line=reaction_id + ": " + (modifiers_str + " + " if modifiers_str else "") + reactants_str +  " => "+ (modifiers_str + " + "  if modifiers_str else "") + products_str + ";\n"
            reaction_line=reaction_line.replace("+  =>","=>")
            output_file.write(reaction_line)
            
            # Check if the reaction is reversible
            if reversible:
                # Write the backward reaction (products as reactants, reactants as products) with modifiers on both sides of the reaction
                output_file.write(reaction_id + "_b: " + (modifiers_str + " + " if modifiers_str else "") + products_str + " + " +  " => " + ( modifiers_str + " + " if modifiers_str else "") + reactants_str + ";\n")

    return(load_pyCOT_from_file(output_file_name)) 

def load_pyCOT_from_file(file):
    input_file_name = file
    
    base_name, extension = os.path.splitext(input_file_name)
    if extension=='.xml':
        return(load_pyCOT_from_Sbml(file))
    elif extension=='.txt':
        return(load_pyCOT_from_Txt(file))
    else:
        print("load_pyCOT_from_file: Wrong extension in input, use .xml or .txt!")
        

def download_webpage(url, filename):
    response = requests.head(url)
    if response.status_code == 200:
        # URL exists, download the webpage
        response = requests.get(url)
        with open(filename, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded {url} to {filename}")
    else:
        print(f"{url} does not exist")


