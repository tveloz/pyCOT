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
import xml.etree.ElementTree as ET
from collections import Counter

from pyCOT.reaction_network import *


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

###############   TURN SBML FILE INTO TXT ###########################################3
#Input: An SBML XML file containing species and reaction definitions.
#Output: A text file for each XML file, containing formatted biochemical reactions in a human-readable format.
#Goal: To convert the SBML file into a more digestible format by translating complex species and reactions into a simpler representation.

def extract_reactions_from_file(xml_file, output_folder, log_file):
    # Obtener el nombre del archivo de salida .txt basado en el nombre del archivo XML
    base_filename = os.path.basename(xml_file).replace('.xml', '.txt')
    output_file = os.path.join(output_folder, base_filename)

    # Cargar el archivo XML
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Extraer el espacio de nombres automáticamente del root
    namespace = {'sbml': root.tag.split('}')[0].strip('{')} if '}' in root.tag else {}

    # Diccionario para almacenar el id y nombre de las especies
    species_dict = {}

    # Lista para almacenar las reacciones formateadas
    reaction_list = []

    # Extraer especies (id -> name) para asignar nombres a los reactantes y productos
    for species in root.findall('.//sbml:listOfSpecies/sbml:species', namespace):
        species_id = species.get('id')
        species_name = species.get('name', species_id)  # Usa species_id como fallback
        # Evitar confusión entre especies y eliminar cualquier ";" en los nombres
        if species_id and species_name:
            # Reemplazamos cualquier ";" en el nombre por "_"
            clean_species_name = species_name.replace(";", "_").replace(":", "_")
            unique_species_name = f"{species_id} ({clean_species_name})"
            species_dict[species_id] = unique_species_name

    # Extraer reacciones y asignar nombres "R0:", "R1:", etc.
    for idx, reaction in enumerate(root.findall('.//sbml:listOfReactions/sbml:reaction', namespace)):
        # Asignar el nombre de la reacción como "R0:", "R1:", etc.
        reaction_name = f"R{idx}:"

        # Obtener y contar los reactantes
        reactants = []
        for reactant in reaction.findall('.//sbml:listOfReactants/sbml:speciesReference', namespace):
            reactant_id = reactant.get('species')  # ID del reactante
            stoichiometry = float(reactant.get('stoichiometry', 1))  # Obtener el coeficiente estequiométrico (por defecto 1)
            if reactant_id and reactant_id in species_dict:
                reactants.append((species_dict[reactant_id], stoichiometry))

        # Obtener y contar los productos
        products = []
        for product in reaction.findall('.//sbml:listOfProducts/sbml:speciesReference', namespace):
            product_id = product.get('species')  # ID del producto
            stoichiometry = float(product.get('stoichiometry', 1))  # Obtener el coeficiente estequiométrico (por defecto 1)
            if product_id and product_id in species_dict:
                products.append((species_dict[product_id], stoichiometry))

        # Extraer y contar los modificadores, incluyendo su estequiometría
        modifiers = []
        for modifier in reaction.findall('.//sbml:listOfModifiers/sbml:modifierSpeciesReference', namespace):
            modifier_id = modifier.get('species')
            stoichiometry = float(modifier.get('stoichiometry', 1))  # Permitir obtener la estequiometría si está presente
            if modifier_id and modifier_id in species_dict:
                modifiers.append((species_dict[modifier_id], stoichiometry))

        # Sumamos modificadores a reactantes y productos, contabilizando correctamente las cantidades
        reactants_counter = Counter()
        products_counter = Counter()

        # Sumar reactantes con su estequiometría
        for reactant, stoichiometry in reactants:
            reactants_counter[reactant] += stoichiometry
        
        # Sumar productos con su estequiometría
        for product, stoichiometry in products:
            products_counter[product] += stoichiometry

        # Modificadores se suman tanto a reactantes como a productos con su propia estequiometría
        for modifier, stoichiometry in modifiers:
            reactants_counter[modifier] += stoichiometry
            products_counter[modifier] += stoichiometry

        # Verificar si hay al menos un reactante o producto para la reacción
        if len(reactants_counter) == 0 and len(products_counter) == 0:
            # Si no hay reactantes ni productos, marcar como reacción "vacía" o "desconocida"
            reactants_str = "No reactantes"
            products_str = "No productos"
        else:
            # Crear la representación de reactantes y productos con multiplicidades, usando floats y reemplazando el ";"
            reactants_str = " + ".join([f"{count:.2f}{species}" if count > 1 else species for species, count in reactants_counter.items()])
            products_str = " + ".join([f"{count:.2f}{species}" if count > 1 else species for species, count in products_counter.items()])

        # Crear la cadena en el formato "R0: reactantes => productos"
        reaction_str = f"{reaction_name} {reactants_str} => {products_str}"
        reaction_list.append(reaction_str)

    # Escribir las reacciones en un archivo .txt
    with open(output_file, 'w') as f:
        for reaction in reaction_list:
            f.write(reaction + '\n')

    # Mensaje para consola y archivo de log
    log_message = f"Procesado {xml_file} - {len(reaction_list)} reacciones guardadas en {output_file}.\n"
    print(log_message)
    
def process_single_file(xml_file):
    # Crear la carpeta de salida con el nombre del archivo de entrada sin la extensión + "_txt"
    base_name = os.path.splitext(os.path.basename(xml_file))[0]
    output_folder = base_name + "_txt"
    os.makedirs(output_folder, exist_ok=True)

       # Procesar el archivo XML
    extract_reactions_from_file(xml_file, output_folder, log_file)

    print(f"El archivo .xml ha sido procesado. Los resultados están en la carpeta: {output_folder}")
def load_pyCOT_from_file(file):
    input_file_name = file
    
    base_name, extension = os.path.splitext(input_file_name)
    if extension=='.xml':
        print("please convert file into txt using function extract_reactions_from_file")
        return(0)
    elif extension=='.txt':
        return(load_pyCOT_from_Txt(file))
    else:
        print("load_pyCOT_from_file: Wrong extension in input, use .xml or .txt!")


# Ejemplo de uso
#xml_file = "BIOMD0000000110.xml"  # Reemplaza con la ruta a tu archivo .xml
#process_single_file(xml_file)

############### DEPRECATED CODE TO WORK WITH SBML FILES #######################
# def load_pyCOT_from_Sbml(file):
    
#     # Input XML file name
#     input_file_name = file
#     base_name, extension = os.path.splitext(input_file_name)
#     if extension!='.xml':
#         print("load_pyCOT_from_Sbml: Wrong extension in input!")
#         return 0
#     # Open the XML file and read its contents
#     with open(input_file_name, "r") as file:
#         xml_content = file.read()
    
#     # Parse the XML content using BeautifulSoup with the "xml" parser
#     soup = bs(xml_content, "xml")
    
#     # Extract file name without extension
#     file_name_without_extension = os.path.splitext(input_file_name)[0]
    
#     # Output file name with .txt extension
#     output_file_name = file_name_without_extension + ".txt"
    
#     # Extract all <reaction> tags
#     reaction_tags = soup.find_all("reaction")
    
#     # Create a new text file for writing
#     with open(output_file_name, "w") as output_file:
#         # Process each <reaction> tag
#         for reaction_tag in reaction_tags:
#             # Extract the reaction id and reversible attribute
#             reaction_id = reaction_tag["id"]
#             reversible = reaction_tag.get("reversible", "false") == "true"
    
#             # Extract reactants, products, and modifiers for the reaction if <listOfReactants>, <listOfProducts>, and <listOfModifiers> tags exist
#             reactants_tag = reaction_tag.find("listOfReactants")
#             products_tag = reaction_tag.find("listOfProducts")
#             modifiers_tag = reaction_tag.find("listOfModifiers")
            
#             # Extract reactants, products, and modifiers with stoichiometry if the tags exist, otherwise set to empty lists with stoichiometry 1
#             if reactants_tag:
#                 reactants = [(reactant["species"], int(reactant.get("stoichiometry", "1"))) for reactant in reactants_tag.find_all("speciesReference")]
#             else:
#                 reactants = []
    
#             if products_tag:
#                 products = [(product["species"], int(product.get("stoichiometry", "1"))) for product in products_tag.find_all("speciesReference")]
#             else:
#                 products = []
    
#             if modifiers_tag:
#                 modifiers = [(modifier["species"], int(modifier.get("stoichiometry", "1"))) for modifier in modifiers_tag.find_all("modifierSpeciesReference")]
#             else:
#                 modifiers = []
    
#             # Filter out empty categories
#             reactants_str = " + ".join([f"{stoich}{species}" for species, stoich in reactants if stoich != 0])
#             products_str = " + ".join([f"{stoich}{species}" for species, stoich in products if stoich != 0])
#             modifiers_str = " + ".join([f"{stoich}{species}" for species, stoich in modifiers if stoich != 0])
#             # Write the reaction id and the reaction as a sum of stoichiometry*species => sum of stoichiometry*species with modifiers on both sides of the reaction, ending with a semicolon
#             reaction_line=reaction_id + ": " + (modifiers_str + " + " if modifiers_str else "") + reactants_str +  " => "+ (modifiers_str + " + "  if modifiers_str else "") + products_str + ";\n"
#             reaction_line=reaction_line.replace("+  =>","=>")
#             output_file.write(reaction_line)
            
#             # Check if the reaction is reversible
#             if reversible:
#                 # Write the backward reaction (products as reactants, reactants as products) with modifiers on both sides of the reaction
#                 output_file.write(reaction_id + "_b: " + (modifiers_str + " + " if modifiers_str else "") + products_str + " + " +  " => " + ( modifiers_str + " + " if modifiers_str else "") + reactants_str + ";\n")

#     return(load_pyCOT_from_file(output_file_name)) 
        
# def download_webpage(url, filename):
#     response = requests.head(url)
#     if response.status_code == 200:
#         # URL exists, download the webpage
#         response = requests.get(url)
#         with open(filename, 'wb') as f:
#             f.write(response.content)
#         print(f"Downloaded {url} to {filename}")
#     else:
#         print(f"{url} does not exist")1
