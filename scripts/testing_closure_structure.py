#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 15:28:44 2023

@author: tveloz
"""

import sys
import os

# Añadir el directorio raíz al PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.reaction_network import ReactionNetwork
from pyCOT.closure_structure import *
from  pyCOT.reactive_features import *
import networkx as nx
from pyCOT.file_manipulation import *
import matplotlib.pyplot as plt
import time
from collections import defaultdict
from collections import Counter # para el gráfico de cantidad de básicos vs repeticiones

# Create an instance of the HelloWorld class
# SpStr= ['a', 'b', 'c', 'd']  # Default: ['a', 'b', 'c', 'd']
# SpBt=bt([True,True,True,True])  # Default: [a, b, c, d]
# RnStr= ['r1', 'r2']  # Updated: ['r1', 'r2']
# RnBt= bt([True, True])  # Default: [r1, r2]
# RnVecS = np.array([[1, 0, 0, 0], [0, 1, 0, 0]])  # Updated for reactions
# RnVecP= np.array([[0, 2, 0, 0], [0, 1, 1, 0]])  # Updated for reactions

# testRN=pyCOT(SpStr,SpBt,RnStr,RnBt,RnVecS,RnVecP)
# print("specs")
# print(testRN.SpStr)
# print("specsBt")
# print(testRN.SpBt)
# print("reacsBt")
# print(testRN.RnBt)
# print("Reac")
# print(testRN.RnStr)
# print("SuppVec")
# print(testRN.RnVecS)
# print("ProdVec")
# print(testRN.RnVecP)
#Example usage:
#file_path = '../../networks/testing/in_out.txt'
# file_path = '../networks/testing/RedPDoSR01.txt'
# file_path = '../networks/testing/RN_AEP_03.txt'
# file_path = '../networks/testing/RN_AEP_04.txt'
# file_path = '../networks/testing/navarino.txt'
#file_path = 'networks/testing/ERCs_test.txt'
# file_path = '../networks/testing/RedPDoSR01.txt'
#file_path = '../networks/biomodels_all/BIOMD0000000011/BIOMD0000000011.xml'
#file_path = 'networks/biomodels_interesting/BIOMD0000000237_manyOrgs.txt'  #ERROR: ZeroDivisionError: division by zero


#file_path = 'networks/biomodels_interesting/BIOMD0000000652_manyOrgs.txt'

#file_path = '../networks/testing/ERCs_test.txt'
#file_path = '../../networks/testing/Synergy_test.txt'
# file_path = '../networks/testing/Farm.txt'
#file_path = '../networks/testing/RedPDoSR01.txt'
#file_path = 'networks/testing/MSORN_test1.txt'
file_path = 'networks/testing/MSORN_test2.txt'
file_path = 'networks/Navarino/RN_IN_02_PyCOT.txt'
file_path = "networks/farm.txt"
#file_path = 'networks/RandomAlife/RN_Ns_40_Norg_20_id_396.txt'
#file_path = 'networks/biomodels_interesting/central_ecoli.txt'
#file_path = '../../networks/biomodels_all/BIOMD0000000086_url.xml'
#file_path = '../../networks/testing/MSORN_test2.txt'
#msorn= load_pyCOT_from_file(file_path)#load_pyCOT_from_Sbml(file_path)
# folder_path='../networks/biomodels_all/'
# subfolders = [f.name for f in os.scandir(folder_path) if f.is_dir()]
# for i in range(40):
#     file_name=subfolders[i]
#     file_path = '../networks/biomodels_all/'
#     file_path=file_path+file_name+"/"+file_name+'.xml'
#     print("Processing path "+file_path)
#     testRN2 = load_pyCOT_from_file(file_path)
#     df_sp=reac_analysis(testRN2)
#     print("analysis done!")
    
#     print(df_sp)
#     array = df_sp.values.astype(float)
#     plt.matshow(array, cmap=plt.cm.viridis)
#     plt.xticks(np.arange(len(df_sp.columns)), df_sp.columns)
#     plt.yticks(np.arange(len(df_sp)), df_sp.index)
#     plt.title(file_name)
#     plt.colorbar()
#     plt.show()
# --------------------------------------------------------
# testRN2 = load_pyCOT_from_file(file_path)
# df_sp=reac_analysis(testRN2)
# #print(df_sp)
# array = df_sp.values.astype(float)
# plt.matshow(array, cmap=plt.cm.viridis)
# plt.xticks(np.arange(len(df_sp.columns)), df_sp.columns)
# plt.yticks(np.arange(len(df_sp)), df_sp.index)
# plt.colorbar()
# plt.show()
# --------------------------------------------------------



# Assuming load_pyCOT_from_file and reac_analysis are defined elsewhere in your code

testRN2 = load_pyCOT_from_file(file_path)
df_sp = reac_analysis(testRN2)

# Uncomment to print the dataframe
print(df_sp)

# Convert the dataframe to a numpy array
array = df_sp.values.astype(float)

# Create a plot with a heatmap using matshow
plt.matshow(array, cmap=plt.cm.viridis)

# Set labels for x and y ticks
plt.xticks(np.arange(len(df_sp.columns)), df_sp.columns)
plt.yticks(np.arange(len(df_sp)), df_sp.index)

# Add colorbar to show the scale
plt.colorbar()

# Add values to the plot and exclude zero values
for i in range(len(df_sp.index)):
    for j in range(len(df_sp.columns)):
        value = array[i, j]
        if value != 0:  # Exclude zero values
            plt.text(j, i, f'{value:.0f}', ha='center', va='center', color='white', fontsize=8, weight='bold')

# Show the plot
plt.show()



erc=ERCs(testRN2)
print()
print("ERCs created")
print()

# --------------------------------------------------------
print("the list of "+str(len(erc))+" ERCs:")
for g in erc:
      print()
      print(g)
print()
# -----------------------------------------------------------

# Open a file in write mode
with open("ERCsDetails.txt", "w") as file:
    # Redirect the standard output to the file
    import sys
    sys.stdout = file
    
    # Your script goes here
    # print("the list of "+str(len(erc))+" ERCs:")
    for g in erc:
        # print()
        print(g)
    print()

# Restore standard output
sys.stdout = sys.__stdout__

                
# Input file path
input_file_path = "ERCsDetails.txt"
# Output file path
output_file_path = "ERCs.txt"

# Open the input file for reading
with open(input_file_path, 'r') as input_file:
    # Open the output file for writing
    with open(output_file_path, 'w') as output_file:
        # Iterate over each line in the input file
        for line in input_file:
            # Strip leading and trailing whitespaces from the line
            line = line.strip()
            # Check if the line is not empty
            if line:
                # Split the line by space characters
                terms = line.split()  # Adjust the delimiter if needed
                # Check if there are terms in the line
                if terms:
                    # Extract the last term using negative indexing
                    last_term = terms[-1]
                    # Remove ']' from the last term if it exists
                    last_term = last_term.rstrip(']')
                    # Add "E" before the last term without space between them
                    modified_term = "E" + last_term
                    # Write the modified term to the output file
                    output_file.write(modified_term + '\n')



print("#Reactions: "+str(len(testRN2.RnStr))+", #ERCs: "+str(len(erc)))
print()
con=get_containments(testRN2,erc)
print("check containment")

# for g in con:
#       print(g)
print("the list of "+str(len(con))+" containments")

dc=get_direct_containments(testRN2,erc,con)
# print("check direct containment")
# for g in dc:
#       print(g)
print("the list of "+str(len(dc))+" direct containments")

# print("check synergies")
syn=get_synergies(testRN2,erc)
# for g in syn:
#       print(g)
print("the list of "+str(len(syn))+" synergies")

d_syn=get_direct_synergies(testRN2,erc,dc,syn)

# for g in d_syn:
#       print(g)
print("the list of "+str(len(d_syn))+" direct synergies")

# Open a file in write mode
with open("Analysis.txt", "w") as file:
    # Write the lengths of con, dc, syn, and d_syn to the file
    file.write(f"{len(con)}\n")
    file.write(f"{len(dc)}\n")
    file.write(f"{len(syn)}\n")
    file.write(f"{len(d_syn)}\n")
    
print()
sorn=second_order_network(testRN2,erc,dc,d_syn)
# print("The list of second order network species")
# for i in range(len(sorn[0][0])):
#       print(str(sorn[0][0][i])+" has lenght "+str(len(sorn[0][1][i])))

print("the list of "+str(len(sorn[1]))+" second order reaction network")
# for i in range(len(sorn[1])):    
#     print("r"+str(i)+": "+str(sorn[1][i]))

print("We have "+str(len(sorn[0][0]))+" second order network species")
# for i in range(len(sorn[0][0])):


##############################################################################


# Assuming MSORN, terminal_species, and reac_analysis are defined elsewhere in your code

msorn = MSORN(sorn)
terminal = terminal_species(msorn)
print("We have " + str(len(terminal)) + " Terminal species")


# Perform reaction analysis
df_sp = reac_analysis(msorn)

# Uncomment to print the dataframe
# print(df_sp)

# Convert the dataframe to a numpy array
array = df_sp.values.astype(float)

# Create a plot with a heatmap using matshow
plt.matshow(array, cmap=plt.cm.viridis)

# Set labels for x and y ticks
plt.xticks(np.arange(len(df_sp.columns)), df_sp.columns)
plt.yticks(np.arange(len(df_sp)), df_sp.index)

# Add colorbar to show the scale
plt.colorbar()

# Add values to the plot and exclude zero values
for i in range(len(df_sp.index)):
    for j in range(len(df_sp.columns)):
        value = array[i, j]
        if value != 0:  # Exclude zero values
            plt.text(j, i, f'{value:.0f}', ha='center', va='center', color='white', fontsize=8, weight='bold')

# Show the plot
plt.show()

print()
start_time = time.time()
cl_pth=[]
for ter in terminal:
    print("computing backward chains of "+ter)
    cl_pth.append(backward_chain_closed(msorn,ter))

end_time = time.time()
execution_time = end_time - start_time
print("Execution time computing closed chains", execution_time, "seconds")
suma=0
for l in cl_pth:
    suma=suma+len(l)
print()
print("number of closed paths = "+str(suma))
file_path = "cl_pth_test.txt"
with open(file_path, 'w') as file:
# Write the string to the file
    for pth in cl_pth:
        for p in pth:
            file.write(str(p)+'\n')

print("Init computing closed backward paths")

start_time = time.time()
cl_pth=[]
cl_pth=Closed_paths_backward(msorn)

end_time = time.time()
execution_time = end_time - start_time
for pth in cl_pth:
    print(str(pth))
print("Execution time computing closed paths", execution_time, "seconds")
print("Total number of closed paths = ", str(len(cl_pth)))

# ----------------------------------------------------------

# Compute the ratios
ratios = [len(con)/len(erc), len(dc)/len(erc), len(syn)/len(erc), len(d_syn)/len(erc)]

# Define labels for each bar
labels = ['containments', 'direct con.', 'synergies', 'direct synergies']

# Define colors for each bar (change as needed)
colors = ['darkblue', 'darkorange', 'darkgreen', 'darkred']

# Plot the bar chart
plt.bar(labels, ratios, color=colors)

# Add values on top of each bar
for i, ratio in enumerate(ratios):
    plt.text(i, ratio, f'{ratio:.2f}', ha='center', va='bottom')

# Add legend
plt.legend(['containments', 'direct containments', 'synergies', 'direct synergies'])

# Set labels and title
plt.xlabel('Ratios of # properties calculated to # ERCs')
plt.ylabel('Ratio')
plt.title('Basic statistical analysis')

# Show plot
plt.show()


# ----------------------------------------------------------

# Define the file path
file_path = "cl_pth_test.txt"

# Read data from the file and count the number of components in each line
component_counts = []
with open(file_path, 'r') as file:
    for line in file:
        components = line.strip().split()
        num_components = len(components)
        component_counts.append(num_components)

# Count the frequency of each component count
count_distribution = Counter(component_counts)

# Extract x and y values for plotting
x_values = list(count_distribution.keys())
y_values = list(count_distribution.values())

color = 'darkblue'

# Plot the distribution
plt.bar(x_values, y_values, align='center', alpha=0.7, color=color)  # Change color here
plt.xlabel('Number of ERCs')
plt.ylabel('Frequency')
plt.title('ERCs Count Distribution')
# plt.grid(True)

# Set integer ticks on both axes
plt.xticks(range(min(x_values), max(x_values) + 1, 1))
plt.yticks(range(min(y_values), max(y_values) + 1, 1))

# Uncomment the following lines to add values to the bars
for i, value in enumerate(y_values):
    plt.text(x_values[i], value + 0.1, str(value), ha='center')

plt.show()


print()


import re
import matplotlib.pyplot as plt

# Step 1: Read data from files
def read_file(filename):
    with open(filename, 'r') as file:
        return [line.strip() for line in file]

ercs = read_file('ERCs.txt')
msorn = read_file('MSORN.txt')

# Step 2: Parse the reactions and perform comparisons
exact_matches = {erc: 0 for erc in ercs}
intersection_matches = {erc: 0 for erc in ercs}

for reaction in msorn:
    match = re.search(r'r\d+:\s*(.*?)\s*=>', reaction)
    if match:
        reactants = match.group(1)
        if '+' in reactants:
            elements = reactants.split('+')
        else:
            elements = [reactants]

        for erc in ercs:
            if erc in elements:
                if len(elements) == 1:
                    exact_matches[erc] += 1
                else:
                    intersection_matches[erc] += 1

# Step 3: Write results to a file
with open('results.txt', 'w') as file:
    for erc in ercs:
        file.write(f"{erc} {exact_matches[erc]} {intersection_matches[erc]}\n")

# Step 4: Plot the results
erc_list = list(exact_matches.keys())
exact_counts = list(exact_matches.values())
intersection_counts = list(intersection_matches.values())

x = range(len(erc_list))

plt.figure(figsize=(10, 5))
plt.bar(x, exact_counts, width=0.4, label='Exact Match', align='center')
plt.bar(x, intersection_counts, width=0.4, label='Intersection Match', align='edge')
plt.xlabel('ERCs')
plt.ylabel('Count')
plt.title('Count of Exact and Intersection Matches')
plt.xticks(x, erc_list, rotation='vertical')
plt.legend()
plt.tight_layout()
plt.savefig('matches_plot.png')
plt.show()







# Closure: The set is closed under the given reaction rules, meaning that any reaction that starts with elements from the set produces elements that are also within the set. In other words, no new species are produced that are not already in the set.

# Self-Maintenance: The set is self-maintaining, meaning that the species within the set can reproduce each other through the given reactions. This implies that each species in the set can be produced from other species within the set.

# These two properties ensure that the set of species can exist over time without external inputs, maintaining its composition through internal reactions.



