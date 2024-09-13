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
from pathlib import Path
import pickle

#file_path = '../../networks/testing/in_out.txt'
# file_path = '../networks/testing/RedPDoSR01.txt'
# file_path = '../networks/testing/RN_AEP_03.txt'
# file_path = '../networks/testing/RN_AEP_04.txt'
#file_path = "networks/Navarino/RN_IN_02_PyCOT.txt"

#file_path = 'networks/biomodels_interesting/bigg_iAF692.txt'

#Definitig the function
# Assuming load_pyCOT_from_file and reac_analysis are defined elsewhere in your code
def Closure_Str_RN_Computing(file_path):
       testRN2 = load_pyCOT_from_file(file_path)
       print("The RN has "+str(len(testRN2.SpStr))+" species and "+str(len(testRN2.RnStr))+" reactions")
       #df=reac_analysis(testRN2)
       #print(df)


       ##########---ERC Caculation---################
       print("Computing ERCs")
       start_time = time.time()
       erc=ERCs(testRN2)
       end_time = time.time()
       execution_time = end_time - start_time
       input_path = Path(file_path)
       
       if input_path.suffix == '.txt':
              output_file = input_path.with_suffix('.pkl')
       else:
              output_file = input_path.with_name(input_path.stem + '.pkl')
       
       #Save the results to the output file using pickle
       with open(output_file, 'wb') as file:
              new_data = {"ERCs": erc}
              pickle.dump(new_data, file)
              new_data = {"ERCs exec time": execution_time}
              pickle.dump(new_data, file)


       # --------------------------------------------------------

       print("#Reactions: "+str(len(testRN2.RnStr))+", #ERCs: "+str(len(erc)))

       ##########---Containments Caculation---################
       print("checking containments")

       start_time = time.time()
       con=get_containments(testRN2,erc)
       end_time = time.time()
       execution_time=end_time-start_time

       # for g in con:
       #       print(g)

       with open(output_file, 'ab') as file:
              new_data = {"Containments": con}
              pickle.dump(new_data, file)
              new_data = {"Containments exec time": execution_time}
              pickle.dump(new_data, file)
       print("We obtained "+str(len(con))+" containments")

       ##########---Direct Containments Caculation---################
       print("Flitering direct containments")
       start_time = time.time()
       dc=get_direct_containments(testRN2,erc,con)
       end_time = time.time()
       execution_time=end_time-start_time
       with open(output_file, 'ab') as file:
              new_data = {"Direct Containments": dc}
              pickle.dump(new_data, file)
              new_data = {"Direct Containments exec time": execution_time}
              pickle.dump(new_data, file)

       print("We obtained "+str(len(dc))+" direct containments")


       ##########---Synergy Caculation---################

       print("Obtaining synergies")
       start_time = time.time()
       syn=get_synergies(testRN2,erc)
       end_time = time.time()
       execution_time=end_time-start_time
       with open(output_file, 'ab') as file:
              new_data = {"Synergies": syn}
              pickle.dump(new_data, file)
              new_data = {"Synergies exec time": execution_time}
              pickle.dump(new_data, file)

       print("We found "+str(len(syn))+" synergies")

       ##########---Direct Synergy Caculation---################
       print("Filtering direct synergies")

       start_time = time.time()
       d_syn=get_direct_synergies(testRN2,erc,dc,syn)

       end_time = time.time()
       execution_time=end_time-start_time
       with open(output_file, 'ab') as file:
              new_data = {"Direct Synergies": d_syn}
              pickle.dump(new_data, file)
              new_data = {"Direct Synergies exec time": execution_time}
              pickle.dump(new_data, file)

       print("We found "+str(len(d_syn))+" direct synergies")
       
       ##########---Second Order RN Caculation---################
       # print("Obtaining SORN")
       # start_time = time.time()
       # sorn=second_order_network(testRN2,erc,dc,d_syn)

       # end_time = time.time()
       # execution_time=end_time-start_time
       # with open(output_file, 'ab') as file:
       #        new_data = {"SORN species": len(sorn[0][0])}
       #        pickle.dump(new_data, file)
       #        new_data = {"SORN reactions": len(sorn[1])}
       #        pickle.dump(new_data, file)
       #        new_data = {"SORN exec time": execution_time}
       #        pickle.dump(new_data, file)
       # print("We found "+str(len(sorn[1]))+" second order reactions")
       # print("We found "+str(len(sorn[0][0]))+" second order network species")
       # for i in range(len(sorn[0][0])):


       ##########--- Minimal Second Order RN and terminal species Caculation---################

       # print("Obtaining MSORN and Terminal species")
       # start_time = time.time()
       # msorn = MSORN(sorn)
       # end_time = time.time()
       # terminal = terminal_species(msorn)
       # execution_time=end_time-start_time
       # with open(output_file, 'ab') as file:
       #        new_data = {"Terminal species": terminal}
       #        pickle.dump(new_data, file)
       #        new_data = {"Terminal species exec time": execution_time}
       #        pickle.dump(new_data, file)

       # print("We have " + str(len(terminal)) + " Terminal species")

       ##########--- Backward chains from terminal species Caculation---################

       # start_time = time.time()
       # cl_pth=[]
       # for ter in terminal:
       #     print("computing backward chains of "+ter)
       #     cl_pth.append(backward_chain_closed(msorn,ter))
       # end_time = time.time()
       # execution_time=end_time-start_time
       # with open(output_file, 'ab') as file:
       #        new_data = {"Number of closed paths": len(cl_pth)}
       #        pickle.dump(new_data, file)
       #        new_data = {"Direct Containments exec time": execution_time}
       #        pickle.dump(new_data, file)

       # print("number of closed paths = "+len(cl_pth))

       ##########--- Closed paths backward from terminal species Caculation---################
       # print("Computing closed paths backwards")
       # start_time = time.time()
       # cl_pth=[]
       # cl_pth=Closed_paths_backward(msorn)

       # end_time = time.time()
       # execution_time = end_time - start_time
       # with open(output_file, 'ab') as file:
       #        new_data = {"Number of Closed Paths": len(cl_pth)}
       #        pickle.dump(new_data, file)
       #        new_data = {"Closed Paths exec time": execution_time}
       #        pickle.dump(new_data, file)
       # print("Execution time computing closed paths backwards", execution_time, "seconds")
       # print("Total number of closed paths = ", str(len(cl_pth)))

       # ----------------------------------------------------------

       ##########--- Closed sets from Closed paths Caculation---################
       # print("Obtaining closed sets from closed paths")

       # start_time = time.time()
       # closed_sets=Closed_paths_to_Sp(cl_pth,sorn)
       # end_time = time.time()
       # execution_time = end_time - start_time
       # with open(output_file, 'ab') as file:
       #        new_data = {"Closed sets": closed_sets}
       #        pickle.dump(new_data, file)
       #        new_data = {"Closed sets exec time": execution_time}
       #        pickle.dump(new_data, file)
       # print("Execution time computing closed paths", execution_time, "seconds")
       # print("Total number of closed sets = ", str(len(closed_sets)))

       ##########--- Semi-self-maintaining sets from Closed sets Caculation---################

       # print("Filtering semi-self-maintaining sets from closed sets")
       # start_time = time.time()
       # ssm_sets=[ssm_set for ssm_set in closed_sets if testRN2.is_semi_self_maintaining (ssm_set)]
       # end_time = time.time()
       # execution_time = end_time - start_time
       # with open(output_file, 'ab') as file:
       #        new_data = {"Self-maintaining sets": ssm_sets}
       #        pickle.dump(new_data, file)
       #        new_data = {"Self-maintaining sets exec time": execution_time}
       #        pickle.dump(new_data, file)
       # print("Execution time com3puting ssm paths", execution_time, "seconds")
       # print("Total number of ssm sets = ", str(len(ssm_sets)))

#Closure_Str_RN_Computing(file_path)
#results={"ERC":erc,"Containments": con,"Direct containments":dc,
#         "Synergies":syn,"Direct Synergies": d_syn,"Sorn": sorn,"Minimal sorn": msorn, 
#         "Execution time": execution_time,"Closed paths": cl_pth,
#         "Closed sets":closed_sets,"Semi-self-maintaining sets": ssm_sets}
#input_path = Path(file_path)
    
#if input_path.suffix == '.txt':
#        output_file = input_path.with_suffix('.pkl')
#else:
#        output_file = input_path.with_name(input_path.stem + '.pkl')
    
    # Save the results to the output file using pickle
#with open(output_file, 'ab') as file:
#        pickle.dump(results, file)
#tbp="Execution time"
#print(tbp)
#print(results[tbp])

# print("Number of species ="+str(len(testRN2.SpStr)))
# print("Number of reactions ="+str(len(testRN2.RnStr)))
# print("Number of ERCs ="+str(len(msorn.SpStr)))
# print("SSM = "+str(SSM)+", not SSM ="+str(nSSM))
# print("relative size SSM ="+str(sizeSSM)+"/"+str(SSM))
# print("relative size nSSM ="+str(sizenSSM)+"/"+str(nSSM))

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

# --------------------------------------------------------
