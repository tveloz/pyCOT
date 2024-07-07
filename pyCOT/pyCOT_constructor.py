from collections.abc import Sequence

import numpy as np
from bitarray import bitarray as bt
import networkx as nx
# import pandas as pd

import py_data_conversion as dc

# Created by Tomas Veloz - github.com/tveloz
############################################################################################################
# Summary: Create/Load/Visualize reaction networks and compute COT-related properties#######################
############################################################################################################


#Table of contents:#########################################################################################

#Constructor################################################################################################
##*Vector of positions to bt transformations################################################################
##*Bt to/from species list transformations #################################################################
##*obtaining sets of species/reactions producing/consuming one another######################################
##*Getting various notions of species connectivity##########################################################
##*verify relational properties of species/reactions########################################################


class pyCOT:
    """
    Class representing pyCOT (Python Chemical Organization Theory Object)

    Attributes:
    - SpBt: Bitarray identification for species
    - SpStr: List of strings (species names) identification
    - RnBt: Bitarray identification for reactions
    - RnMsupp: Vector (numpy.array) identification support of reactions 
    - RnMprod: Vector (numpy.array) identification for products of reactions
    - RnStr: List of strings (reaction names) identification

    Methods (not to be updated just read below):
    - __init__: Constructor method to initialize the class with the provided parameters.
    - get_id_from_bt: Function that returns a vector from bitarray representation.
    - set_bt_from_id: Function that returns bitarray from vector representation.
    """

    def __init__(self, SpStr: Sequence[str], SpBt: bt, RnStr: Sequence[str], RnBt: bt,
                 RnMsupp: np.ndarray, RnMprod: np.ndarray):
        """
        Constructor for pyCOT class.

        Parameters:
        - SpBt: Bitarray identification for species
        - SpStr: Sequence of strings (species names) identification
        - RnBt: Bitarray identification for reactions
        - RnMsupp: Matrix (numpy.darray) identification support of reactions 
        - RnMprod: Matrix (numpy.darray) identification for products of reactions
        - RnStr: Sequence of strings (reaction names) identification
        """
        # Species
        self.SpStr = SpStr
        self.SpBt = SpBt

        # Reactions
        self.RnStr = RnStr
        self.RnBt = RnBt
        self.RnMsupp = RnMsupp
        self.RnMprod = RnMprod



  ##########################################################################################################
  #############bt From/To str representations###############################################################
  ##########################################################################################################

    def get_bt_from_species(self, SpStr):
        """Function that returns bitarray from a List of strings (species names) identification."""
        if isinstance(SpStr, str):
            SpStr = [SpStr]
        for i in range(len(SpStr)):
            if SpStr[i] not in self.SpStr:
                print("get_bt_from_species ERROR: input is not a list of recognized species")
                return []
        bitarray = bt(len(self.SpStr))
        bitarray.setall(0)
        for i in range(len(SpStr)):
            for j in range(len(bitarray)):
                if SpStr[i] == self.SpStr[j]:
                    bitarray[j] = True
        return bitarray
    
    def get_bt_from_reactions(self, RnStr):
        """Function that returns bitarray from a List of strings (reaction names) identification."""
        if isinstance(RnStr, str):
            RnStr = [RnStr]
        for i in range(len(RnStr)):
            if RnStr[i] not in self.RnStr:
                print("get_bt_from_reactions ERROR: input is not a list of recognized species")
                return []
        bitarray = bt(len(self.RnStr))
        bitarray.setall(0)
        for i in range(len(RnStr)):
            for j in range(len(self.RnStr)):
                if RnStr[i] == self.RnStr[j]:
                    bitarray[j] = True
        return bitarray

    def get_species_from_bt(self, bitarray):
        """Function that returns List of strings (species names) identification from bitarray.""" 
        species_list = self.SpStr
        if not isinstance(bitarray, bt):
            print("get_species_from_bt ERROR: input is not a bitarray ")
            return None
        if len(species_list) != len(bitarray):
            print("get_species_from_bt ERROR: bitarray input has different length than species set size, can't continue")
        else:
            selected_species = []
            for i in range(len(bitarray)):
                if bitarray[i]:
                    selected_species.append(species_list[i])
            return selected_species

    def get_reactions_from_bt(self, bitarray):
        """Function that returns List of strings (reactions names) identification from bitarray."""
        selected_reactions = []
        if not isinstance(bitarray, bt):
            print("get_reactions_from_bt ERROR: input is not a bitarray ")
            return None
        if len(self.RnStr) != len(bitarray):
            print("get_reactions_from_bt ERROR: bitarray input has different length than reactions set size, can't continue")
            return None
        else:
            for i in range(len(bitarray)):
                if bitarray[i]:
                    selected_reactions.append(self.RnStr[i])
            return selected_reactions


    def get_species_abstraction_from_vector(self, vec, t):
        """Function that returns a the list of species string with value larger than t in a vector"""
        bitarray = dc.get_bt_abstraction_from_vector(vec, t)
        return(self.get_species_from_bt(bitarray))

    def get_supp_bt_from_reaction(self, reaction_name, t=0):
        """Function that returns bitarray (reaction's supports)  from a string (reaction names) identification ¿in t=0?""" #G: To check description about "t"
        if reaction_name not in self.RnStr:
            print("get_supp_bt_from_reaction ERROR: Reaction '{reaction_name}' not found in the reactions set.")
            return None
        else:
            reaction_index = self.RnStr.index(reaction_name)
            support_vec = self.RnMsupp[reaction_index]
            support_bitarray = dc.get_bt_abstraction_from_vector(support_vec, t)
        return support_bitarray

    def get_prod_bt_from_reaction(self,reaction_name,t=0):
        """Function that returns bitarray (reaction's products) from a string (reaction names) identification ¿in t=0?""" #G: To check description about "t"
        if reaction_name not in self.RnStr:
            print("get_prod_bt_from_reaction ERROR: Reaction '{reaction_name}' not found in the reactions set.")
            return None
        else:
            reaction_index = self.RnStr.index(reaction_name)
            product_vec = self.RnMprod[reaction_index]
            product_bitarray = dc.get_bt_abstraction_from_vector(
                product_vec, t)
        return product_bitarray

    ########################################################################################################
    #### obtaining sets of species/reactions producing/consuming one another################################
    ########################################################################################################              
    
    def get_reactions_from_species(self, SpStr,t=0):
    # Get the bitarray for the given set of species
        if not isinstance(SpStr, list):
            SpStr = [SpStr]
        for i in range(len(SpStr)):
            if SpStr[i] not in self.SpStr:
                print("get_reactions_from_species ERROR: input" +
                      SpStr[i]+" is not a list of recognized species")
                return []

        species_bitarray = self.get_bt_from_species(SpStr)
    # Initialize an empty bitarray for triggered reactions
        triggered_reactions_bitarray = bt(len(self.RnStr))
        triggered_reactions_bitarray.setall(0)
    # Iterate through reactants and check if the species can trigger the reaction
        for i in range(len(self.RnStr)):
            supp = self.RnMsupp[i]
            supp_bt = self.get_bt_abstraction_from_vector(supp, t)
            # checking if supp_bt is contained in species_bitarray
            if (supp_bt & species_bitarray) == supp_bt:
                triggered_reactions_bitarray[i] = True
            else:
                triggered_reactions_bitarray[i] = False
        return self.get_reactions_from_bt(triggered_reactions_bitarray)

    def get_supp_from_reactions(self, RnStr):
        if not isinstance(RnStr, list):
            RnStr = [RnStr]
        for i in range(len(RnStr)):
            if RnStr[i] not in self.RnStr:
                print("ERROR in get_supp_from_reactions: input" +
                      RnStr[i]+" is not a list of recognized reactions")
                return []

        reactions_bitarray = self.get_bt_from_reactions(RnStr)
        specs = bt(len(self.SpStr))
        specs.setall(0)
        for i in range(len(self.RnStr)):
            if reactions_bitarray[i]:
                supp = self.get_supp_bt_from_reaction(self.RnStr[i])
                specs = specs | supp
        return self.get_species_from_bt(specs)

    def get_prod_from_reactions(self, RnStr):
        if not isinstance(RnStr, list):
            RnStr = [RnStr]
        for i in range(len(RnStr)):
            if RnStr[i] not in self.RnStr:
                print("ERROR in get_prod_from_reactions: input" +
                      RnStr[i]+" is not a list of recognized reactions")
                return []
        reactions_bitarray = self.get_bt_from_reactions(RnStr)
        #print("get_bt_from_reactions "+str(reactions_bitarray))
        specs = bt(len(self.SpStr))
        specs.setall(0)
        for i in range(len(self.RnStr)):
            if reactions_bitarray[i]:
                prod = self.get_prod_bt_from_reaction(self.RnStr[i])
                #print("bit "+str(i)+" adds "+str(prod))
                specs = specs | prod
                #print("specs updated "+str(specs))
        return self.get_species_from_bt(specs)

    def get_species_from_reactions(self, RnStr):
        if not isinstance(RnStr, list):
            RnStr = [RnStr]
        for i in range(len(RnStr)):
            if RnStr[i] not in self.RnStr:
                print("get_prod_from_reactions ERROR: input" +
                      RnStr[i]+" is not a list of recognized reactions")
                return []
        prod = self.get_prod_from_reactions(RnStr)
        supp = self.get_supp_from_reactions(RnStr)
        return(list(set(prod).union(set(supp))))

    def get_prod_from_species(self, SpStr):
        if not isinstance(SpStr, list):
            SpStr = [SpStr]
        for i in range(len(SpStr)):
            if SpStr[i] not in self.SpStr:
                print("get_prod_from_species ERROR: input" +
                      SpStr[i]+" is not a list of recognized species")
                return []
        reactions = self.get_reactions_from_species(SpStr)
        prod = self.get_prod_from_reactions(reactions)
        return prod

    def get_reactions_consuming_species(self, SpStr):
        if not isinstance(SpStr, list):
            SpStr = [SpStr]
        for i in range(len(SpStr)):
            if SpStr[i] not in self.SpStr:
                print("get_reactions_consuming_species ERROR: input" +
                      SpStr[i]+" is not a list of recognized species")
                return []

        reactions_list = []
        for i in range(len(self.RnStr)):
            r_supp = self.get_supp_from_reactions(self.RnStr[i])
            if len(set(SpStr).intersection(set(r_supp))) > 0 | (len(r_supp) == 0):
                reactions_list.append(self.RnStr[i])
        return(reactions_list)

    def get_reactions_producing_species(self, SpStr):
        if not isinstance(SpStr, list):
            SpStr = [SpStr]
        for i in range(len(SpStr)):
            if SpStr[i] not in self.SpStr:
                print("get_reactions_producing_species ERROR: input" +
                      SpStr[i]+" is not a list of recognized species")
                return []

        reactions_list = []
        for i in range(len(self.RnStr)):
            r_prod = self.get_prod_from_reactions(self.RnStr[i])
         #   print(str(self.RnStr[i])+" produces "+str(r_prod))
            if len(set(SpStr).intersection(set(r_prod))) > 0:
              #      print(str(self.RnStr[i])+" intersection found")
                reactions_list.append(self.RnStr[i])
        return(reactions_list)

    # ######################################################################################################
    # ######################################################################################################
    #********************** IN PROGRESS ******************************######################################
    # #### obtaining vector of support/prod from reactions/species##########################################
    # ######################################################################################################

    # def get_indexes_for_reactions(self,RnStr):
    #     reactions=self.RnStr
    #     rn_vec=[]
    #     for i in range(len(reactions)):
    #         if reactions[i] in RnStr:
    #             rn_vec.append(i)
    #     return rn_vec             

    # def get_indexes_for_species(self,SpStr):
    #     species=self.SpStr
    #     sp_vec=[]
    #     for i in range(len(species)):
    #         if species[i] in SpStr:
    #             sp_vec.append(i)
    #     return sp_vec

    # def get_supp_matrix_from_reactions(self, RnStr):
    #     supp_Ms = self.RnMsupp
    #     sp_active = self.get_species_from_reactions(RnStr)
    #     sub_Msupp = np.zeros((len(RnStr), len(sp_active)))
    #     reactions_index_list = self.get_indexes_for_reactions(RnStr)
    #     species_index_list = self.get_indexes_for_species(sp_active)
    #     # print("active reactions")
    #     # for r in reactions_index_list:
    #     #     print(self.RnStr[r])
    #     # print("active species")
    #     # for s in species_index_list:
    #     #     print(self.SpStr[s])
    #     for i, ri in enumerate(reactions_index_list):
    #         for j, sj in enumerate(species_index_list):
    #             sub_Msupp[i, j] = supp_Ms[ri][sj]
    #     return sub_Msupp

    # def get_prod_matrix_from_reactions(self, RnStr):
    #     supp_Mp = self.RnMprod
    #     sp_active = self.get_species_from_reactions(RnStr)
    #     sub_Mprod = np.zeros((len(RnStr), len(sp_active)))
    #     reactions_index_list = self.get_indexes_for_reactions(RnStr)
    #     species_index_list = self.get_indexes_for_species(sp_active)
    #     # print("active reactions")
    #     # for r in reactions_index_list:
    #     #     print(self.RnStr[r])
    #     # print("active species")
    #     # for s in species_index_list:
    #     #     print(self.SpStr[s])
    #     for i, ri in enumerate(reactions_index_list):
    #         for j, sj in enumerate(species_index_list):
    #             sub_Mprod[i, j] = supp_Mp[ri][sj]
    #     return sub_Mprod

    # #######Sub and Super element constructor####################
    # # def sub_pyCOT(self,SpStr=SpStr,RnStr=None):
    # #     if RnStr=None:
    # #         RnStr=

    ########################################################################################################
    #############Getting various notions of species connectivity############################################
    ########################################################################################################

    def get_connected_species_to_species(self, SpStr):
        if not isinstance(SpStr, list):
            SpStr = [SpStr]
        for i in range(len(SpStr)):
            if SpStr[i] not in self.SpStr:
                print("get_connected_species_to_species ERROR: input" +
                      SpStr[i]+" is not a list of recognized species")
                return []

        new = list(set(SpStr).union(set(self.get_inflow())))
        result = []
        while len(new) != 0:
            print("iter " + str(new))
            supp = []
            prod = []
            reacs = []
            result = list(set(result).union(set(new)))
            for i in range(len(new)):
                r_prod = self.get_reactions_producing_species(new[i])
                print("adding prod"+str(r_prod))
                r_supp = self.get_reactions_consuming_species(new[i])
                print("adding supp"+str(r_supp))
                reacs = list(set(reacs).union(set(r_supp)))
                reacs = list(set(reacs).union(set(r_prod)))
            supp = self.get_supp_from_reactions(reacs)
            prod = self.get_prod_from_reactions(reacs)
            new = list(set(supp).union(set(new)))
            new = list(set(prod).union(set(new)))
            new = list(set(new)-set(result))
        return result

    def get_immediately_connected_species_to_species(self, SpStr):
        if not isinstance(SpStr, list):
            new = [SpStr]
        else:
            new = SpStr.copy()
        for i in range(len(new)):
            if new[i] not in self.SpStr:
                print("get_immediately_connected_species_to_species ERROR: input" + 
                      new[i]+" is not a list of recognized species")
                return []

        supp = []
        prod = []
        reacs = []
        for i in range(len(new)):
            r_prod = self.get_reactions_producing_species(new[i])
            r_supp = self.get_reactions_consuming_species(new[i])
            reacs = list(set(reacs).union(set(r_supp)))
            reacs = list(set(reacs).union(set(r_prod)))
        supp = self.get_supp_from_reactions(reacs)
        prod = self.get_prod_from_reactions(reacs)
        new = list(set(supp).union(set(new)))
        new = list(set(prod).union(set(new)))
        return(new)

    def get_forward_connected_species_to_species(self, SpStr):
        if not isinstance(SpStr, list):
            new = [SpStr]
        else:
            new = SpStr.copy()
        for i in range(len(new)):
            if new[i] not in self.SpStr:
                print("get_forward_connected_species_to_species ERROR: input" +
                      new[i]+" is not a list of recognized species")
                return []

        result = []
        #print("initial list "+str(new))
        while len(new) != 0:
            prod = []
            reacs = []
            result = list(set(result).union(set(new)))
            for i in range(len(new)):
                r_supp = self.get_reactions_consuming_species(new[i])
                reacs = list(set(reacs).union(set(r_supp)))
            prod = self.get_prod_from_reactions(reacs)
            new = list(set(prod).union(set(new)))
            new = list(set(new)-set(result))
        return(result)

    def get_immediately_forward_connected_species_to_species(self, SpStr):
        if not isinstance(SpStr, list):
            new = [SpStr]
        else:
            new = SpStr.copy()
        for i in range(len(new)):
            if new[i] not in self.SpStr:
                print("get_immediately_forward_connected_species_to_species ERROR: input" +
                      new[i]+" is not a list of recognized species")
                return []

        prod = []
        reacs = []
        for i in range(len(new)):
            r_supp = self.get_reactions_consuming_species(new[i])
            reacs = list(set(reacs).union(set(r_supp)))
        prod = self.get_prod_from_reactions(reacs)
        new = list(set(prod).union(set(new)))
        return(new)

    def get_backward_connected_species_to_species(self, SpStr):
        if not isinstance(SpStr, list):
            new = [SpStr].copy()
        else:
            new = SpStr.copy()
        for i in range(len(new)):
            if new[i] not in self.SpStr:
                print("get_backward_connected_species_to_species ERROR: input" +
                      new[i]+" is not a list of recognized species")
                return []

        result = []
        while len(new) != 0:
            supp = []
            reacs = []
            result = list(set(result).union(set(new)))
            for i in range(len(new)):
                r_prod = self.get_reactions_producing_species(new[i])
                reacs = list(set(reacs).union(set(r_prod)))
            supp = self.get_supp_from_reactions(reacs)
            new = list(set(supp).union(set(new)))
            new = list(set(new)-set(result))
        return(result)

    def get_immediately_backward_connected_species_to_species(self, SpStr):
        if not isinstance(SpStr, list):
            new=[SpStr]
        else:
            new=SpStr.copy()
        for i in range(len(new)):
            if new[i] not in self.SpStr:
                print("get_immediately_backward_connected_species_to_species ERROR: input" + 
                       new[i]+" is not a list of recognized species")
                return []

        supp=[]
        reacs=[]
        for i in range(len(new)):
            r_prod=self.get_reactions_producing_species(new[i])
            reacs=list(set(reacs).union(set(r_prod)))
        supp=self.get_supp_from_reactions(reacs)
        new=list(set(supp).union(set(new)))
        return new

    def get_immediately_strictly_backward_connected_species_to_species(self, SpStr):
        if not isinstance(SpStr, list):
            new = [SpStr]
        else:
            new = SpStr.copy()
        for i in range(len(new)):
            if new[i] not in self.SpStr: # It said Sptr
                print("get_immediately_strictly_backward_connected_species_to_species ERROR: input" +
                      new[i]+" is not a list of recognized species")
                return []

        supp = []
        reacs = []
        for i in range(len(new)):
            r_prod = self.get_reactions_producing_species(new[i])
            reacs = list(set(reacs).union(set(r_prod)))
        supp = self.get_supp_from_reactions(reacs)
        return supp
          
    ########################################################################################################
    ############# get inflow and outflow####################################################################
    ########################################################################################################

    def get_inflow(self):
        return(self.get_prod_from_species([]))

    def get_outflow(self):
        result = bt(len(self.SpBt))
        result.setall(0)
        reacs = self.RnStr
        for i in range(len(reacs)):
            supp = self.get_supp_bt_from_reaction(reacs[i])
            prod = self.get_prod_bt_from_reaction(reacs[i])
            if not prod.any():
                result = result | supp
        return(self.get_species_from_bt(result))

    ########################################################################################################
    ############# verify relational properties of species/reactions#########################################
    ########################################################################################################

    def is_closed(self, SpStr):
        species_bitarray = self.get_bt_from_species(SpStr)
        reactions_list = self.get_reactions_from_species(SpStr)
        prod_of_reactions = self.get_prod_from_reactions(reactions_list)
        prod_bitarray = self.get_bt_from_species(prod_of_reactions)
        return (prod_bitarray | species_bitarray) == species_bitarray

    def is_semi_self_maintaining(self, SpStr):
        reactions_list = self.get_reactions_from_species(SpStr)
        prod_of_reactions = self.get_prod_from_reactions(reactions_list)
        prod_bitarray = self.get_bt_from_species(prod_of_reactions)
        supp_of_reactions = self.get_supp_from_reactions(reactions_list)
        supp_bitarray = self.get_bt_from_species(supp_of_reactions)
        return (supp_bitarray & prod_bitarray) == supp_bitarray
     
    def is_connected(self,SpStr):
        if not isinstance(SpStr,list):
            SpStr=[SpStr]
        if len(SpStr)==0:
            return True
        else:
            connected=self.get_connected_species_to_species(SpStr[0])
            if set(SpStr).issubset(set(connected)):
                print("True becasue conn= "+str(connected))
                return True
            else:
                print("False becasue conn= "+str(connected))
                return False
              

 ###########################################################################################################
 ############# Transform and operate the pyCOT object as a bipartite graph in networkx #####################
 ###########################################################################################################

    def pyCOT_to_Graph(self, SpStr=None, RnStr=None):
        if SpStr == None:
            SpStr = self.SpStr
        if RnStr == None:
            RnStr = self.RnStr
        
        G = nx.DiGraph()

        for s in SpStr:
            G.add_node(s, bipartite=1, type="species")
        for r in RnStr:
            G.add_node(r, bipartite=0, type="reaction")
        for j in range(len(self.RnMsupp)):
            reaction = self.RnMsupp[j]
            for i in range(len(self.SpStr)):
                if reaction[i] > 0:
                    G.add_edge(self.SpStr[i], self.RnStr[j])
        for j in range(len(self.RnMprod)):
            reaction = self.RnMprod[j]
            for i in range(len(self.SpStr)):
                if reaction[i] > 0:
                    G.add_edge(self.RnStr[j], self.SpStr[i])
        return G
