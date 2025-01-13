from collections.abc import Sequence

import numpy as np
from bitarray import bitarray as bt
import networkx as nx
import itertools
# import pandas as pd

import pyCOT.data_conversion as dc

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


class ReactionNetwork:
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

        species_bitarray = dc.get_bt_from_names(SpStr, self.SpStr)
    # Initialize an empty bitarray for triggered reactions
        triggered_reactions_bitarray = bt(len(self.RnStr))
        triggered_reactions_bitarray.setall(0)
    # Iterate through reactants and check if the species can trigger the reaction
        for i in range(len(self.RnStr)):
            supp = self.RnMsupp[i]
            supp_bt = dc.get_bt_abstraction_from_vector(supp, t)
            # checking if supp_bt is contained in species_bitarray
            if (supp_bt & species_bitarray) == supp_bt:
                triggered_reactions_bitarray[i] = True
            else:
                triggered_reactions_bitarray[i] = False
        return dc.get_reactions_from_bt(triggered_reactions_bitarray, self.RnStr)

    def get_supp_from_reactions(self, RnStr):
        #Input: list of reactions
        #Output:  List of reactants of those reactions
        if not isinstance(RnStr, list):
            RnStr = [RnStr]
        reactions_bitarray = dc.get_bt_from_names(RnStr, self.RnStr)
        specs = bt(len(self.SpStr))
        specs.setall(0)
        for i in range(len(self.RnStr)):
            if reactions_bitarray[i]:
                supp = dc.get_species_bt_from_reaction(self.RnStr[i], self.RnStr, self.RnMsupp)
                specs = specs | supp
        return dc.get_species_from_bt(specs, self.SpStr)

    def get_prod_from_reactions(self, RnStr):
        #Input: List of reactions
        #Output: List of products of those reactions
        if not isinstance(RnStr, list):
            RnStr = [RnStr]
        for i in range(len(RnStr)):
            if RnStr[i] not in self.RnStr:
                print("ERROR in get_prod_from_reactions: input" +
                      RnStr[i]+" is not a list of recognized reactions")
                return []
        reactions_bitarray = dc.get_bt_from_names(RnStr, self.RnStr)
        specs = bt(len(self.SpStr))
        specs.setall(0)
        for i in range(len(self.RnStr)):
            if reactions_bitarray[i]:
                prod = dc.get_species_bt_from_reaction(self.RnStr[i], self.RnStr, self.RnMprod)
                #print("bit "+str(i)+" adds "+str(prod))
                specs = specs | prod
                #print("specs updated "+str(specs))
        return dc.get_species_from_bt(specs, self.SpStr)

    def get_species_from_reactions(self, RnStr):
        #Input: List of reactions
        #Output:  List of reactants and products of those reactions
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
        #Input: List of species
        #Output: List of products of reactions triggered by the input
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
    def get_req_from_species(self, SpStr):
        #Input: List of species
        #Output: List of species that are consumed but not produced by the input
        reactions_list = self.get_reactions_from_species(SpStr)
        prod_of_reactions = self.get_prod_from_reactions(reactions_list)
        supp_of_reactions = self.get_supp_from_reactions(reactions_list)
        return [sp for sp in supp_of_reactions if sp not in prod_of_reactions]

    def get_reactions_consuming_species(self, SpStr):
        #Input: List of species
        #Output:  List of reactions triggered by the input
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
        #Input: List of species
        #Output:  List of reactions producing species in the input
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
        #Input: List of species
        #Output:  List of species connected to the input
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
        #Input: List of species
        #Output:  List of species directly connected to the input
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
        #Input: List of species
        #Output:  List of species that are transitively obtained as products of the input 
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
         #Input: List of species
        #Output:  List of species that are directly obtained as products of the input  
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
         #Input: List of species
        #Output:  List of species from which input species can be obtained as transitive product
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
         #Input: List of species
        #Output:  List of species from which input species can be obtained as product
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
         #Input: List of species
        #Output:  List of species from which input species can be obtained as product, minus the input
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
        #Input: 
        #Output:  List of species that can be produced without any reactants
        return(self.get_prod_from_species([]))

    def get_outflow(self):
        #Input: 
        #Output:  List of species that can be consumed without generating products
        result = bt(len(self.SpBt))
        result.setall(0)
        reacs = self.RnStr
        for i in range(len(reacs)):
            supp = dc.get_species_bt_from_reaction(reacs[i], self.RnStr, self.RnMsupp)
            prod = dc.get_species_bt_from_reaction(reacs[i], self.RnStr, self.RnMprod)
            if not prod.any():
                result = result | supp
        return dc.get_species_from_bt(result, self.SpStr)

    ########################################################################################################
    ############# verify relational properties of species/reactions#########################################
    ########################################################################################################

    def is_closed(self, SpStr):
        #Input: List of species
        #Output:  Binary that verifies closure property
        species_bitarray = dc.get_bt_from_names(SpStr, self.SpStr)
        reactions_list = self.get_reactions_from_species(SpStr)
        prod_of_reactions = self.get_prod_from_reactions(reactions_list)
        prod_bitarray = dc.get_bt_from_names(prod_of_reactions, self.SpStr)
        return (prod_bitarray | species_bitarray) == species_bitarray

    def is_semi_self_maintaining(self, SpStr):
        #Input: List of species
        #Output:  Binary that verifies semi_self_maintaining property
        reactions_list = self.get_reactions_from_species(SpStr)
        prod_of_reactions = self.get_prod_from_reactions(reactions_list)
        prod_bitarray = dc.get_bt_from_names(prod_of_reactions, self.SpStr)
        supp_of_reactions = self.get_supp_from_reactions(reactions_list)
        supp_bitarray = dc.get_bt_from_names(supp_of_reactions, self.SpStr)
        return (supp_bitarray & prod_bitarray) == supp_bitarray
     
    def is_connected(self,SpStr):
        #Input: List of species
        #Output:  Binary that verifies all species are connected
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
