#!/usr/bin/env python3
# Test script for get_reactions_partially_intersecting_support method


###COMMENT, NOT WORKING PROPERLY YET### 
###note that identifies the following synergy:
# Testing ERC with label E3
# Closure: [Species(index=0, name='A', quantity=10.0)]
# Potential synergies with:
#   - ERC E0 with closure [Species(index=0, name='A', quantity=10.0), Species(index=1, name='B', quantity=10.0)]
#Hence Partially activated is not getting the containment case properly..
from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import Hierarchy_ERC

def main():
    """Main function to test get_reactions_partially_intersecting_support method."""
    # Create a new reaction network
    rn = ReactionNetwork()
    
    # Add species
    rn.add_species("A")
    rn.add_species("B")
    rn.add_species("C")
    rn.add_species("D")
    rn.add_species("E")
    rn.add_species("F")
    
    # Set quantities for testing is_active_reaction (optional)
    for species_name in ["A", "B", "C", "D", "E", "F"]:
        species = rn.get_species(species_name)
        species.quantity = 10.0
    
    # Add reactions
    # Reaction 1: A + B -> C
    rn.add_reaction("R1", ["A", "B"], ["C"], [1.0, 1.0], [1.0], 1.0)
    
    # Reaction 2: B + D -> E
    rn.add_reaction("R2", ["B", "D"], ["E"], [1.0, 1.0], [1.0], 1.0)
    
    # Reaction 3: C + E -> F
    rn.add_reaction("R3", ["C", "E"], ["F"], [1.0, 1.0], [1.0], 1.0)
    
    # Reaction 4: A -> D (single reactant)
    rn.add_reaction("R4", ["A"], ["D"], [1.0], [1.0], 1.0)
    
    # Reaction 5: -> A (inflow reaction - no reactants)
    rn.add_reaction("R5", None, ["A"], [], [1.0], 1.0)
    
    # Print all reactions for reference
    print("All reactions in the network:")
    for reaction in rn.reactions():
        print(f"  {reaction.name()}: {[rn.get_species_by_index(i).name for i in reaction.support_indices()]} -> {[rn.get_species_by_index(i).name for i in reaction.products_indices()]}")
    print()
    
    # Test cases
    test_cases = [
        ["A"],                # Should partially intersect with R1 (needs B too)
        ["B"],                # Should partially intersect with R1 and R2
        ["A", "C"],           # Should partially intersect with R3 (needs E too)
        ["B", "C", "D"],      # Should partially intersect with R3 (needs E)
        ["A", "B", "C", "D"]  # No partially intersecting reactions (R1 and R2 fully satisfied)
    ]
    
    for i, species_set in enumerate(test_cases):
        print(f"Test Case {i+1}: Species = {species_set}")
        
        # Get partially intersecting reactions
        partial_reactions = rn.get_reactions_partially_intersecting_support(species_set)
        
        if partial_reactions:
            print(f"  Reactions partially activated by {species_set}:")
            for reaction in partial_reactions:
                # Print the reaction, its reactants, and what's missing
                reactants = [rn.get_species_by_index(i).name for i in reaction.support_indices()]
                missing = [r for r in reactants if r not in species_set]
                print(f"  - {reaction.name()}: {reactants} -> {[rn.get_species_by_index(i).name for i in reaction.products_indices()]}")
                print(f"    Missing reactants: {missing}")
        else:
            print(f"  No reactions are partially activated by {species_set}")
        
        # For comparison, get fully activated reactions
        full_reactions = rn.get_reactions_from_species(species_set)
        print(f"  Reactions fully activated by {species_set}:")
        for reaction in full_reactions:
            print(f"  - {reaction.name()}: {[rn.get_species_by_index(i).name for i in reaction.support_indices()]} -> {[rn.get_species_by_index(i).name for i in reaction.products_indices()]}")
        print()

def test_potential_synergies():
    """Test the get_potential_synergies method."""
    # Create the same reaction network as in main()
    rn = ReactionNetwork()
    
    # Add species
    for species in ["A", "B", "C", "D", "E", "F"]:
        rn.add_species(species)
        species_obj = rn.get_species(species)
        species_obj.quantity = 10.0
    
    # Add reactions
    rn.add_reaction("R1", ["A", "B"], ["C"], [1.0, 1.0], [1.0], 1.0)
    rn.add_reaction("R2", ["B", "D"], ["E"], [1.0, 1.0], [1.0], 1.0)
    rn.add_reaction("R3", ["C", "E"], ["F"], [1.0, 1.0], [1.0], 1.0)
    rn.add_reaction("R4", ["A"], ["D"], [1.0], [1.0], 1.0)
    rn.add_reaction("R5", None, ["A"], [], [1.0], 1.0)
    
    # Create hierarchy of ERCs
    hierarchy = Hierarchy_ERC(rn)
    for erc in hierarchy.ercs:
        print(erc.label, str([r.name for r in erc.all_generators]))
    # Test potential synergies for each ERC
    print("\nTesting potential synergies:")
    for erc in hierarchy.ercs:
        print(f"\nTesting ERC with label {erc.label}")
        print(f"Closure: {erc.get_closure(rn)}")
        synergies = hierarchy.get_potential_synergies(rn,erc)
        if synergies:
            print("Potential synergies with:")
            for syn_erc in synergies:
                print(f"  - ERC {syn_erc.label} with closure {syn_erc.get_closure(rn)}")
        else:
            print("No potential synergies found")

if __name__ == "__main__":
    main()
    test_potential_synergies()