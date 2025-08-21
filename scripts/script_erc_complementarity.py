from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import get_complementarity, is_complementary_type1, is_complementary_type2, is_complementary_type3
from pyCOT.io.functions import read_txt
import time

# Load network and compute ERCs
file_path = 'networks/testing/ERCs_test2.txt'
file_path = 'networks/testing/Farm.txt'
#file_path = 'networks/Navarino/RN_IN_05.txt'
#file_path = 'networks/biomodels_interesting/bigg_iAF692.txt'  # Adjust path as needed

print("Loading network and computing ERCs...")
RN = read_txt(file_path)
hierarchy = ERC_Hierarchy(RN)
hierarchy.build_hierarchy_graph()

print(f"\nNetwork statistics:")
print(f"- Reactions: {len(RN.reactions())}")
print(f"- Species: {len(RN.species())}")
print(f"- ERCs found: {len(hierarchy.ercs)}")

# Test complementarity detection between ERC pairs
print("\n" + "="*60)
print("TESTING COMPLEMENTARITY DETECTION")
print("="*60)

complementarity_results = {
    'type1': [],
    'type2': [],
    'type3': []
}

print("\nTesting complementarities between all ERC pairs...")
start_time = time.time()

for i, erc1 in enumerate(hierarchy.ercs):
    for j, erc2 in enumerate(hierarchy.ercs[i+1:], i+1):
        
        # Get all complementarities
        complementarities = get_complementarity(erc1, erc2, hierarchy, RN)
        
        if complementarities:
            print(f"\n🔗 Complementarities between {erc1.label} and {erc2.label}:")
            
            for comp in complementarities:
                complementarity_results[f'type{comp.comp_type}'].append(comp)
                
                if comp.comp_type == 1:
                    print(f"  ✨ Type 1 (Requirement Reduction): {comp}")
                    
                elif comp.comp_type == 2:
                    print(f"  🔄 Type 2 (Requirement Change): {comp}")
                    
                elif comp.comp_type == 3:
                    print(f"  📈 Type 3 (Product Expansion): {comp}")
                 
end_time = time.time()

# Summary of results
print("\n" + "="*60)
print("COMPLEMENTARITY DETECTION SUMMARY")
print("="*60)

print(f"Computation time: {end_time - start_time:.2f} seconds")
print(f"\nComplementarity counts:")
print(f"- Type 1 (Requirement Reduction): {len(complementarity_results['type1'])}")
print(f"- Type 2 (Requirement Change): {len(complementarity_results['type2'])}")
print(f"- Type 3 (Product Expansion): {len(complementarity_results['type3'])}")

total_complementarities = sum(len(comps) for comps in complementarity_results.values())

# Show detailed examples if any complementarities were found
if total_complementarities > 0:
    print(f"\n📋 Detailed examples:")
    
    if complementarity_results['type1']:
        print(f"\nType 1 Complementarity Example:")
        comp = complementarity_results['type1'][0]
        print(f"  Complementarity: {comp}")
       
    
    if complementarity_results['type2']:
        print(f"\nType 2 Complementarity Example:")
        comp = complementarity_results['type2'][0]
        print(f"  Complementarity: {comp}")
       
    if complementarity_results['type3']:
        print(f"\nType 3 Complementarity Example:")
        comp = complementarity_results['type3'][0]
        print(f"  Complementarity: {comp}")
        print(f"  Same requirements but different products")
       

else:
    print(f"\n❌ No complementarities found in this network.")
    print(f"   This could mean:")
    print(f"   - ERCs are relatively independent")
    print(f"   - Network structure doesn't support complementary interactions")
    print(f"   - ERCs might have containment relationships that prevent complementarity")

# Additional analysis: Show ERC details for better understanding
print(f"\n" + "="*60)
print("ERC DETAILS FOR REFERENCE")
print("="*60)

print(f"\nShowing details for first few ERCs:")
for i, erc in enumerate(hierarchy.ercs[:5]):  # Show first 5 ERCs
    print(f"\nERC {erc.label}:")
    closure_names = species_list_to_names(erc.get_closure(RN))
    required_species = erc.get_required_species(RN)
    produced_species = erc.get_produced_species(RN)
    
    print(f"  Closure: {closure_names}")
    print(f"  Required species: {required_species}")
    print(f"  Produced species: {produced_species}")

print(f"\n✅ Complementarity analysis complete!")