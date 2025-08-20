# ========================================
# 1. LIBRARY LOADING AND CONFIGURATION
# ========================================
import os
import sys

# Add the project root directory to the PYTHONPATH
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import get_complementarity, is_complementary_type1, is_complementary_type2, is_complementary_type3
from pyCOT.io.functions import read_txt
import time

# Load network and compute ERCs
file_path = 'Txt/Farm.txt'
# file_path = 'networks/testing/ERCs_test2.txt'
# file_path = 'networks/testing/Farm.txt'
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
                    reduction = getattr(comp, "reduction_info", {}).get("reduction", 0)
                    # reduction = comp.reduction_info.get('reduction', 0)
                    print(f"     - Reduction in requirements: {reduction} species")
                    
                elif comp.comp_type == 2:
                    print(f"  🔄 Type 2 (Requirement Change): {comp}")
                    req_change = getattr(comp, "reduction_info", {}).get("requirement_change", set())
                    # req_change = comp.reduction_info.get('requirement_change', set())
                    if req_change:
                        print(f"     - New requirements: {req_change}")
                    
                elif comp.comp_type == 3:
                    print(f"  📈 Type 3 (Product Expansion): {comp}")
                    novel_products = getattr(comp, "reduction_info", {}).get('novel_products_direct', set())
                    # novel_products = comp.reduction_info.get('novel_products_direct', set())
                    if novel_products:
                        print(f"     - Novel products: {novel_products}")

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
# if total_complementarities > 0:
#     print(f"\n📋 Detailed examples:")
    
#     if complementarity_results['type1']:
#         print(f"\nType 1 Complementarity Example:")
#         comp = complementarity_results['type1'][0]
#         info = getattr(comp, "reduction_info", None)
#         # info = comp.reduction_info
#         print(f"  Complementarity: {comp}")
#         print(f"  ERC {comp.erc1.label} requirements: {info.get('req1', set())}")
#         print(f"  ERC {comp.erc2.label} requirements: {info.get('req2', set())}")
#         print(f"  Joint requirements: {info.get('joint_req', set())}")
#         print(f"  Reduction: {info.get('reduction', 0)} species")
#         satisfied_by_1 = info.get('satisfied_by_1', set())
#         satisfied_by_2 = info.get('satisfied_by_2', set())
#         if satisfied_by_1:
#             print(f"  {comp.erc1.label} satisfies {comp.erc2.label}'s needs: {satisfied_by_1}")
#         if satisfied_by_2:
#             print(f"  {comp.erc2.label} satisfies {comp.erc1.label}'s needs: {satisfied_by_2}")
    
#     if complementarity_results['type2']:
#         print(f"\nType 2 Complementarity Example:")
#         comp = complementarity_results['type2'][0]
#         info = comp.reduction_info
#         print(f"  Complementarity: {comp}")
#         print(f"  Same total requirements but different requirement set")
#         print(f"  Original requirements: {info.get('req1', set())}")
#         print(f"  Joint requirements: {info.get('joint_req', set())}")
#         req_change = info.get('requirement_change', set())
#         if req_change:
#             print(f"  Requirement changes: {req_change}")
    
#     if complementarity_results['type3']:
#         print(f"\nType 3 Complementarity Example:")
#         comp = complementarity_results['type3'][0]
#         info = comp.reduction_info
#         print(f"  Complementarity: {comp}")
#         print(f"  Same requirements but different products")
#         print(f"  Requirements: {info.get('req1', set())}")
#         print(f"  Original products: {info.get('prod1', set())}")
#         print(f"  Joint products: {info.get('joint_produced', set())}")
#         novel_products = info.get('novel_products_direct', set())
#         if novel_products:
#             print(f"  Novel products: {novel_products}")
#         if info.get('has_synergy', False):
#             print(f"  🔄 This complementarity involves synergistic effects")

# else:
#     print(f"\n❌ No complementarities found in this network.")
#     print(f"   This could mean:")
#     print(f"   - ERCs are relatively independent")
#     print(f"   - Network structure doesn't support complementary interactions")
#     print(f"   - ERCs might have containment relationships that prevent complementarity")
# Show detailed examples if any complementarities were found
if total_complementarities > 0:
    print(f"\n📋 Detailed examples:")
    
    if complementarity_results['type1']:
        print(f"\nType 1 Complementarity Example:")
        comp = complementarity_results['type1'][0]
        info = getattr(comp, "reduction_info", {})
        print(f"  Complementarity: {comp}")
        print(f"  ERC {comp.erc1.label} requirements: {info.get('req1', set())}")
        print(f"  ERC {comp.erc2.label} requirements: {info.get('req2', set())}")
        print(f"  Joint requirements: {info.get('joint_req', set())}")
        print(f"  Reduction: {info.get('reduction', 0)} species")
        satisfied_by_1 = info.get('satisfied_by_1', set())
        satisfied_by_2 = info.get('satisfied_by_2', set())
        if satisfied_by_1:
            print(f"  {comp.erc1.label} satisfies {comp.erc2.label}'s needs: {satisfied_by_1}")
        if satisfied_by_2:
            print(f"  {comp.erc2.label} satisfies {comp.erc1.label}'s needs: {satisfied_by_2}")
    
    if complementarity_results['type2']:
        print(f"\nType 2 Complementarity Example:")
        comp = complementarity_results['type2'][0]
        info = getattr(comp, "requirement_change", {})
        print(f"  Complementarity: {comp}")
        print(f"  Same total requirements but different requirement set")
        print(f"  Requirement changes: {info if info else 'N/A'}")
    
    if complementarity_results['type3']:
        print(f"\nType 3 Complementarity Example:")
        comp = complementarity_results['type3'][0]
        info = getattr(comp, "product_expansion", {})
        print(f"  Complementarity: {comp}")
        print(f"  Same requirements but different products")
        print(f"  Product expansion: {info if info else 'N/A'}")

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