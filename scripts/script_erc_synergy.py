from pyCOT.rn_rustworkx import ReactionNetwork
from pyCOT.ERC_Hierarchy import ERC, ERC_Hierarchy, species_list_to_names
from pyCOT.ERC_Synergy_Complementarity import get_basic_synergies, get_maximal_synergies, get_fundamental_synergies
from pyCOT.io.functions import read_txt
import time

# Load network and compute ERCs
file_path = 'networks/testing/ERCs_test2.txt'
file_path = 'networks/Navarino/RN_IN_05.txt'
#file_path = 'networks/testing/Farm.txt'
#file_path = 'networks/biomodels_interesting/bigg_iAF692.txt'  # Adjust path as needed

print("Loading network and computing ERCs...")
RN = read_txt(file_path)
hierarchy = ERC_Hierarchy(RN)
hierarchy.build_hierarchy_graph()

print(f"\nNetwork statistics:")
print(f"- Reactions: {len(RN.reactions())}")
print(f"- Species: {len(RN.species())}")
print(f"- ERCs found: {len(hierarchy.ercs)}")

# Test synergy detection between ERC pairs
print("\n" + "="*60)
print("TESTING SYNERGY DETECTION")
print("="*60)

synergy_results = {
    'basic': [],
    'maximal': [],
    'fundamental': []
}

print("\nTesting synergies between all ERC pairs...")
start_time = time.time()

for i, erc1 in enumerate(hierarchy.ercs):
    for j, erc2 in enumerate(hierarchy.ercs[i+1:], i+1):
        
        # Get basic synergies
        basic_syn = get_basic_synergies(erc1, erc2, hierarchy, RN)
        if basic_syn:
            synergy_results['basic'].extend(basic_syn)
            print(f"\n🔄 Basic Synergies between {erc1.label} and {erc2.label}:")
            for syn in basic_syn:
                print(f"  {syn}")
        
        # Get maximal synergies
        maximal_syn = get_maximal_synergies(erc1, erc2, hierarchy, RN)
        if maximal_syn:
            synergy_results['maximal'].extend(maximal_syn)
            print(f"\n⭐ Maximal Synergies between {erc1.label} and {erc2.label}:")
            for syn in maximal_syn:
                print(f"  {syn}")
        
        # Get fundamental synergies
        fundamental_syn = get_fundamental_synergies(erc1, erc2, hierarchy, RN)
        if fundamental_syn:
            synergy_results['fundamental'].extend(fundamental_syn)
            print(f"\n💎 Fundamental Synergies between {erc1.label} and {erc2.label}:")
            for syn in fundamental_syn:
                print(f"  {syn}")

end_time = time.time()

# Summary of results
print("\n" + "="*60)
print("SYNERGY DETECTION SUMMARY")
print("="*60)

print(f"Computation time: {end_time - start_time:.2f} seconds")
print(f"\nSynergy counts:")
print(f"- Basic synergies: {len(synergy_results['basic'])}")
print(f"- Maximal synergies: {len(synergy_results['maximal'])}")  
print(f"- Fundamental synergies: {len(synergy_results['fundamental'])}")

# Show some detailed examples if any synergies were found
if any(synergy_results.values()):
    print(f"\n📋 Detailed examples:")
    
    if synergy_results['basic']:
        print(f"\nBasic Synergy Example:")
        syn = synergy_results['basic'][0]
        print(f"  Synergy: {syn}")
        print(f"  Reactant ERCs:")
        for reactant in syn.reactants:
            closure_names = species_list_to_names(reactant.get_closure(RN))
            print(f"    {reactant.label}: {closure_names}")
        product_closure = species_list_to_names(syn.product.get_closure(RN))
        print(f"  Product ERC:")
        print(f"    {syn.product.label}: {product_closure}")
    
    if synergy_results['maximal']:
        print(f"\nMaximal Synergy Example:")
        syn = synergy_results['maximal'][0]
        print(f"  Synergy: {syn}")
        print(f"  This is a maximal synergy - no larger product ERC can be produced")
        print(f"  by the same reactant combination.")
    
    if synergy_results['fundamental']:
        print(f"\nFundamental Synergy Example:")
        syn = synergy_results['fundamental'][0]
        print(f"  Synergy: {syn}")
        print(f"  This is a fundamental synergy - no more basic ERC combination")
        print(f"  can produce the same result.")

else:
    print(f"\n❌ No synergies found in this network.")
    print(f"   This could mean:")
    print(f"   - ERCs are relatively independent")
    print(f"   - Network structure doesn't support synergistic interactions")
    print(f"   - ERCs might have containment relationships that prevent synergies")

print(f"\n✅ Synergy analysis complete!")