#!/usr/bin/env python3
"""
Simple Usage Example for Synergy Algorithm Comparison

This script demonstrates how to use the comparison tool with different
reaction networks and analyze the results.
"""

import sys
import os
from pathlib import Path

# Add the parent directory to the path to import our modules
# Adjust this path based on your project structure
sys.path.append(str(Path(__file__).parent))

from pyCOT.synergy_comparison import SynergyAlgorithmComparator


def example_1_simple_network():
    """Test with a simple reaction network."""
    
    print("🧪 EXAMPLE 1: Simple Reaction Network")
    print("="*60)
    
    # Define a simple reaction network
    simple_network = """
    R1: a+b=>c;
    R2: c+d=>e;
    R3: e=>f+g;
    R4: a+d=>h;
    R5: h+b=>i;
    """
    
    # Create comparator and run comparison
    comparator = SynergyAlgorithmComparator(verbose=True)
    result = comparator.compare_algorithms(simple_network, "Simple Network")
    
    return result


def example_2_file_network():
    """Test with a reaction network from file."""
    
    print("🧪 EXAMPLE 2: File-based Reaction Network")
    print("="*60)
    
    # Check if the test file exists
    test_file = "ERCs_test2.txt"
    if not os.path.exists(test_file):
        print(f"❌ Test file {test_file} not found. Skipping this example.")
        return None
    
    # Create comparator and run comparison
    comparator = SynergyAlgorithmComparator(verbose=True)
    result = comparator.compare_algorithms(test_file, "File Network")
    
    return result


def example_3_multiple_networks():
    """Test with multiple reaction networks."""
    
    print("🧪 EXAMPLE 3: Multiple Network Comparison")
    print("="*60)
    
    # Define multiple test networks
    test_networks = {
        "Linear Chain": """
        R1: a=>b;
        R2: b=>c;
        R3: c=>d;
        R4: a+x=>y;
        """,
        
        "Branched Network": """
        R1: a+b=>c+d;
        R2: c+e=>f;
        R3: d+e=>g;
        R4: f+g=>h;
        R5: a+x=>y+z;
        """,
        
        "Complex Network": """
        R1: a1+b1=>a2+b2;
        R2: a2+b1=>a3+b3;
        R3: a3+b3=>x3+a2;
        R4: a1+b4=>a5+b5;
        R5: a5+b5=>x5+b2;
        R6: a6+b1=>a7+b7;
        R7: a7+b7=>x7+a6;
        """
    }
    
    # Run test suite
    comparator = SynergyAlgorithmComparator(verbose=True)
    results = comparator.run_test_suite(test_networks)
    
    return results


def example_4_detailed_investigation():
    """Example showing detailed investigation of differences."""
    
    print("🧪 EXAMPLE 4: Detailed Investigation")
    print("="*60)
    
    # Use a network that might show differences
    test_network = """
    R1: a1+b1=>a2+b2;
    R2: a2+b1=>a3+b3;
    R3: a3+b3=>x3+a2;
    R4: a1+b4=>a5+b5;
    R5: a5+b5=>x5+b2;
    R6: a6+b1=>a7+b7;
    """
    
    # Run comparison
    comparator = SynergyAlgorithmComparator(verbose=True)
    result = comparator.compare_algorithms(test_network, "Investigation Test")
    
    if result and not result.perfect_match:
        print("\n🔬 DIFFERENCES DETECTED - Running detailed investigation...")
        
        # Load network for investigation
        rn = comparator.load_reaction_network(test_network)
        from pyCOT.ERC_Hierarchy import ERC
        ercs = ERC.ERCs(rn)
        
        # Run detailed investigation
        comparator.print_detailed_investigation(result, rn, ercs)
    
    return result


def benchmark_performance():
    """Benchmark performance on networks of different sizes."""
    
    print("🧪 PERFORMANCE BENCHMARK")
    print("="*60)
    
    # Generate networks of increasing complexity
    def generate_network(size):
        reactions = []
        for i in range(size):
            reactions.append(f"R{i+1}: a{i}+b{i}=>a{i+1}+b{i+1};")
        return "\n".join(reactions)
    
    benchmark_networks = {
        f"Size_{size}": generate_network(size)
        for size in [5, 10, 15, 20]
    }
    
    # Run benchmark
    comparator = SynergyAlgorithmComparator(verbose=False)  # Less verbose for benchmarking
    results = comparator.run_test_suite(benchmark_networks)
    
    # Analyze performance trends
    print("\n📈 PERFORMANCE ANALYSIS:")
    sizes = [5, 10, 15, 20]
    valid_results = [r for r in results if not r.get('error', False)]
    
    for i, result in enumerate(valid_results):
        size = sizes[i] if i < len(sizes) else "Unknown"
        print(f"Size {size:2}: BF={result['bf_time']:.4f}s, HA={result['ha_time']:.4f}s, "
              f"Speedup={result['speedup']:.2f}x, Synergies={result['bf_count']}")
    
    return results


def main():
    """Run all examples."""
    
    print("🚀 SYNERGY ALGORITHM COMPARISON EXAMPLES")
    print("="*80)
    
    try:
        # Run examples
        print("\n" + "="*80)
        example_1_simple_network()
        
        print("\n" + "="*80)
        example_2_file_network()
        
        print("\n" + "="*80)
        example_3_multiple_networks()
        
        print("\n" + "="*80)
        example_4_detailed_investigation()
        
        print("\n" + "="*80)
        benchmark_performance()
        
        print("\n✅ All examples completed!")
        
    except KeyboardInterrupt:
        print("\n⚠️  Examples interrupted by user")
    except Exception as e:
        print(f"\n❌ Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()