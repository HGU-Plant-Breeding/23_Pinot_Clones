#!/usr/bin/env python3
"""
Script Name: run_permutation_test.py
Description: Performs a Permutation (Monte Carlo) test to evaluate the statistical 
             significance of the overlap between a 'Query' set and a 'Database' set.
             
             Crucial Concept: "Universe"
             The randomization is constrained to a 'Universe' BED file. Random regions 
             are sampled *from* this file, ensuring the background model is biologically 
             relevant (e.g., "only selecting from heterozygous promoters", not the whole genome).

Author: [Your Name/Lab Name]
Date: 2024
Dependencies: python3, pandas, pybedtools (requires bedtools in PATH)
"""

import sys
import argparse
import random
import pandas as pd

# Check for pybedtools dependency
try:
    import pybedtools
except ImportError:
    print("Error: The 'pybedtools' library is required.", file=sys.stderr)
    print("Please install it: pip install pybedtools", file=sys.stderr)
    print("Note: 'bedtools' must also be installed in your system PATH.", file=sys.stderr)
    sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Run a permutation test to assess overlap significance constrained to a genomic universe.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('query', help="Query BED file (e.g., Candidate Gene Promoters)")
    parser.add_argument('database', help="Database BED file (e.g., TEs, Indels)")
    parser.add_argument('universe', help="Universe BED file (e.g., All Promoters). Defines the background.")
    parser.add_argument('iterations', type=int, help="Number of random permutations to run (e.g., 10000).")
    parser.add_argument('--output-counts', default="random_counts.txt", help="Filename to save the distribution of random counts (for plotting).")

    args = parser.parse_args()

    print("--- Starting Permutation Test ---")
    print(f"  Query Set    : {args.query}")
    print(f"  Database Set : {args.database}")
    print(f"  Universe Set : {args.universe}")
    print(f"  Iterations   : {args.iterations:,}")
    print("-----------------------------------")

    # --- Load BED Files ---
    try:
        print("Loading BED files...")
        query_bed = pybedtools.BedTool(args.query)
        database_bed = pybedtools.BedTool(args.database)
        universe_bed = pybedtools.BedTool(args.universe)

        # Convert universe to list for efficient Python random sampling
        universe_list = list(universe_bed)
        print(" -> Files loaded successfully.")
    except Exception as e:
        print(f" ! ERROR: Could not load BED files. Details: {e}", file=sys.stderr)
        sys.exit(1)

    num_query = len(query_bed)
    num_universe = len(universe_list)

    if num_query > num_universe:
        print(f" ! ERROR: Query size ({num_query}) is larger than Universe size ({num_universe}).", file=sys.stderr)
        print("   The Query must be a subset of the Universe.", file=sys.stderr)
        sys.exit(1)

    print(f"  Query size:    {num_query:,}")
    print(f"  Universe size: {num_universe:,}")

    # --- Observed Overlap ---
    print("\nCalculating observed overlap...")
    # u=True returns the original feature once if any overlap is found
    observed_intersect = query_bed.intersect(database_bed, u=True)
    observed_count = len(observed_intersect)
    print(f" -> Observed overlap count: {observed_count}")

    # --- Permutations ---
    print(f"\nRunning {args.iterations:,} permutations...")
    random_counts = []

    # Pre-calculate progress milestones
    milestone = max(1, args.iterations // 10)

    for i in range(args.iterations):
        if (i + 1) % milestone == 0:
            print(f"  -> {(i + 1) / args.iterations * 100:.0f}% complete...")

        # 1. Randomly sample 'num_query' regions FROM THE UNIVERSE
        random_sample_regions = random.sample(universe_list, num_query)
        random_sample_bed = pybedtools.BedTool(random_sample_regions)

        # 2. Count overlaps for this random sample
        random_intersect_count = len(random_sample_bed.intersect(database_bed, u=True))
        random_counts.append(random_intersect_count)

    # --- P-value Calculation ---
    print("\nCalculating statistics...")
    random_series = pd.Series(random_counts)
    
    # P-value definition: Proportion of random trials >= observed.
    # We add 1 to both num and denom (pseudo-count) to avoid p=0.
    num_extreme_or_more = (random_series >= observed_count).sum()
    p_value = (num_extreme_or_more + 1) / (args.iterations + 1)

    # --- Report ---
    print("\n" + "="*40)
    print(" PERMUTATION TEST RESULTS")
    print("="*40)
    print(f"  Observed Overlaps       : {observed_count}")
    print(f"  Avg Random Overlap      : {random_series.mean():.2f}")
    print(f"  Random sets >= Observed : {num_extreme_or_more}")
    print(f"  Empirical P-value       : {p_value:.6g}")
    print("="*40)

    # Save output
    random_series.to_csv(args.output_counts, index=False, header=["count"])
    print(f"\nRandom distribution saved to: '{args.output_counts}'")

if __name__ == "__main__":
    main()
