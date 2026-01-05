#!/usr/bin/env python3
"""
Script Name: filter_methylation_matrix.py
Description: Filters a binary methylation matrix to retain only high-confidence, 
             variable sites (Polymorphic Epialleles) for downstream analysis.
             
             Filters Applied:
             1. Missing Data: Removes sites with > X% missing values ('.').
             2. Minor Allele Frequency (MAF): Removes invariant or rare sites 
                (MAF < threshold).

             *Performance Optimized*: Uses vectorized pandas operations instead of loops.

Author: Paolo Callipo
Date: 2025
Dependencies: python3, pandas, numpy
"""

import sys
import argparse
import pandas as pd
import numpy as np

def filter_matrix(input_file, output_file, missing_thresh, maf_thresh):
    print(f"Loading matrix: {input_file}...")
    
    # Read matrix (Rows=Sites, Cols=Samples)
    # Treating '.' as NaN immediately allows for fast numeric operations
    df = pd.read_csv(input_file, sep="\t", index_col=0, na_values=".")
    
    total_sites = len(df)
    print(f"Total sites loaded: {total_sites:,}")

    # --- 1. Filter by Missingness ---
    print(f"Filtering missing data > {missing_thresh*100}%...")
    
    # Calculate fraction of NaNs per row
    missing_fraction = df.isna().mean(axis=1)
    
    # Keep rows where missingness is <= threshold
    df = df.loc[missing_fraction <= missing_thresh]
    
    print(f"Sites remaining after missingness filter: {len(df):,} ({len(df)/total_sites:.1%})")

    if len(df) == 0:
        print("Error: No sites remaining. Check your data or thresholds.")
        sys.exit(1)

    # --- 2. Filter by MAF (Minor Allele Frequency) ---
    print(f"Filtering MAF < {maf_thresh} (Invariant/Rare sites)...")
    
    # Calculate stats on the remaining data
    # count() gives number of non-NA values per row
    # sum() gives number of 1s (Methylated) per row
    valid_counts = df.count(axis=1)
    methylated_counts = df.sum(axis=1, skipna=True)
    
    # Frequency of '1' state
    freq_1 = methylated_counts / valid_counts
    
    # MAF is the smaller of (Freq_1, 1-Freq_1)
    # We use np.minimum to compare the series element-wise
    maf = np.minimum(freq_1, 1 - freq_1)
    
    # Keep rows where MAF >= threshold
    df_final = df.loc[maf >= maf_thresh]
    
    final_count = len(df_final)
    print(f"Sites remaining after MAF filter: {final_count:,} ({final_count/total_sites:.1%})")

    # --- 3. Save ---
    print(f"Saving to {output_file}...")
    
    # Fill NaNs back to '.' for consistent output format if needed, 
    # or keep empty/NA depending on downstream tool requirements.
    # Here we revert to '.' to match the original script format.
    df_final.fillna(".", inplace=True)
    
    # Convert floats back to integers/strings (0.0 -> "0")
    # This creates a clean output
    df_final = df_final.applymap(lambda x: str(int(x)) if x != "." else ".")
    
    df_final.to_csv(output_file, sep="\t")
    print("Done.")

def main():
    parser = argparse.ArgumentParser(
        description="Filter methylation matrix by missingness and allele frequency.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_matrix", help="Path to population matrix (tsv).")
    parser.add_argument("output_matrix", help="Path for filtered output.")
    parser.add_argument("--missing", type=float, default=0.1, 
                        help="Maximum allowed missing data fraction (0.0 - 1.0).")
    parser.add_argument("--maf", type=float, default=0.05, 
                        help="Minimum Minor Allele Frequency (0.0 - 0.5).")

    args = parser.parse_args()

    filter_matrix(args.input_matrix, args.output_matrix, args.missing, args.maf)

if __name__ == "__main__":
    main()
