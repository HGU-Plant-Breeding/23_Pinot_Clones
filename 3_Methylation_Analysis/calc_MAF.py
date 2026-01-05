#!/usr/bin/env python3
"""
Script Name: calc_maf_spectrum.py
Description: Calculates the Minor Allele Frequency (MAF) for every site in a 
             methylation matrix.
             
             Output is used to generate the Site Frequency Spectrum (SFS) plots.
             MAF = min(Frequency_Methylated, Frequency_Unmethylated)

Input Format:
    Tab-separated matrix. Row 1 = Samples. Col 1 = Site IDs.
    Values: 0, 1, or . (missing)

Author: [Your Name/Lab Name]
Date: 2024
Dependencies: python3, pandas, numpy
"""

import sys
import argparse
import pandas as pd
import numpy as np

def calculate_maf(input_file, output_file):
    print(f"Loading matrix: {input_file}...")
    
    # Read matrix, treating '.' as NaN
    try:
        df = pd.read_csv(input_file, sep="\t", index_col=0, na_values=".")
    except Exception as e:
        sys.exit(f"Error reading file: {e}")

    print(f"Calculating MAF for {len(df):,} sites...")

    # --- Vectorized Calculation (Fast) ---
    
    # 1. Count valid data per row
    n_valid = df.count(axis=1)
    
    # 2. Sum methylated entries ('1')
    # skipna=True treats NaNs as 0 in sum
    n_meth = df.sum(axis=1, skipna=True)
    
    # 3. Calculate Frequency
    # Handle division by zero (sites with no valid data) by filling with 0
    freq_meth = (n_meth / n_valid).fillna(0)
    
    # 4. Calculate MAF
    # MAF is the smaller of (Freq, 1-Freq)
    maf = np.minimum(freq_meth, 1.0 - freq_meth)

    # --- Save Output ---
    print(f"Saving results to {output_file}...")
    
    # Create output DataFrame
    out_df = pd.DataFrame({
        "ID": df.index,
        "MAF": maf
    })
    
    # Filter out lines that had 0 valid coverage if desired, 
    # or keep them (MAF will be 0.0). Keeping them for now.
    
    out_df.to_csv(output_file, sep="\t", index=False, float_format='%.4f')
    print("Done.")

def main():
    parser = argparse.ArgumentParser(
        description="Calculate MAF per site from methylation matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_matrix", help="Input matrix (tsv).")
    parser.add_argument("output_file", help="Output file (tsv) with ID and MAF columns.")

    args = parser.parse_args()

    calculate_maf(args.input_matrix, args.output_file)

if __name__ == "__main__":
    main()
