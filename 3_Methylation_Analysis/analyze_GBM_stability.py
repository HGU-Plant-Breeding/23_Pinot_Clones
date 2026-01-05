#!/usr/bin/env python3
"""
Script Name: analyze_gbm_stability.py
Description: Analyzes the stability of Gene Body Methylation (gbM) states across clones.
             It processes the 'Class Matrix' (from GBM_calc.py) to identify:
             1. gbM Shifts: Genes that are BM in at least one clone AND UM in at least one other.
             2. BM Fraction: The proportion of classifiable states that are Methylated.

Methodology:
    - Inputs: gbm_class_matrix.tsv (Rows=Genes, Cols=Clones, Values=BM/IM/UM/LowCov)
    - Shift Definition: (Has_BM == True) AND (Has_UM == True).
    - Intermediate (IM) and LowCov states are ignored for the shift definition.

Author: Paolo Callipo
Date: 2025
Dependencies: python3, pandas, numpy
"""

import sys
import argparse
import pandas as pd
import numpy as np

def analyze_stability(input_matrix, output_file):
    print(f"Loading class matrix: {input_matrix}...")
    
    try:
        # Load matrix (Index=GeneID)
        cls = pd.read_csv(input_matrix, sep="\t", index_col=0)
    except FileNotFoundError:
        sys.exit(f"Error: File '{input_matrix}' not found.")

    total_genes = len(cls)
    print(f"Total genes loaded: {total_genes:,}")

    # --- 1. Identify States per Gene ---
    # Boolean masks for presence of specific states across the row (clones)
    has_BM = (cls == "BM").any(axis=1)
    has_UM = (cls == "UM").any(axis=1)
    
    # Definition of Shift: Gene exists as BOTH BM and UM in the population
    gbm_shift = has_BM & has_UM

    # --- 2. Calculate BM Fraction ---
    # Helper to calculate frequency of BM among valid calls (BM+UM)
    def calculate_bm_fraction(row):
        # Filter for only definitive states
        valid_states = row[row.isin(["BM", "UM"])]
        if len(valid_states) == 0:
            return np.nan
        # Mean of boolean (True=1=BM) gives fraction
        return (valid_states == "BM").mean()

    print("Calculating stability metrics...")
    bm_fraction = cls.apply(calculate_bm_fraction, axis=1)
    
    # Count how many clones provided a valid classification (BM or UM)
    n_classifiable = cls.apply(lambda r: r.isin(["BM", "UM"]).sum(), axis=1)

    # --- 3. Construct Summary Table ---
    summary = pd.DataFrame({
        "BM_fraction": bm_fraction,
        "is_gbM_shift": gbm_shift,
        "n_valid_clones": n_classifiable,
        "has_BM": has_BM,
        "has_UM": has_UM
    })

    # --- 4. Save Output ---
    print(f"Saving results to: {output_file}")
    summary.to_csv(output_file, sep="\t", float_format='%.4f')

    # --- 5. Print Report ---
    num_shifts = gbm_shift.sum()
    percent_shift = (num_shifts / total_genes) * 100
    
    print("\n" + "="*40)
    print(" gbM STABILITY REPORT")
    print("="*40)
    print(f"Total Genes:      {total_genes:,}")
    print(f"gbM Shift Candidates: {num_shifts:,} ({percent_shift:.2f}%)")
    print(f"  (Genes BM in ≥1 clone AND UM in ≥1 clone)")
    print("-" * 40)
    print("Breakdown of Shifts:")
    print(summary["is_gbM_shift"].value_counts())
    print("="*40)

def main():
    parser = argparse.ArgumentParser(
        description="Analyze gbM stability and identify epigenetic shifts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_matrix", help="Path to gbm_class_matrix.tsv (from GBM_calc.py).")
    parser.add_argument("--output", default="gene_gbm_stability.tsv", help="Output filename.")

    args = parser.parse_args()

    analyze_stability(args.input_matrix, args.output)

if __name__ == "__main__":
    main()
