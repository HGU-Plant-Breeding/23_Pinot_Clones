#!/usr/bin/env python3
"""
Script Name: merge_methylation_matrix.py
Description: Combines multiple binarized methylation files (ID, Score) into a 
             single population matrix.
             
Logic:
    - Reads multiple input files.
    - Performs an OUTER JOIN on the Site ID (Chrom_Pos).
    - Fills missing data (sites present in some clones but not others) with '.'.
    - Renames columns to match sample names.

Output Format:
    Tab-separated matrix. Row 1 = Sample IDs. Col 1 = Site IDs.

Author: [Your Name/Lab Name]
Date: 2024
Dependencies: python3, pandas
"""

import sys
import argparse
import os
import pandas as pd

def get_sample_name(filepath):
    """
    Extracts a clean sample name from the filename.
    Removes common extensions like .bed, .tsv, .txt, .gz
    """
    base = os.path.basename(filepath)
    # Remove extensions repeatedly until clean (e.g., sample.CG.binary.bed -> sample.CG.binary -> sample.CG -> sample)
    # Adjust logic here if you want to keep specific tags like "CG"
    name = os.path.splitext(base)[0]
    if name.endswith('.bed') or name.endswith('.tsv') or name.endswith('.txt'):
        name = os.path.splitext(name)[0]
    return name

def combine_samples(output_file, input_files):
    print(f"Merging {len(input_files)} samples into population matrix...")
    
    df_list = []

    for file_path in input_files:
        sample_name = get_sample_name(file_path)
        # print(f"  - Loading {sample_name}...")

        try:
            # Assumes input has a header "ID\tScore" (from binarize_methylation.py)
            # using sep='\t' and default header='infer' (row 0)
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            
            # Rename the single data column to the sample name
            # (Assumes there is only one data column per file)
            if df.shape[1] == 1:
                df.columns = [sample_name]
            else:
                # Fallback if multiple cols, take the last one or explicit 'Score'
                print(f"Warning: {file_path} has multiple columns. Using the last one as score.")
                df = df.iloc[:, -1:]
                df.columns = [sample_name]
                
            df_list.append(df)
            
        except Exception as e:
            print(f"Error reading {file_path}: {e}", file=sys.stderr)
            continue

    if not df_list:
        sys.exit("Error: No valid data loaded.")

    print("Concatenating data (Outer Join)...")
    # Outer join ensures we keep all sites found in ANY clone
    combined_df = pd.concat(df_list, axis=1, join='outer')

    print("Filling missing values...")
    combined_df.fillna('.', inplace=True)

    print(f"Saving matrix to {output_file}...")
    combined_df.to_csv(output_file, sep='\t')
    
    print("Done.")

def main():
    parser = argparse.ArgumentParser(
        description="Merge individual binarized methylation files into a population matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("output_file", help="Output filename (e.g., matrix.tsv).")
    parser.add_argument("input_files", nargs='+', help="List of input files (space separated).")

    args = parser.parse_args()

    combine_samples(args.output_file, args.input_files)

if __name__ == '__main__':
    main()
