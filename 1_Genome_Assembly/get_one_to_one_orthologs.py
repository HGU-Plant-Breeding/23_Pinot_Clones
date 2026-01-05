#!/usr/bin/env python3
"""
Script Name: get_one_to_one_orthologs.py
Description: Parses OrthoFinder's Orthogroups.tsv to identify single-copy 
             orthologs (1:1 relationships) between two haplotypes.
             It outputs a TSV file with two columns: Gene_HapA, Gene_HapB.
Author: Paolo Callipo
Date: 2025
Dependencies: python3
Input: Orthogroups.tsv from OrthoFinder.
"""

import sys
import argparse
import os

def normalize_id(protein_id):
    """
    Removes the transcript suffix (.1, .2) from protein IDs if present.
    Example: 'GeneID.1' -> 'GeneID'
    """
    last_dot_index = protein_id.rfind('.')
    if last_dot_index != -1 and protein_id[last_dot_index+1:].isdigit():
        return protein_id[:last_dot_index]
    return protein_id

def extract_orthologs(ortho_file, output_file, col_a, col_b):
    print(f"Parsing {os.path.basename(ortho_file)} for 1-to-1 orthologs...")
    
    count = 0
    try:
        with open(ortho_file, 'r') as f, open(output_file, 'w') as out_f:
            # Write Header
            out_f.write("HapA_GeneID\tHapB_GeneID\n")
            
            # Skip file header
            header = f.readline()
            
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')
                
                # Ensure line has enough columns
                if len(parts) <= max(col_a, col_b):
                    continue

                hap_a_str = parts[col_a].strip()
                hap_b_str = parts[col_b].strip()

                # Logic: We want exactly one gene in each column.
                # OrthoFinder separates multiple genes with commas.
                # So we check that both fields are not empty AND contain no commas.
                if hap_a_str and hap_b_str and ',' not in hap_a_str and ',' not in hap_b_str:
                    gene_a = normalize_id(hap_a_str)
                    gene_b = normalize_id(hap_b_str)
                    out_f.write(f"{gene_a}\t{gene_b}\n")
                    count += 1
                    
        print(f"Found {count:,} one-to-one ortholog pairs.")
        print(f"Output saved to {output_file}")

    except FileNotFoundError:
        print(f"Error: File not found at '{ortho_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Extract 1-to-1 orthologs from OrthoFinder output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('orthogroups_file', help="Path to Orthogroups.tsv")
    parser.add_argument('--output', default="one_to_one_orthologs.tsv", help="Output filename.")
    parser.add_argument('--col-a', type=int, default=1, help="0-based column index for Haplotype A (usually 1).")
    parser.add_argument('--col-b', type=int, default=2, help="0-based column index for Haplotype B (usually 2).")
    
    args = parser.parse_args()
    
    extract_orthologs(args.orthogroups_file, args.output, args.col_a, args.col_b)

if __name__ == "__main__":
    main()
