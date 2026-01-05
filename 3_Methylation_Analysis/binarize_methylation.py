#!/usr/bin/env python3
"""
Script Name: binarize_methylation.py
Description: Transforms continuous methylation data into discrete binary states 
             for population epigenetic analysis.
             
Logic:
    - ID generation: Chromosome_Position
    - Missing Data (.): Coverage < MIN_COV or Score between thresholds (Ambiguous).
    - Unmethylated (0): Score <= LOW_THRESH (default 30%).
    - Methylated (1): Score >= HIGH_THRESH (default 70%).

Input Format:
    Tab-delimited file with at least 4 columns: CHROM, POS, COVERAGE, SCORE.

Author: Paolo Callipo
Date: 2025
Dependencies: python3
"""

import sys
import argparse
import os

def transform_methylation(input_file, output_file, min_cov, low_thresh, high_thresh):
    print(f"Binarizing methylation data from {os.path.basename(input_file)}...")
    print(f"Thresholds: Unmethylated <= {low_thresh}%, Methylated >= {high_thresh}%, Min Cov >= {min_cov}")

    processed_count = 0
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        # Write header for the output
        fout.write("ID\tScore\n")
        
        for line in fin:
            if not line.strip() or line.startswith('#') or line.startswith('track'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue

            chrom = fields[0]
            pos = fields[1]
            
            try:
                coverage = float(fields[2])
                methyl_score = float(fields[3])
            except ValueError:
                # Skip header lines or malformed numbers silently
                continue

            processed_count += 1
            
            # Determine State
            if coverage < min_cov:
                score = "."
            elif methyl_score <= low_thresh:
                score = "0"
            elif methyl_score >= high_thresh:
                score = "1"
            else:
                score = "." # Intermediate/Ambiguous

            # Construct ID and write
            id_field = f"{chrom}_{pos}"
            fout.write(f"{id_field}\t{score}\n")

    print(f"Finished. Processed {processed_count} sites.")
    print(f"Output saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert methylation levels to binary states (0/1/.).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_file", help="Input BED file (Chrom, Pos, Cov, Score).")
    parser.add_argument("output_file", help="Output ID-Score file.")
    
    # Optional thresholds
    parser.add_argument("--min-cov", type=float, default=4.0, 
                        help="Minimum coverage required to call a state.")
    parser.add_argument("--low", type=float, default=30.0, 
                        help="Threshold for Unmethylated state (<=).")
    parser.add_argument("--high", type=float, default=70.0, 
                        help="Threshold for Methylated state (>=).")

    args = parser.parse_args()

    transform_methylation(args.input_file, args.output_file, args.min_cov, args.low, args.high)

if __name__ == '__main__':
    main()
