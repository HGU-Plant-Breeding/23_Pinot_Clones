#!/usr/bin/env python3
"""
Script Name: process_methylation_bed.py
Description: Parses and processes raw methylation BED files (e.g., from modkit).
             It handles context-specific logic:
             - CG:  Merges symmetrical sites (+/- strands). Unpaired sites are skipped.
             - CHG: Merges symmetrical sites (with 1bp offset). Unpaired sites are kept.
             - CHH: No merging. Extracts counts directly.

             *Memory Optimized*: Processes file line-by-line instead of loading into RAM.

Author: Paolo Callipo
Date: 2025
Dependencies: python3
"""

import sys
import argparse

def parse_line(line_str):
    """Helper to split line and convert numerical fields."""
    parts = line_str.strip().split('\t')
    if len(parts) < 11: return None
    return parts

def process_file(input_file, output_file, context):
    print(f"Processing {context} sites from {input_file}...")
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        
        # --- CHH STRATEGY: Simple Filter ---
        if context == 'CHH':
            for line in fin:
                parts = parse_line(line)
                if not parts: continue
                
                # Check context in column 4 (index 3)
                if "CHH" in parts[3]:
                    # Out: Chrom, Start, Coverage, Score
                    fout.write(f"{parts[0]}\t{parts[1]}\t{parts[4]}\t{parts[10]}\n")
            return

        # --- CG / CHG STRATEGY: Streaming Merge ---
        # We use an iterator to look ahead without loading the whole file
        iterator = iter(fin)
        try:
            prev_line = next(iterator)
        except StopIteration:
            return

        for curr_line in iterator:
            p_parts = parse_line(prev_line)
            c_parts = parse_line(curr_line)

            if not p_parts or not c_parts:
                prev_line = curr_line
                continue

            merged = False
            
            # Extract common data
            # Col 5 (index 4) = Coverage, Col 11 (index 10) = Score
            # Col 6 (index 5) = Strand (+/-)
            
            # --- MERGING LOGIC ---
            # Check strands: Prev must be +, Curr must be -
            if p_parts[5] == '+' and c_parts[5] == '-':
                
                is_pair = False
                
                # CG Condition: Ends match exactly (0-based BED)
                if context == 'CG' and p_parts[2] == c_parts[1]:
                    is_pair = True
                
                # CHG Condition: Offset by 1bp
                elif context == 'CHG' and (int(p_parts[2]) + 1) == int(c_parts[1]):
                    is_pair = True

                if is_pair:
                    try:
                        cov1, score1 = float(p_parts[4]), float(p_parts[10])
                        cov2, score2 = float(c_parts[4]), float(c_parts[10])
                        
                        total_cov = cov1 + cov2
                        w_score = (score1 * cov1 + score2 * cov2) / total_cov if total_cov > 0 else 0.0
                        
                        # Output using the + strand position
                        fout.write(f"{p_parts[0]}\t{p_parts[1]}\t{total_cov}\t{w_score}\n")
                        merged = True
                    except ValueError:
                        pass # Skip malformed numbers

            # --- HANDLING UNMERGED ---
            if merged:
                # If merged, we consumed both lines. 
                # Advance iterator for 'prev_line'
                try:
                    prev_line = next(iterator)
                except StopIteration:
                    prev_line = None
            else:
                # If NOT merged:
                # CG behavior (Original script): Skip unpaired lines.
                # CHG behavior (Original script): Write unpaired lines.
                if context == 'CHG':
                    fout.write(f"{p_parts[0]}\t{p_parts[1]}\t{p_parts[4]}\t{p_parts[10]}\n")
                
                # Move window: Current becomes Previous
                prev_line = curr_line

        # Handle the very last line if CHG
        if prev_line and context == 'CHG':
            parts = parse_line(prev_line)
            if parts:
                fout.write(f"{parts[0]}\t{parts[1]}\t{parts[4]}\t{parts[10]}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Process and merge methylation bed files (CG/CHG/CHH).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_file", help="Input BED file (e.g., from modkit).")
    parser.add_argument("output_file", help="Output processed BED file.")
    parser.add_argument("--context", choices=['CG', 'CHG', 'CHH'], required=True, 
                        help="Methylation context. CG/CHG triggers merging logic; CHH extracts data.")

    args = parser.parse_args()

    process_file(args.input_file, args.output_file, args.context)

if __name__ == '__main__':
    main()
