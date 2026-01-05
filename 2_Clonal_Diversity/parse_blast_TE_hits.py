#!/usr/bin/env python3
"""
Script Name: parse_blast_te_hits.py
Description: Parses BLAST output to classify insertions as TEs based on the 
             best hit (by bit-score).
             
             It filters hits based on:
             1. Minimum Percent Identity (default 80%)
             2. Minimum Query Coverage (default 80%)

Input Requirement:
    BLAST output must be tabular (-outfmt 6) with SPECIFIC columns:
    qseqid sseqid pident length qlen slen bitscore

Author: Paolo Callipo  
Date: 2025
Dependencies: python3
"""

import sys
import argparse
from collections import defaultdict

def parse_blast(blast_file, min_identity, min_coverage, output_file):
    print(f"Parsing BLAST results from {blast_file}...")
    print("Strategy: Keep best hit per query based on Bit-Score.")
    
    best_hits = {}

    try:
        with open(blast_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'): continue

                parts = line.split('\t')
                
                # Check for correct column count (7 columns required)
                if len(parts) < 7:
                    if line_num == 1:
                        print(f"Warning: Line 1 has {len(parts)} columns. Expected 7.", file=sys.stderr)
                        print("Format must be: qseqid sseqid pident length qlen slen bitscore", file=sys.stderr)
                    continue

                try:
                    qseqid = parts[0]
                    sseqid = parts[1]
                    pident = float(parts[2])
                    length = int(parts[3]) # Alignment length
                    qlen = int(parts[4])   # Query length
                    # slen = int(parts[5]) # Subject length (unused here)
                    bitscore = float(parts[6])

                    # Logic: Update if this is the first time seeing this query
                    # OR if this hit has a higher bitscore than the previous best
                    if qseqid not in best_hits or bitscore > best_hits[qseqid]['bitscore']:
                        best_hits[qseqid] = {
                            'subject_id': sseqid,
                            'percent_identity': pident,
                            'align_length': length,
                            'query_length': qlen,
                            'bitscore': bitscore
                        }
                except ValueError:
                    continue
    except FileNotFoundError:
        sys.exit(f"Error: BLAST file '{blast_file}' not found.")

    # --- Filter and Write ---
    total_insertions = len(best_hits)
    te_families = defaultdict(int)
    count_classified = 0

    print(f"Total queries with hits: {total_insertions}")
    print(f"Filtering (Identity >= {min_identity}%, Coverage >= {min_coverage*100}%)...")

    with open(output_file, "w") as f_out:
        f_out.write("Insertion_ID\tTE_Family\n")

        for query_id, hit in best_hits.items():
            query_length = hit['query_length']
            
            # Avoid division by zero
            if query_length == 0: continue
            
            query_coverage = hit['align_length'] / query_length

            if hit['percent_identity'] >= min_identity and query_coverage >= min_coverage:
                # It passes filters -> Classify as TE
                f_out.write(f"{query_id}\t{hit['subject_id']}\n")
                te_families[hit['subject_id']] += 1
                count_classified += 1

    # --- Summary ---
    percent_te = (count_classified / total_insertions * 100) if total_insertions > 0 else 0
    
    print("\n--- Summary ---")
    print(f"Classified as TEs: {count_classified}")
    print(f"Percentage:        {percent_te:.2f}%")
    print(f"Output saved to:   {output_file}")

    print("\n--- Top 5 TE Families ---")
    for family, count in sorted(te_families.items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"{family}: {count}")

def main():
    parser = argparse.ArgumentParser(
        description="Filter BLAST hits to identify TE insertions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("blast_file", help="BLAST output file (cols: qseqid sseqid pident length qlen slen bitscore)")
    parser.add_argument("--output", default="te_insertion_families.tsv", help="Output filename")
    parser.add_argument("--min-identity", type=float, default=80.0, help="Minimum percent identity")
    parser.add_argument("--min-coverage", type=float, default=0.80, help="Minimum query coverage (0.0 - 1.0)")

    args = parser.parse_args()
    
    parse_blast(args.blast_file, args.min_identity, args.min_coverage, args.output)

if __name__ == '__main__':
    main()
