#!/usr/bin/env python3
"""
Script Name: calc_dmc_density.py
Description: Calculates the density of DMCs (events per 100kb) per feature type 
             (Exon, Intron, Intergenic) for each chromosome.
             
Logic:
    1. Loads genomic features into Interval Trees.
    2. Streams DMCs (single positions).
    3. Assigns DMC to a feature based on priority: Exon > Intron > Intergenic.
    4. Normalizes counts by the total length of that feature type.

Author: Paolo Callipo
Date: 2025
Dependencies: python3, pandas, intervaltree
"""

import sys
import argparse
import os
from collections import defaultdict
import pandas as pd

try:
    from intervaltree import Interval, IntervalTree
except ImportError:
    sys.exit("Error: 'intervaltree' library is required. Install via: pip install intervaltree")

def build_interval_trees(bed_file_path):
    print(f"  -> Building interval trees from {os.path.basename(bed_file_path)}...")
    trees = defaultdict(IntervalTree)
    try:
        with open(bed_file_path, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith('#'): continue
                try:
                    parts = line.strip().split('\t')
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    trees[chrom].add(Interval(start, end))
                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        sys.exit(f"Error: File not found: {bed_file_path}")
    return trees

def calculate_feature_lengths(bed_file_path):
    print(f"  -> Calculating feature lengths from {os.path.basename(bed_file_path)}...")
    lengths = defaultdict(int)
    try:
        with open(bed_file_path, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith('#'): continue
                try:
                    parts = line.strip().split('\t')
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    lengths[chrom] += end - start
                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        sys.exit(f"Error: File not found: {bed_file_path}")
    return lengths

def main():
    parser = argparse.ArgumentParser(
        description="Calculate DMC density per genomic feature.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--dmcs", required=True, help="Input DMCs BED file (Chrom, Pos).")
    parser.add_argument("--exons", required=True, help="Exons BED file.")
    parser.add_argument("--introns", required=True, help="Introns BED file.")
    parser.add_argument("--intergenic", required=True, help="Intergenic BED file.")
    parser.add_argument("--fai", required=True, help="Fasta Index (.fai) file to define chromosomes.")
    parser.add_argument("--output", default="dmc_density_stats.tsv", help="Output filename.")

    args = parser.parse_args()

    print("--- Starting DMC Density Calculation ---")

    # 1. Feature Lengths (Denominator)
    exon_lengths = calculate_feature_lengths(args.exons)
    intron_lengths = calculate_feature_lengths(args.introns)
    intergenic_lengths = calculate_feature_lengths(args.intergenic)

    # 2. Build Trees
    exon_trees = build_interval_trees(args.exons)
    intron_trees = build_interval_trees(args.introns)
    intergenic_trees = build_interval_trees(args.intergenic)

    # 3. Stream DMCs (Numerator)
    print(f"  -> Streaming DMCs from {os.path.basename(args.dmcs)}...")
    dmc_counts = defaultdict(lambda: defaultdict(int))
    
    with open(args.dmcs, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            try:
                chrom, start = line.strip().split('\t')[:2]
                pos = int(start)

                # Prioritized Check: Exon > Intron > Intergenic
                if exon_trees[chrom].at(pos):
                    dmc_counts[chrom]['Exon'] += 1
                elif intron_trees[chrom].at(pos):
                    dmc_counts[chrom]['Intron'] += 1
                elif intergenic_trees[chrom].at(pos):
                    dmc_counts[chrom]['Intergenic'] += 1
                
            except (ValueError, IndexError):
                continue

    # 4. Output
    print("  -> Calculating final densities...")
    output_data = []

    try:
        all_chromosomes = [line.split('\t')[0] for line in open(args.fai)]
    except FileNotFoundError:
        sys.exit(f"Error: FAI file '{args.fai}' not found.")

    for chrom in all_chromosomes:
        for category in ["Exon", "Intron", "Intergenic"]:
            length = 0
            count = dmc_counts[chrom][category]

            if category == "Exon": length = exon_lengths[chrom]
            elif category == "Intron": length = intron_lengths[chrom]
            elif category == "Intergenic": length = intergenic_lengths[chrom]

            density = 0.0
            if length > 0:
                density = (count / length) * 100000

            output_data.append({
                'Chromosome': chrom,
                'Category': category,
                'Total_Length_bp': length,
                'DMC_Count': count,
                'Density_per_100kb': density
            })

    final_df = pd.DataFrame(output_data)
    final_df.to_csv(args.output, sep='\t', index=False, float_format='%.4f')
    print(f"Success! Results saved to {args.output}")

if __name__ == "__main__":
    main()
