#!/usr/bin/env python3
"""
Script Name: calc_variant_density.py
Description: Calculates the density of Variants (SNPs, Indels, or SVs) per 
             feature type (Exon, Intron, Intergenic) for each chromosome.

Universal Logic:
    1. Loads genomic features into Interval Trees.
    2. Streams variants (reading Start AND End).
    3. Assigns variant to a feature based on OVERLAP priority:
       Exon > Intron > Intergenic.
       (e.g., If an SV spans an Intron and an Exon, it is counted as Exon).
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
                    chrom, start, end = line.strip().split('\t')[:3]
                    trees[chrom].add(Interval(int(start), int(end)))
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
                    chrom, start, end = line.strip().split('\t')[:3]
                    lengths[chrom] += int(end) - int(start)
                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        sys.exit(f"Error: File not found: {bed_file_path}")
    return lengths

def main():
    parser = argparse.ArgumentParser(
        description="Calculate Variant (SNP/SV) density per genomic feature.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--variants", required=True, help="Input Variant BED file (Chr, Start, End).")
    parser.add_argument("--exons", required=True, help="Exons BED file.")
    parser.add_argument("--introns", required=True, help="Introns BED file.")
    parser.add_argument("--intergenic", required=True, help="Intergenic BED file.")
    parser.add_argument("--fai", required=True, help="Fasta Index (.fai) file to define chromosomes.")
    parser.add_argument("--output", default="variant_density_stats.tsv", help="Output filename.")

    args = parser.parse_args()

    print("--- Starting Variant Density Calculation ---")

    # 1. Feature Lengths (Denominator)
    exon_lengths = calculate_feature_lengths(args.exons)
    intron_lengths = calculate_feature_lengths(args.introns)
    intergenic_lengths = calculate_feature_lengths(args.intergenic)

    # 2. Build Trees
    exon_trees = build_interval_trees(args.exons)
    intron_trees = build_interval_trees(args.introns)
    intergenic_trees = build_interval_trees(args.intergenic)

    # 3. Stream Variants (Numerator)
    print(f"  -> Streaming variants from {os.path.basename(args.variants)}...")
    variant_counts = defaultdict(lambda: defaultdict(int))
    
    with open(args.variants, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'): continue
            try:
                chrom, start, end = line.strip().split('\t')[:3]
                start, end = int(start), int(end)
                
                # Handling single-base inputs (like VCF->BED sometimes produces Start=End)
                if start == end:
                    end += 1

                # Prioritized Overlap Check: Exon > Intron > Intergenic
                # intervaltree.overlap(start, end) returns a set of overlapping intervals
                if exon_trees[chrom].overlap(start, end):
                    variant_counts[chrom]['Exon'] += 1
                elif intron_trees[chrom].overlap(start, end):
                    variant_counts[chrom]['Intron'] += 1
                elif intergenic_trees[chrom].overlap(start, end):
                    variant_counts[chrom]['Intergenic'] += 1
                
            except (ValueError, IndexError):
                continue

    # 4. Output
    print("  -> calculating final densities...")
    output_data = []

    try:
        all_chromosomes = [line.split('\t')[0] for line in open(args.fai)]
    except FileNotFoundError:
        sys.exit(f"Error: FAI file '{args.fai}' not found.")

    for chrom in all_chromosomes:
        for category in ["Exon", "Intron", "Intergenic"]:
            length = 0
            count = variant_counts[chrom][category]

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
                'Variant_Count': count,
                'Density_per_100kb': density
            })

    final_df = pd.DataFrame(output_data)
    final_df.to_csv(args.output, sep='\t', index=False, float_format='%.4f')
    print(f"Success! Results saved to {args.output}")

if __name__ == "__main__":
    main()
