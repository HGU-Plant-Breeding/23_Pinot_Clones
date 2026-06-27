#!/usr/bin/env python3
"""
Script: 05_calc_mp_rate.py
Description:
    Calculates the rate of MPs (MPs bins / callable bins) per genomic feature
    (Exon, Intron, Intergenic) for each chromosome.

    This replaces the density approach (VMRs per 100kb of feature length) with
    a coverage-aware rate that uses the actual interrogated space as denominator.

    Denominator: all bins in the continuous matrix (passed coverage + missingness filters)
    Numerator:   all bins in the MP binary matrix (passed MP calling filters)

    Both numerator and denominator are assigned to features using majority overlap
    (>50% of the 200bp bin overlaps a feature), with priority fallback:
    Exon > Intron > Intergenic.

Input:
    - Continuous matrix: CG.continuous.tsv (all callable bins)
    - MP binary matrix: CG.vmr.binary.tsv (MP bins only)
    - Exon, Intron, Intergenic BED files
    - FAI file for chromosome list

Output:
    TSV with columns:
    Chromosome, Category, Callable_Bins, VMR_Bins, VMR_Rate

Usage:
    python 05_calc_mp_rate.py --continuous CG.continuous.tsv \
        --vmrs CG.vmr.binary.tsv --exons exons.bed \
        --introns introns.bed --intergenic intergenic.bed \
        --fai genome.fai --output CG.vmr_rate.tsv

Author: Paolo Callipo
Date: 2026
Dependencies: python3, intervaltree
"""

import sys
import argparse
import os
from collections import defaultdict

try:
    from intervaltree import Interval, IntervalTree
except ImportError:
    sys.exit("Error: 'intervaltree' library is required. Install via: pip install intervaltree")


def build_interval_trees(bed_file_path):
    """Build interval trees from a BED file, keyed by chromosome."""
    print(f"  -> Building interval trees from {os.path.basename(bed_file_path)}...")
    trees = defaultdict(IntervalTree)
    try:
        with open(bed_file_path, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                try:
                    parts = line.strip().split('\t')
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    if end > start:
                        trees[chrom].add(Interval(start, end))
                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        sys.exit(f"Error: File not found: {bed_file_path}")
    return trees


def calculate_overlap(tree, chrom, start, end):
    """Calculate total bp overlap between a region and intervals in the tree."""
    if chrom not in tree:
        return 0
    overlaps = tree[chrom].overlap(start, end)
    total = 0
    for iv in overlaps:
        ov_start = max(start, iv.begin)
        ov_end = min(end, iv.end)
        total += max(0, ov_end - ov_start)
    return total


def assign_feature(exon_trees, intron_trees, intergenic_trees, chrom, start, end):
    """
    Assign a 200bp bin to a genomic feature.
    Strategy: majority overlap (>50%), fallback to any overlap (Exon > Intron > Intergenic).
    """
    bin_size = end - start
    half = bin_size / 2.0

    exon_ov = calculate_overlap(exon_trees, chrom, start, end)
    intron_ov = calculate_overlap(intron_trees, chrom, start, end)
    intergenic_ov = calculate_overlap(intergenic_trees, chrom, start, end)

    if exon_ov > half:
        return 'Exon'
    if intron_ov > half:
        return 'Intron'
    if intergenic_ov > half:
        return 'Intergenic'

    if exon_ov > 0:
        return 'Exon'
    if intron_ov > 0:
        return 'Intron'

    return 'Intergenic'


def stream_bins(filepath, exon_trees, intron_trees, intergenic_trees, label):
    """
    Stream a matrix file (continuous or binary), assign each bin to a feature.
    Returns: dict {chrom: {feature: count}}
    """
    print(f"  -> Streaming {label} from {os.path.basename(filepath)}...")
    counts = defaultdict(lambda: defaultdict(int))
    total = 0

    try:
        with open(filepath, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith('chrom'):
                    continue
                try:
                    parts = line.strip().split('\t')
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])

                    feature = assign_feature(exon_trees, intron_trees, intergenic_trees,
                                             chrom, start, end)
                    counts[chrom][feature] += 1
                    total += 1

                except (ValueError, IndexError):
                    continue
    except FileNotFoundError:
        sys.exit(f"Error: File not found: {filepath}")

    print(f"     Total bins: {total:,}")
    for feat in ['Exon', 'Intron', 'Intergenic']:
        feat_total = sum(counts[c][feat] for c in counts)
        print(f"     {feat}: {feat_total:,} ({100*feat_total/max(total,1):.1f}%)")

    return counts, total


def main():
    parser = argparse.ArgumentParser(
        description="Calculate VMR rate (VMR bins / callable bins) per genomic feature.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--continuous", required=True,
                        help="Continuous matrix (all callable bins): CG.continuous.tsv")
    parser.add_argument("--vmrs", required=True,
                        help="VMR binary matrix (VMR bins only): CG.vmr.binary.tsv")
    parser.add_argument("--exons", required=True, help="Exons BED file.")
    parser.add_argument("--introns", required=True, help="Introns BED file.")
    parser.add_argument("--intergenic", required=True, help="Intergenic BED file.")
    parser.add_argument("--fai", required=True, help="Fasta Index (.fai) file for chromosome list.")
    parser.add_argument("--output", default="vmr_rate.tsv", help="Output filename.")

    args = parser.parse_args()

    print("--- Starting VMR Rate Calculation ---")
    print("  Denominator: all callable bins (continuous matrix)")
    print("  Numerator:   VMR bins (binary matrix)\n")

    # Build feature interval trees
    exon_trees = build_interval_trees(args.exons)
    intron_trees = build_interval_trees(args.introns)
    intergenic_trees = build_interval_trees(args.intergenic)

    # Count callable bins per chrom/feature (denominator)
    print("\n--- Counting callable bins (denominator) ---")
    callable_counts, total_callable = stream_bins(
        args.continuous, exon_trees, intron_trees, intergenic_trees, "callable bins"
    )

    # Count VMR bins per chrom/feature (numerator)
    print("\n--- Counting VMR bins (numerator) ---")
    vmr_counts, total_vmrs = stream_bins(
        args.vmrs, exon_trees, intron_trees, intergenic_trees, "VMR bins"
    )

    print(f"\n  Overall VMR rate: {total_vmrs}/{total_callable} = "
          f"{100*total_vmrs/max(total_callable,1):.2f}%")

    # Load chromosome order from FAI
    try:
        all_chromosomes = [line.split('\t')[0] for line in open(args.fai) if line.strip()]
    except FileNotFoundError:
        sys.exit(f"Error: FAI file '{args.fai}' not found.")

    # Write output
    print("\n  -> Writing output...")
    with open(args.output, 'w') as fout:
        fout.write("Chromosome\tCategory\tCallable_Bins\tVMR_Bins\tVMR_Rate\n")

        for chrom in all_chromosomes:
            for category in ["Exon", "Intron", "Intergenic"]:
                callable_n = callable_counts[chrom][category]
                vmr_n = vmr_counts[chrom][category]
                rate = vmr_n / callable_n if callable_n > 0 else 0.0

                fout.write(f"{chrom}\t{category}\t{callable_n}\t{vmr_n}\t{rate:.6f}\n")

    print(f"\nSuccess! Results saved to {args.output}")


if __name__ == "__main__":
    main()
