#!/usr/bin/env python3
"""
Script: 04_calc_vmr_maf.py
Description:
    Calculates the Minor Allele Frequency (MAF) for every MP in the
    binarized methylation matrix (output of 03_build_matrix.py).
    
    MAF = min(freq_methylated, freq_unmethylated) among samples with valid data.
    Output is used to generate Site Frequency Spectrum (SFS) plots.

Input Format:
    Tab-separated matrix from 03_build_matrix.py .binary.tsv:
    chrom  bin_start  bin_end  sample1  sample2  ...
    Values: 0, 1, or . (missing/uncertain)

Output:
    chrom  bin_start  bin_end  n_valid  n_meth  n_unmeth  freq_meth  MAF

Usage:
    python 04_calc_vmr_maf.py CG.binary.tsv CG.maf.tsv

Author: Paolo Callipo
Date: 2026
"""

import sys
import argparse
import numpy as np


def calculate_maf(input_file, output_file):
    print(f"Loading matrix: {input_file}...")
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        # Read header
        header = fin.readline().strip().split('\t')
        sample_names = header[3:]  # skip chrom, bin_start, bin_end
        n_samples = len(sample_names)
        
        print(f"  Samples: {n_samples}")
        
        # Write output header
        fout.write("chrom\tbin_start\tbin_end\tn_valid\tn_meth\tn_unmeth\tfreq_meth\tMAF\n")
        
        total_vmrs = 0
        maf_values = []
        
        for line in fin:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            chrom = parts[0]
            bstart = parts[1]
            bend = parts[2]
            genotypes = parts[3:]
            
            # Count 0s, 1s, and missing
            n_meth = genotypes.count('1')
            n_unmeth = genotypes.count('0')
            n_valid = n_meth + n_unmeth
            
            if n_valid == 0:
                continue
            
            freq_meth = n_meth / n_valid
            maf = min(freq_meth, 1.0 - freq_meth)
            
            fout.write(f"{chrom}\t{bstart}\t{bend}\t{n_valid}\t{n_meth}\t{n_unmeth}\t"
                       f"{freq_meth:.4f}\t{maf:.4f}\n")
            
            maf_values.append(maf)
            total_vmrs += 1
    
    # Summary stats
    if maf_values:
        maf_arr = np.array(maf_values)
        print(f"\n  Total VMRs:     {total_vmrs:,}")
        print(f"  MAF mean:       {maf_arr.mean():.4f}")
        print(f"  MAF median:     {np.median(maf_arr):.4f}")
        print(f"  Singletons:     {np.sum(maf_arr < 0.05):,} "
              f"({100*np.sum(maf_arr < 0.05)/len(maf_arr):.1f}%)")
        print(f"  Common (>0.1):  {np.sum(maf_arr > 0.1):,} "
              f"({100*np.sum(maf_arr > 0.1)/len(maf_arr):.1f}%)")
    
    print(f"\n  Output: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate MAF per VMR from binarized methylation matrix.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_matrix", help="Binarized VMR matrix (.binary.tsv)")
    parser.add_argument("output_file", help="Output MAF file (.maf.tsv)")
    
    args = parser.parse_args()
    calculate_maf(args.input_matrix, args.output_file)


if __name__ == "__main__":
    main()
