#!/usr/bin/env python3
"""
Script: 02_bin_methylation.py
Description:
    Takes the output of 01_parse_modkit.py (chr, start, Nmod, Nvalid per site)
    and aggregates sites into fixed-size non-overlapping genomic bins.

    For each bin, it sums:
      - Total Nmod (methylated reads across all cytosines in the bin)
      - Total Nvalid (total coverage across all cytosines in the bin)
      - Number of cytosines in the bin (Nsites)

    Bin methylation = (Sum_Nmod / Sum_Nvalid) * 100

    Output: chr <tab> bin_start <tab> bin_end <tab> Nsites <tab> Sum_Nmod <tab> Sum_Nvalid <tab> methylation%

Usage:
    python 02_bin_methylation.py sample.CG.bed chrom.sizes output.CG.bins.bed --bin_size 100

Author: Paolo Callipo
Date: 2026
"""

import sys
import argparse


def load_chrom_sizes(chrom_sizes_file):
    """Load chromosome sizes from a .fai or 2-column chrom.sizes file."""
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def process_bins(input_file, chrom_sizes_file, output_file, bin_size, min_sites):
    """
    Streaming binning: reads sorted input, aggregates into bins.
    Assumes input is sorted by chrom, start (which modkit output is).
    """
    chrom_sizes = load_chrom_sizes(chrom_sizes_file)
    
    total_bins = 0
    written_bins = 0
    total_sites = 0
    
    # Current bin state
    curr_chrom = None
    curr_bin_start = None
    curr_bin_end = None
    curr_nmod = 0
    curr_nvalid = 0
    curr_nsites = 0
    
    def write_bin(fout):
        nonlocal written_bins
        if curr_nsites >= min_sites and curr_nvalid > 0:
            meth_pct = (curr_nmod / curr_nvalid) * 100.0
            fout.write(f"{curr_chrom}\t{curr_bin_start}\t{curr_bin_end}\t"
                       f"{curr_nsites}\t{curr_nmod}\t{curr_nvalid}\t{meth_pct:.2f}\n")
            written_bins += 1
    
    def reset_bin(chrom, pos):
        nonlocal curr_chrom, curr_bin_start, curr_bin_end, curr_nmod, curr_nvalid, curr_nsites
        bin_idx = pos // bin_size
        curr_chrom = chrom
        curr_bin_start = bin_idx * bin_size
        curr_bin_end = min(curr_bin_start + bin_size, chrom_sizes.get(chrom, curr_bin_start + bin_size))
        curr_nmod = 0
        curr_nvalid = 0
        curr_nsites = 0
    
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        # Write header
        fout.write("#chrom\tbin_start\tbin_end\tNsites\tNmod\tNvalid\tmeth_pct\n")
        
        for line in fin:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            try:
                chrom = parts[0]
                pos = int(parts[1])
                nmod = int(parts[2])
                nvalid = int(parts[3])
            except ValueError:
                continue
            
            total_sites += 1
            bin_idx = pos // bin_size
            bin_start = bin_idx * bin_size
            
            # Check if we're still in the same bin
            if chrom == curr_chrom and bin_start == curr_bin_start:
                # Accumulate
                curr_nmod += nmod
                curr_nvalid += nvalid
                curr_nsites += 1
            else:
                # Write previous bin (if any)
                if curr_chrom is not None:
                    total_bins += 1
                    write_bin(fout)
                
                # Start new bin
                reset_bin(chrom, pos)
                curr_nmod = nmod
                curr_nvalid = nvalid
                curr_nsites = 1
        
        # Write last bin
        if curr_chrom is not None:
            total_bins += 1
            write_bin(fout)
    
    # Report
    print(f"\nDone binning {input_file}", file=sys.stderr)
    print(f"  Bin size:       {bin_size} bp", file=sys.stderr)
    print(f"  Min sites/bin:  {min_sites}", file=sys.stderr)
    print(f"  Total sites:    {total_sites:,}", file=sys.stderr)
    print(f"  Total bins:     {total_bins:,}", file=sys.stderr)
    print(f"  Written bins:   {written_bins:,} (passed min_sites filter)", file=sys.stderr)
    print(f"  Filtered bins:  {total_bins - written_bins:,}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Bin methylation sites into fixed genomic windows, summing raw read counts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_file", 
                        help="Parsed methylation file from 01_parse_modkit.py (chr, start, Nmod, Nvalid).")
    parser.add_argument("chrom_sizes", 
                        help="Chromosome sizes file (.fai or 2-col: chrom<tab>size).")
    parser.add_argument("output_file", 
                        help="Output binned methylation file.")
    parser.add_argument("--bin_size", type=int, default=100,
                        help="Bin size in bp.")
    parser.add_argument("--min_sites", type=int, default=1,
                        help="Minimum number of cytosines in a bin to keep it. "
                             "Recommended: 3 for CG/CHG, 5 for CHH.")

    args = parser.parse_args()
    process_bins(args.input_file, args.chrom_sizes, args.output_file, 
                 args.bin_size, args.min_sites)


if __name__ == '__main__':
    main()
