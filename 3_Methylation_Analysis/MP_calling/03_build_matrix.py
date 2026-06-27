#!/usr/bin/env python3
"""
Script: 03_build_matrix_v2.py
Description:
    Binarize-first approach for MP calling.
    
    Logic:
      1. Load all sample bins, key by (chrom, bin_start, bin_end)
      2. Apply per-sample coverage filter (min total coverage per bin)
      3. Apply population missingness filter (min samples with valid data)
      4. Binarize each bin per sample using context-specific thresholds:
           < low_thresh  -> 0 (unmethylated)
           > high_thresh -> 1 (methylated)
           in between    -> . (uncertain)
           no data       -> . (missing)
      5. Call MPs: a bin is a MP only if it has BOTH 0s AND 1s across samples
      6. Filter by minimum informative cells (samples with 0 or 1)
      7. Generate distribution plots
      8. Output matrices

    This avoids the problem of calling VMRs on continuous values then losing
    information during binarization. A bin where all samples are 80-100%
    methylated is NOT a VMR because they all agree on state (all = 1).

Usage:
    python 03_build_matrix_v2.py --binned_dir ./02_binned --context CG \\
        --output_prefix ./03_matrix_v2/CG

Author: Paolo Callipo
Date: 2026
"""

import os
import sys
import glob
import argparse
import numpy as np
from collections import defaultdict


def load_sample(filepath, min_cov):
    """
    Load one sample's binned file.
    Returns dict: (chrom, bin_start, bin_end) -> {'nmod': int, 'nvalid': int, 'meth': float}
    Bins below min_cov are returned as None.
    """
    bins = {}
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            try:
                chrom = parts[0]
                bstart = int(parts[1])
                bend = int(parts[2])
                nmod = int(parts[4])
                nvalid = int(parts[5])
                
                key = (chrom, bstart, bend)
                
                if nvalid >= min_cov:
                    meth = (nmod / nvalid) * 100.0
                    bins[key] = {'nmod': nmod, 'nvalid': nvalid, 'meth': meth}
                else:
                    bins[key] = None
                    
            except (ValueError, ZeroDivisionError):
                continue
    return bins


def build_and_binarize(binned_dir, context, min_cov, min_samples,
                       low_thresh, high_thresh, min_informative):
    """
    Load all samples, merge into matrix, binarize, then call VMRs.
    
    Returns:
      - sample_names
      - all_keys: all bins passing coverage + missingness
      - continuous: dict {key: [meth_or_None per sample]}
      - binary: dict {key: [str per sample]}  ('0', '1', '.')
      - vmr_keys: bins with both 0s and 1s AND >= min_informative
      - stats: dict with counts at each stage
    """
    pattern = os.path.join(binned_dir, f"*.{context}.bins.bed")
    files = sorted(glob.glob(pattern))
    
    if not files:
        print(f"ERROR: No files found matching {pattern}", file=sys.stderr)
        sys.exit(1)
    
    sample_names = []
    all_data = []
    all_keys = set()
    
    print(f"Loading {len(files)} samples for {context}...", file=sys.stderr)
    
    for filepath in files:
        fname = os.path.basename(filepath)
        sample_name = fname.replace(f".{context}.bins.bed", "")
        sample_names.append(sample_name)
        
        data = load_sample(filepath, min_cov)
        all_data.append(data)
        all_keys.update(data.keys())
        
        valid = sum(1 for v in data.values() if v is not None)
        print(f"  {sample_name}: {valid:,} valid bins", file=sys.stderr)
    
    total_raw = len(all_keys)
    print(f"\nTotal unique bins across all samples: {total_raw:,}", file=sys.stderr)
    
    # Build continuous matrix + apply missingness filter
    continuous = {}
    filtered_keys = []
    
    for key in sorted(all_keys):
        row = []
        valid_count = 0
        for data in all_data:
            entry = data.get(key)
            if entry is not None:
                row.append(entry['meth'])
                valid_count += 1
            else:
                row.append(None)
        
        if valid_count >= min_samples:
            continuous[key] = row
            filtered_keys.append(key)
    
    after_miss = len(filtered_keys)
    print(f"After missingness filter (>={min_samples}/{len(sample_names)}): "
          f"{after_miss:,}", file=sys.stderr)
    
    # Binarize ALL filtered bins
    print(f"\nBinarizing (low={low_thresh}%, high={high_thresh}%)...", file=sys.stderr)
    binary = {}
    for key in filtered_keys:
        row = []
        for val in continuous[key]:
            if val is None:
                row.append('.')
            elif val < low_thresh:
                row.append('0')
            elif val > high_thresh:
                row.append('1')
            else:
                row.append('.')
        binary[key] = row
    
    # Count binary stats on all filtered bins
    all_zeros = sum(row.count('0') for row in binary.values())
    all_ones = sum(row.count('1') for row in binary.values())
    all_dots = sum(row.count('.') for row in binary.values())
    all_cells = len(filtered_keys) * len(sample_names)
    
    print(f"  All filtered bins binarized:", file=sys.stderr)
    print(f"    0s: {all_zeros:,} ({100*all_zeros/all_cells:.1f}%)", file=sys.stderr)
    print(f"    1s: {all_ones:,} ({100*all_ones/all_cells:.1f}%)", file=sys.stderr)
    print(f"    .s: {all_dots:,} ({100*all_dots/all_cells:.1f}%)", file=sys.stderr)
    
    # Call VMRs: bins with BOTH 0s AND 1s present
    print(f"\nCalling VMRs (bins with both 0 and 1 states)...", file=sys.stderr)
    vmr_candidates = []
    for key in filtered_keys:
        row = binary[key]
        has_zero = '0' in row
        has_one = '1' in row
        if has_zero and has_one:
            vmr_candidates.append(key)
    
    print(f"  VMR candidates (both states present): {len(vmr_candidates):,}", file=sys.stderr)
    
    # Filter by min informative cells
    print(f"  Filtering: >={min_informative} informative cells...", file=sys.stderr)
    vmr_keys = []
    for key in vmr_candidates:
        row = binary[key]
        informative = row.count('0') + row.count('1')
        if informative >= min_informative:
            vmr_keys.append(key)
    
    print(f"  Final VMRs: {len(vmr_keys):,}", file=sys.stderr)
    
    # VMR stats
    if vmr_keys:
        vmr_cells = len(vmr_keys) * len(sample_names)
        vmr_zeros = sum(binary[k].count('0') for k in vmr_keys)
        vmr_ones = sum(binary[k].count('1') for k in vmr_keys)
        vmr_dots = sum(binary[k].count('.') for k in vmr_keys)
        
        print(f"\n  VMR binary matrix stats:", file=sys.stderr)
        print(f"    Total cells:  {vmr_cells:,}", file=sys.stderr)
        print(f"    0 (unmeth):   {vmr_zeros:,} ({100*vmr_zeros/vmr_cells:.1f}%)", file=sys.stderr)
        print(f"    1 (meth):     {vmr_ones:,} ({100*vmr_ones/vmr_cells:.1f}%)", file=sys.stderr)
        print(f"    . (miss/unc): {vmr_dots:,} ({100*vmr_dots/vmr_cells:.1f}%)", file=sys.stderr)
    
    stats = {
        'total_raw': total_raw,
        'after_missingness': after_miss,
        'vmr_candidates': len(vmr_candidates),
        'vmr_final': len(vmr_keys),
    }
    
    return sample_names, filtered_keys, continuous, binary, vmr_keys, stats


def plot_distributions(continuous, binary, filtered_keys, vmr_keys, context, 
                       output_path, low_thresh, high_thresh):
    """
    Generate distribution plots:
      1. Methylation distribution with binarization thresholds
      2. Same in log scale
      3. Per-bin composition: how many 0s vs 1s vs .s per VMR
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plots", file=sys.stderr)
        return
    
    # Collect all methylation values from filtered bins
    all_values = []
    for key in filtered_keys:
        values = [v for v in continuous[key] if v is not None]
        all_values.extend(values)
    
    # Collect MAF for VMRs
    vmr_mafs = []
    for key in vmr_keys:
        row = binary[key]
        n_zeros = row.count('0')
        n_ones = row.count('1')
        n_inform = n_zeros + n_ones
        if n_inform > 0:
            freq = n_ones / n_inform
            maf = min(freq, 1.0 - freq)
            vmr_mafs.append(maf)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'{context} Context — Binarize-First Approach '
                 f'(n={len(filtered_keys):,} bins, {len(vmr_keys):,} VMRs)', fontsize=13)
    
    # Plot 1: Methylation distribution
    ax1 = axes[0]
    ax1.hist(all_values, bins=100, range=(0, 100), color='steelblue', edgecolor='none', alpha=0.8)
    ax1.axvline(low_thresh, color='red', linestyle='--', linewidth=1.5,
                label=f'Unmeth threshold ({low_thresh}%)')
    ax1.axvline(high_thresh, color='darkred', linestyle='--', linewidth=1.5,
                label=f'Meth threshold ({high_thresh}%)')
    ax1.axvspan(low_thresh, high_thresh, alpha=0.1, color='grey', label='Uncertain zone')
    ax1.set_xlabel('Methylation (%)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Bin methylation distribution')
    ax1.legend(fontsize=7)
    
    # Plot 2: Log scale
    ax2 = axes[1]
    ax2.hist(all_values, bins=100, range=(0, 100), color='steelblue', edgecolor='none', alpha=0.8)
    ax2.axvline(low_thresh, color='red', linestyle='--', linewidth=1.5)
    ax2.axvline(high_thresh, color='darkred', linestyle='--', linewidth=1.5)
    ax2.axvspan(low_thresh, high_thresh, alpha=0.1, color='grey')
    ax2.set_xlabel('Methylation (%)')
    ax2.set_ylabel('Frequency (log)')
    ax2.set_title('Bin methylation distribution (log scale)')
    ax2.set_yscale('log')
    
    # Plot 3: MAF distribution of VMRs
    ax3 = axes[2]
    if vmr_mafs:
        # Bin into n_samples/2 bins for SFS-like plot
        ax3.hist(vmr_mafs, bins=50, range=(0, 0.5), color='darkorange', edgecolor='none', alpha=0.8)
    ax3.set_xlabel('Minor Allele Frequency (MAF)')
    ax3.set_ylabel('Number of VMRs')
    ax3.set_title('VMR frequency spectrum')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Distribution plot saved: {output_path}", file=sys.stderr)


def write_matrix(filepath, keys, sample_names, matrix, value_format="continuous"):
    """Write matrix to TSV."""
    with open(filepath, 'w') as f:
        f.write("chrom\tbin_start\tbin_end\t" + "\t".join(sample_names) + "\n")
        
        for key in keys:
            chrom, bstart, bend = key
            row = matrix[key]
            
            if value_format == "continuous":
                vals = []
                for v in row:
                    if v is None:
                        vals.append(".")
                    else:
                        vals.append(f"{v:.2f}")
                f.write(f"{chrom}\t{bstart}\t{bend}\t" + "\t".join(vals) + "\n")
            else:
                f.write(f"{chrom}\t{bstart}\t{bend}\t" + "\t".join(row) + "\n")


def write_stats(filepath, context, sample_names, params, stats):
    """Write filtering statistics."""
    with open(filepath, 'w') as f:
        f.write(f"Pipeline: Binarize-First VMR Calling (v2)\n")
        f.write(f"Context: {context}\n")
        f.write(f"Samples: {len(sample_names)}\n")
        f.write(f"Sample names: {', '.join(sample_names)}\n\n")
        f.write(f"Parameters:\n")
        for k, v in params.items():
            f.write(f"  {k}: {v}\n")
        f.write(f"\nFiltering cascade:\n")
        f.write(f"  1. Total unique bins:            {stats['total_raw']:,}\n")
        f.write(f"  2. After coverage+missingness:   {stats['after_missingness']:,}\n")
        f.write(f"  3. VMR candidates (0+1 present): {stats['vmr_candidates']:,}\n")
        f.write(f"  4. Final VMRs (min informative): {stats['vmr_final']:,}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Binarize-first VMR calling: binarize, then identify bins with "
                    "both methylated and unmethylated states across samples.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--binned_dir", required=True,
                        help="Directory with binned .bins.bed files from 02_bin_methylation.py")
    parser.add_argument("--context", required=True, choices=['CG', 'CHG', 'CHH'],
                        help="Methylation context to process")
    parser.add_argument("--output_prefix", required=True,
                        help="Output prefix (e.g., ./03_matrix_v2/CG)")
    
    # Filtering
    parser.add_argument("--min_cov", type=int, default=5,
                        help="Minimum total coverage per bin per sample")
    parser.add_argument("--min_samples", type=int, default=18,
                        help="Minimum samples with valid data per bin (out of 23)")
    parser.add_argument("--min_informative", type=int, default=21,
                        help="Minimum samples with 0 or 1 (not .) per VMR in binary matrix")
    
    # Binarization thresholds
    parser.add_argument("--low_thresh", type=float, default=None,
                        help="Below this = unmethylated (0). Defaults: CG=30, CHG=25, CHH=5")
    parser.add_argument("--high_thresh", type=float, default=None,
                        help="Above this = methylated (1). Defaults: CG=70, CHG=50, CHH=15")
    
    args = parser.parse_args()
    
    # Context-specific defaults
    defaults = {
        'CG':  {'low_thresh': 30, 'high_thresh': 70},
        'CHG': {'low_thresh': 25, 'high_thresh': 50},
        'CHH': {'low_thresh': 5,  'high_thresh': 15},
    }
    
    ctx_def = defaults[args.context]
    low_thresh = args.low_thresh if args.low_thresh is not None else ctx_def['low_thresh']
    high_thresh = args.high_thresh if args.high_thresh is not None else ctx_def['high_thresh']
    
    # Create output directory
    out_dir = os.path.dirname(args.output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    
    print(f"{'='*60}", file=sys.stderr)
    print(f"  Pipeline:         Binarize-First VMR Calling (v2)", file=sys.stderr)
    print(f"  Context:          {args.context}", file=sys.stderr)
    print(f"  Min coverage:     {args.min_cov}x per bin per sample", file=sys.stderr)
    print(f"  Min samples:      {args.min_samples}/23 (coverage)", file=sys.stderr)
    print(f"  Min informative:  {args.min_informative}/23 (binary)", file=sys.stderr)
    print(f"  Binarize:         <{low_thresh}% = 0, >{high_thresh}% = 1", file=sys.stderr)
    print(f"  VMR definition:   bins with BOTH 0 and 1 states", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)
    
    # Run pipeline
    sample_names, filtered_keys, continuous, binary, vmr_keys, stats = \
        build_and_binarize(args.binned_dir, args.context, args.min_cov, 
                          args.min_samples, low_thresh, high_thresh, 
                          args.min_informative)
    
    # Plots
    print(f"\nGenerating distribution plots...", file=sys.stderr)
    plot_distributions(continuous, binary, filtered_keys, vmr_keys, args.context,
                       f"{args.output_prefix}.distributions.png",
                       low_thresh, high_thresh)
    
    # Write outputs
    print(f"\nWriting output files...", file=sys.stderr)
    
    # Continuous matrix (all filtered bins)
    write_matrix(f"{args.output_prefix}.continuous.tsv",
                 filtered_keys, sample_names, continuous, "continuous")
    
    # VMR continuous (just the VMR bins with continuous values)
    write_matrix(f"{args.output_prefix}.vmr.continuous.tsv",
                 vmr_keys, sample_names, continuous, "continuous")
    
    # VMR binary (the final output for phylogenetics)
    write_matrix(f"{args.output_prefix}.vmr.binary.tsv",
                 vmr_keys, sample_names, binary, "binary")
    
    # Stats
    params = {
        'min_coverage': args.min_cov,
        'min_samples': args.min_samples,
        'min_informative': args.min_informative,
        'low_thresh': low_thresh,
        'high_thresh': high_thresh,
    }
    write_stats(f"{args.output_prefix}.stats.txt", args.context, 
                sample_names, params, stats)
    
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"  Output files:", file=sys.stderr)
    print(f"    {args.output_prefix}.continuous.tsv       (all filtered bins)", file=sys.stderr)
    print(f"    {args.output_prefix}.vmr.continuous.tsv   (VMR bins, continuous)", file=sys.stderr)
    print(f"    {args.output_prefix}.vmr.binary.tsv       (VMR bins, 0/1/.)", file=sys.stderr)
    print(f"    {args.output_prefix}.distributions.png", file=sys.stderr)
    print(f"    {args.output_prefix}.stats.txt", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)


if __name__ == '__main__':
    main()
