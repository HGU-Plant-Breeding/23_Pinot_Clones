#!/usr/bin/env python3
"""
Script Name: classify_syri_genome.py
Description: Classifies genomic regions into Hemizygous, Heterozygous, and Homozygous 
             states based on SyRI output and variant density.
Author: Paolo Callipo
Date: 2025
Dependencies: 
    - python3 (intervaltree)
    - bedtools
    - sort
"""

import argparse
import sys
import os
import subprocess
import shutil
from collections import defaultdict

# check for external python dependencies
try:
    from intervaltree import Interval, IntervalTree
except ImportError:
    print("Error: The 'intervaltree' library is required.\n"
          "Please install it: pip install intervaltree", file=sys.stderr)
    sys.exit(1)


def run_command(command, check=True, shell=False):
    """
    Helper to run shell commands safely with error handling.
    """
    try:
        # Force standard sorting locale to prevent bedtools sorting errors
        env = os.environ.copy()
        env['LC_ALL'] = 'C' 
        
        result = subprocess.run(
            command,
            check=check,
            shell=shell,
            capture_output=True,
            text=True,
            env=env
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        cmd_str = ' '.join(e.cmd) if isinstance(e.cmd, list) else e.cmd
        print(f"Error executing command: {cmd_str}", file=sys.stderr)
        print(f"Stderr: {e.stderr}", file=sys.stderr)
        sys.exit(1)


def parse_syri_to_bed_files(syri_file, intermediate_dir, hemi_threshold, target_genome):
    """
    Parses syri.out, separating large structural/unique events (Hemizygous)
    from small variants (for Density calculation).
    """
    print(f"[Step 1/5] Parsing {os.path.basename(syri_file)} for Target Genome '{target_genome}'...")

    # SyRI columns: RefChr(0), RefStart(1), RefEnd(2) | QryChr(3), QryStart(4), QryEnd(5)
    if target_genome == 'A':
        chr_col, start_col, end_col = 0, 1, 2
    elif target_genome == 'B':
        chr_col, start_col, end_col = 3, 4, 5
    else:
        # Should be caught by argparse, but double checking
        sys.exit("Error: Invalid target genome.")

    paths = {
        'hemi': os.path.join(intermediate_dir, 'hemizygous_features.bed'),
        'density': os.path.join(intermediate_dir, 'small_variants_for_density.bed')
    }

    with open(syri_file, 'r') as infile, \
         open(paths['hemi'], 'w') as hemi_out, \
         open(paths['density'], 'w') as density_out:

        for line in infile:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue

            # Ensure coordinates exist for the selected genome (A or B)
            if parts[start_col].isdigit():
                anno_type = parts[10]
                chrom = parts[chr_col]
                # SyRI is 1-based, BED is 0-based. Convert Start.
                start = int(parts[start_col]) - 1
                end = int(parts[end_col])
                length = end - start
                
                bed_line = f"{chrom}\t{start}\t{end}\n"

                # Logic:
                # Hemizygous = NOTAL (Not Aligned) OR Large Indels/HDRs
                # Density variants = SNPs OR Small Indels/HDRs
                
                is_large_sv = (anno_type in ['HDR', 'INS', 'DEL'] and length >= hemi_threshold)
                is_small_sv = (anno_type in ['HDR', 'INS', 'DEL'] and length < hemi_threshold)

                if anno_type == 'NOTAL' or is_large_sv:
                    hemi_out.write(bed_line)
                elif anno_type == 'SNP' or is_small_sv:
                    density_out.write(bed_line)
                    
    return paths


def classify_genome(args):
    """
    Main workflow logic.
    """
    # Create output directory structures
    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    output_basename = os.path.splitext(os.path.basename(args.output))[0]
    intermediate_dir = os.path.join(output_dir, f"{output_basename}_tmp_files")
    os.makedirs(intermediate_dir, exist_ok=True)

    print(f"\n--- Processing started ---")
    if args.keep_intermediate:
        print(f"Intermediate files stored in: {intermediate_dir}")

    # 1. Parse SyRI
    bed_paths = parse_syri_to_bed_files(args.syri_out, intermediate_dir, args.hemi_threshold, args.target_genome)

    # 2. Define Hemizygous Regions
    print("[Step 2/5] Defining hemizygous regions...")
    genome_bed = os.path.join(intermediate_dir, 'genome.bed')
    
    # Create a BED representing the whole genome size
    with open(args.chr_sizes, 'r') as infile, open(genome_bed, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                outfile.write(f"{parts[0]}\t0\t{parts[1]}\n")

    final_hemizygous_bed = os.path.join(intermediate_dir, 'final_hemizygous.bed')
    
    # Merge overlapping hemizygous features
    # Note: quoted paths ("{path}") handle spaces in filenames safely
    run_command(f"sort -k1,1 -k2,2n \"{bed_paths['hemi']}\" | bedtools merge -i stdin > \"{final_hemizygous_bed}\"", shell=True)

    # Subtract hemizygous regions from genome to get analyzable regions
    analyzable_bed = os.path.join(intermediate_dir, 'analyzable.bed')
    run_command(f"sort -k1,1 -k2,2n \"{genome_bed}\" | bedtools subtract -a stdin -b \"{final_hemizygous_bed}\" > \"{analyzable_bed}\"", shell=True)

    # 3. Create Windows
    print(f"[Step 3/5] Creating {args.window_size/1000:.1f}kb windows...")
    windows_bed = os.path.join(intermediate_dir, 'windows.bed')
    run_command(f"bedtools makewindows -b \"{analyzable_bed}\" -w {args.window_size} > \"{windows_bed}\"", shell=True)

    # 4. Calculate Density & Classify
    print("\n[Step 4/5] Building variant index and classifying windows...")

    variant_trees = defaultdict(IntervalTree)
    variant_count_total = 0
    
    print("--> Loading small variants into memory...")
    with open(bed_paths['density'], 'r') as f:
        for line in f:
            try:
                chrom, start, end = line.strip().split('\t')[:3]
                variant_trees[chrom].add(Interval(int(start), int(end)))
                variant_count_total += 1
            except (ValueError, IndexError):
                continue
    print(f"--> Loaded {variant_count_total} variants.")

    final_homozygous_bed = os.path.join(intermediate_dir, 'final_homozygous.bed')
    final_heterozygous_bed = os.path.join(intermediate_dir, 'final_heterozygous.bed')

    windows_processed = 0
    
    with open(windows_bed, 'r') as infile, \
         open(final_homozygous_bed, 'w') as hom_out, \
         open(final_heterozygous_bed, 'w') as het_out:

        for line in infile:
            try:
                chrom, start, end = line.strip().split('\t')[:3]
                start, end = int(start), int(end)
                windows_processed += 1
                
                # Count variants intersecting this window
                overlapping_variants = variant_trees[chrom][start:end]
                variant_count = len(overlapping_variants)
                
                # Normalize density (variants per bp * window_size)
                window_len = end - start
                if window_len > 0:
                    normalized_density = (variant_count / window_len) * args.window_size
                else:
                    normalized_density = 0
                
                bed_line = f"{chrom}\t{start}\t{end}\n"
                
                if normalized_density <= args.hom_threshold:
                    hom_out.write(bed_line)
                else:
                    het_out.write(bed_line)
            except (ValueError, IndexError):
                continue

    print(f"--> Processed {windows_processed} windows.")

    # 5. Merge and Output
    print("\n[Step 5/5] Merging results and writing final BED file...")
    merged_hom_bed = os.path.join(intermediate_dir, 'merged_hom.bed')
    merged_het_bed = os.path.join(intermediate_dir, 'merged_het.bed')

    run_command(f"sort -k1,1 -k2,2n \"{final_homozygous_bed}\" | bedtools merge -i stdin > \"{merged_hom_bed}\"", shell=True)
    run_command(f"sort -k1,1 -k2,2n \"{final_heterozygous_bed}\" | bedtools merge -i stdin > \"{merged_het_bed}\"", shell=True)

    # Write final colored BED for IGV/UCSC
    with open(args.output, 'w') as outfile:
        outfile.write('track name="Genome Classification" description="Genomic regions by allelic state" itemRgb="On"\n')

        # RGB Colors: Orange (Hemi), Blue (Het), Grey (Hom)
        color_map = {
            'Hemizygous': '230,126,34', 
            'Heterozygous': '52,152,219', 
            'Homozygous': '155,155,155'
        }
        final_files = {
            'Hemizygous': final_hemizygous_bed, 
            'Homozygous': merged_hom_bed, 
            'Heterozygous': merged_het_bed
        }

        for classification, file_path in final_files.items():
             if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                with open(file_path, 'r') as infile:
                    for line in infile:
                        if not line.strip(): continue
                        parts = line.strip().split('\t')
                        chrom, start, end = parts[0], parts[1], parts[2]
                        color = color_map[classification]
                        # Standard BED9 format
                        outfile.write(f"{chrom}\t{start}\t{end}\t{classification}\t0\t.\t{start}\t{end}\t{color}\n")

    # Clean up
    if not args.keep_intermediate:
        print("Cleaning up intermediate files...")
        shutil.rmtree(intermediate_dir)
    else:
        print("Skipping cleanup (--keep-intermediate was set).")

    print(f"\nDone. The final BED file is at: {args.output}")


def main():
    parser = argparse.ArgumentParser(
        description="Classify a genome into Hemizygous, Heterozygous, and Homozygous regions based on SyRI output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/Output
    parser.add_argument('--syri-out', required=True, 
                        help="Path to the input syri.out file.")
    parser.add_argument('--chr-sizes', required=True, 
                        help="Tab-separated file of chromosome names and sizes. MUST match the target genome.")
    parser.add_argument('--output', required=True, 
                        help="Path for the final output BED file.")
    
    # Parameters
    parser.add_argument('--window-size', type=int, default=10000, 
                        help="Window size (bp) for density calculations.")
    parser.add_argument('--hemi-threshold', type=int, default=1000, 
                        help="Min size (bp) for INDELs/HDRs to be considered Hemizygous.")
    parser.add_argument('--hom-threshold', type=int, default=5, 
                        help="Max variants per window (normalized) to be called Homozygous.")
    parser.add_argument('--target-genome', choices=['A', 'B'], default='A', 
                        help="Analyze Genome 'A' (Ref, cols 1-3) or 'B' (Query, cols 4-6).")
    
    # Flags
    parser.add_argument('--keep-intermediate', action='store_true', 
                        help="If set, temporary files will not be deleted.")

    args = parser.parse_args()

    # Pre-flight check for external tools
    if not shutil.which("bedtools"):
        sys.exit("Error: 'bedtools' not found in PATH. Please install it (conda install bedtools).")
    if not shutil.which("sort"):
        sys.exit("Error: 'sort' command not found.")

    classify_genome(args)

if __name__ == '__main__':
    main()
