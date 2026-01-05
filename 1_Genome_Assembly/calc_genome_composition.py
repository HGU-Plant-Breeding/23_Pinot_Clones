#!/usr/bin/env python3
"""
Script Name: calc_genome_composition.py
Description: Calculates the total size and percentage composition of genomic 
             regions based on the 4th column of a BED file.
             Useful for summarizing genome partitions (e.g., Hetrozygous vs Hemyzigous vs Homozygous).
Author: Paolo Callipo
Date: 2025
Dependencies: python3
Input Format: BED file (at least 4 columns). 
              Col 2: Start, Col 3: End, Col 4: Region Type (Category).
"""

import sys
import argparse
import os

def calculate_composition(bed_file):
    # Dictionary to store the cumulative size for each region type
    stats = {}
    total_analyzed_size = 0

    print(f"Processing file: {os.path.basename(bed_file)} ...")

    try:
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                # Skip any empty lines or header/comment lines
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('track'):
                    continue

                # Split the line into columns (assuming tab-separated)
                parts = line.split('\t')

                # Ensure the line has enough columns to avoid errors
                if len(parts) < 4:
                    # Optional: Warn user about malformed lines, or just skip
                    continue

                try:
                    # Extract the start, end, and region type
                    # BED is 0-based: Start (col 1), End (col 2), Name/Type (col 3)
                    start = int(parts[1])
                    end = int(parts[2])
                    region_type = parts[3]

                    # Calculate the length of the region
                    length = end - start

                    # Basic validation
                    if length < 0:
                        print(f"Warning: Negative length at line {line_num}. Skipping.", file=sys.stderr)
                        continue

                    # Add the length to the running total for that region type
                    stats[region_type] = stats.get(region_type, 0) + length

                    # Add to the total size
                    total_analyzed_size += length

                except ValueError:
                    print(f"Warning: Non-integer coordinates at line {line_num}. Skipping.", file=sys.stderr)
                    continue

    except FileNotFoundError:
        print(f"ERROR: The file '{bed_file}' was not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

    return stats, total_analyzed_size

def print_report(stats, total_size):
    if total_size == 0:
        print("\nNo data found or file was empty.")
        return

    print("\n" + "="*75)
    print(f" STATISTICS SUMMARY")
    print("="*75)
    print(f"Total Analyzed Size: {total_size:,} bp\n")

    # Print a header for the table
    header = f"{'Region Type':<25} | {'Cumulative Size (bp)':>20} | {'% of Genome':>15}"
    print(header)
    print("-" * len(header))

    # Sort the results by region type for consistent output
    # You could also sort by size: sorted(stats.items(), key=lambda x: x[1], reverse=True)
    for region_type, cumulative_size in sorted(stats.items(), key=lambda x: x[1], reverse=True):
        # Calculate the percentage
        percentage = (cumulative_size / total_size) * 100

        # Print the formatted row
        print(f"{region_type:<25} | {cumulative_size:>20,} | {percentage:>14.2f}%")

    print("-" * len(header))

def main():
    parser = argparse.ArgumentParser(
        description="Calculate genomic composition statistics from a BED file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('input_bed', help="Path to the input BED file. Requires at least 4 columns (Chr, Start, End, Type).")
    
    args = parser.parse_args()

    stats, total_size = calculate_composition(args.input_bed)
    print_report(stats, total_size)

if __name__ == '__main__':
    main()
