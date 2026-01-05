#!/usr/bin/env python3
"""
Script Name: parse_variants_for_plotsr.py
Description: Parses a variant file (e.g., modified SyRI output) to generate a
             simplified list of variant positions. This output is specifically
             formatted for use with the modified version of PlotSR's --snp-density 
             feature, allowing visualization of small variant density heatmaps.
Author: Paolo Callipo
Date: 2025
Dependencies: 
    - python3
Input Format Expectation:
    Tab-separated or whitespace-separated file with a header containing
    'TYPE', 'CHR_A', 'START_A', 'CHR_B', 'START_B' columns.
Output Format:
    genome_id\tchromosome\tposition (tab-separated) for each selected variant.
"""

import argparse
import sys
import os

def parse_variants(input_file, output_file, genome_a_id, genome_b_id, variant_types):
    """
    Parses a variant file to extract variant positions for two genomes.
    The output is formatted as 'genome_id\tchromosome\tposition'.

    Args:
        input_file (str): Path to the input variant file.
        output_file (str): Path to write the formatted output.
        genome_a_id (str): Identifier for the first genome (e.g., 'A', 'Hap1').
        genome_b_id (str): Identifier for the second genome (e.g., 'B', 'Hap2').
        variant_types (list): A list of variant types to include (e.g., ['SNP', 'INS']).
    """
    print(f"Starting parsing of '{os.path.basename(input_file)}'.")
    print(f"Including variant types: {', '.join(variant_types)}")

    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            # Read the header line to get column indices, making it robust to column order changes
            header = infile.readline().strip().split()
            try:
                # Find the 0-based index for each required column
                type_idx = header.index('TYPE')
                chra_idx = header.index('CHR_A')
                starta_idx = header.index('START_A')
                chrb_idx = header.index('CHR_B')
                startb_idx = header.index('START_B')
            except ValueError as e:
                print(f"Error: Missing required column in header - {e}. Cannot proceed.", file=sys.stderr)
                print(f"Expected columns are: TYPE, CHR_A, START_A, CHR_B, START_B", file=sys.stderr)
                sys.exit(1)

            # Process the rest of the file
            input_variants_processed = 0 # This will count input variant lines
            for line in infile:
                if not line.strip():
                    continue

                parts = line.strip().split()

                # Ensure the line has enough columns to avoid index errors
                # +1 because index is 0-based, length is 1-based
                required_cols = max(type_idx, chra_idx, starta_idx, chrb_idx, startb_idx) + 1
                if len(parts) < required_cols:
                    print(f"Warning: Skipping malformed line (too few columns): {line.strip()}", file=sys.stderr)
                    continue

                variant_type = parts[type_idx]

                if variant_type in variant_types:
                    try:
                        # Extract data for Genome A
                        chr_a = parts[chra_idx]
                        start_a = parts[starta_idx]

                        # Extract data for Genome B
                        chr_b = parts[chrb_idx]
                        start_b = parts[startb_idx]

                        # Write the formatted lines to the output file
                        outfile.write(f"{genome_a_id}\t{chr_a}\t{start_a}\n")
                        outfile.write(f"{genome_b_id}\t{chr_b}\t{start_b}\n")
                        input_variants_processed += 1
                    except IndexError as e:
                        print(f"Warning: Skipping malformed line (index error): {line.strip()} - {e}", file=sys.stderr)
                        continue
                    except ValueError as e:
                        print(f"Warning: Skipping line with non-integer position: {line.strip()} - {e}", file=sys.stderr)
                        continue


        print(f"Successfully processed {input_variants_processed} variant entries from the input file.")
        print(f"Output written to '{output_file}'.")

    except FileNotFoundError:
        print(f"Error: Input file not found at '{input_file}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    """Main function to handle command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse a variant file to generate a SNP position list compatible with PlotSR's --snp-density feature.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        '--input',
        required=True,
        help="Path to the input variant file (whitespace-separated, typically tab-separated)."
    )
    parser.add_argument(
        '--output',
        required=True,
        help="Path for the output SNP list file. Format: genome_id\\tchromosome\\tposition."
    )
    parser.add_argument(
        '--genome-a-id',
        required=True,
        help="Identifier for the genome in the 'A' columns (e.g., 'A', 'Hap1')."
    )
    parser.add_argument(
        '--genome-b-id',
        required=True,
        help="Identifier for the genome in the 'B' columns (e.g., 'B', 'Hap2')."
    )
    parser.add_argument(
        '--variant-types',
        nargs='+',
        default=['SNP'],
        help="Space-separated list of variant 'TYPE's from the input file to include.\n"
             "Default is 'SNP'. Example: --variant-types SNP INS DEL"
    )

    args = parser.parse_args()

    parse_variants(args.input, args.output, args.genome_a_id, args.genome_b_id, args.variant_types)

if __name__ == '__main__':
    main()
