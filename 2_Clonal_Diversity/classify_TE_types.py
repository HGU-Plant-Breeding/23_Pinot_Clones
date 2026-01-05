#!/usr/bin/env python3
"""
Script Name: classify_sv_types.py
Description: Classifies Structural Variants (INS/DEL) into three categories based on a hierarchy:
             1. Transposable Elements (TEs): Matches ID against a provided TE annotation list.
             2. Centromeric/Satellite: Matches SV length to known repeat unit multiples.
             3. Unknown: If neither of the above.

             *Note: Centromeric unit lengths are configured for Vitis vinifera (Grapevine).*

Author: Paolo Callipo
Date: 2025
Dependencies: python3
"""

import sys
import argparse
import os
from collections import defaultdict

# --- CONFIGURATION ---
# Lengths of centromeric/satellite repeat units in bp.
# These specific values (107, 79, 135, 187) are typical for Grapevine.
CENTROMERIC_UNITS = [107, 79, 135, 187] 
MAX_MULTIPLE = 20
TOLERANCE = 1
# ---------------------

def load_te_families(filename):
    """
    Loads a two-column file of SV_ID -> TE_Family into a dictionary.
    Expected format: ID <tab> Family
    """
    families = {}
    print(f"Loading TE families from {os.path.basename(filename)}...")
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#') or 'Insertion_ID' in line or 'Deletion_ID' in line:
                    continue
                
                try:
                    parts = line.split('\t')
                    if len(parts) < 2: continue
                    
                    sv_id = parts[0]
                    family = parts[1]
                    
                    # Simplify the family name (e.g., "LTR/Gypsy#LTR" -> "LTR/Gypsy")
                    simple_family = "/".join(family.split('#')[1:]) if '#' in family else family
                    families[sv_id] = simple_family
                except ValueError:
                    continue
    except FileNotFoundError:
        sys.exit(f"Error: TE family file '{filename}' not found.")
        
    return families

def generate_centromeric_sizes(units, max_multiple, tolerance):
    """
    Generates a set of target sizes for centromeric repeats based on
    multiples of the base units +/- tolerance.
    """
    sizes = set()
    for unit in units:
        for i in range(1, max_multiple + 1):
            target = unit * i
            for t in range(-tolerance, tolerance + 1):
                sizes.add(target + t)
    return sizes

def main():
    parser = argparse.ArgumentParser(
        description="Classify SVs into TE, Centromeric, or Unknown categories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("sv_list", help="Input file with SVs. Format: SV_ID <tab> TYPE <tab> LENGTH")
    parser.add_argument("--te-ins", required=True, help="TSV file mapping Insertion IDs to TE families.")
    parser.add_argument("--te-del", required=True, help="TSV file mapping Deletion IDs to TE families.")
    parser.add_argument("--output", default="sv_classification_summary.tsv", help="Output summary filename.")
    
    args = parser.parse_args()

    # --- 1. SETUP ---
    te_ins_families = load_te_families(args.te_ins)
    te_del_families = load_te_families(args.te_del)

    print("Generating centromeric target sizes...")
    centromeric_sizes = generate_centromeric_sizes(CENTROMERIC_UNITS, MAX_MULTIPLE, TOLERANCE)
    print(f"Generated {len(centromeric_sizes)} target sizes based on units: {CENTROMERIC_UNITS}")

    ins_counts = defaultdict(int)
    del_counts = defaultdict(int)

    # --- 2. PROCESS ALL SVs ---
    print(f"Processing SV list from {args.sv_list}...")
    try:
        with open(args.sv_list, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'): continue
                
                parts = line.split('\t')
                if len(parts) < 3: continue

                sv_id, sv_type, sv_len_str = parts[0], parts[1], parts[2]
                
                try:
                    sv_len = abs(int(sv_len_str))
                except ValueError:
                    continue # Skip if length is not an integer

                # Classification Logic
                if sv_type == "INS":
                    if sv_id in te_ins_families:
                        family = te_ins_families[sv_id]
                        ins_counts[family] += 1
                    elif sv_len in centromeric_sizes:
                        ins_counts['Centromeric/Satellite'] += 1
                    else:
                        ins_counts['Unknown'] += 1

                elif sv_type == "DEL":
                    if sv_id in te_del_families:
                        family = te_del_families[sv_id]
                        del_counts[family] += 1
                    elif sv_len in centromeric_sizes:
                        del_counts['Centromeric/Satellite'] += 1
                    else:
                        del_counts['Unknown'] += 1
    except FileNotFoundError:
        sys.exit(f"Error: SV list file '{args.sv_list}' not found.")

    # --- 3. SAVE RESULTS ---
    print(f"Saving results to {args.output}...")
    with open(args.output, "w") as f:
        f.write("SV_Type\tCategory\tCount\n")
        
        # Sort for cleaner output
        for category, count in sorted(ins_counts.items(), key=lambda x: x[1], reverse=True):
            f.write(f"INS\t{category}\t{count}\n")
            
        for category, count in sorted(del_counts.items(), key=lambda x: x[1], reverse=True):
            f.write(f"DEL\t{category}\t{count}\n")

    print("--- Classification Summary ---")
    # Print a quick preview to stdout
    print(f"{'Type':<5} | {'Category':<30} | {'Count'}")
    print("-" * 45)
    for k, v in list(sorted(ins_counts.items(), key=lambda x: x[1], reverse=True))[:5]:
        print(f"{'INS':<5} | {k:<30} | {v}")
    print("...")

if __name__ == "__main__":
    main()
