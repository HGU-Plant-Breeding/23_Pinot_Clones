#!/usr/bin/env python3
"""
Script Name: classify_gene_orthology.py
Description: Takes a list of genes and checks their orthology status using 
             OrthoFinder results. Determines if genes have orthologs in the 
             other haplotype or paralogs in the same haplotype.
Author: Paolo Callipo
Date: 2025
Dependencies: python3
"""

import argparse
import sys
import os

def normalize_gene_id(gene_id: str) -> str:
    """
    Normalizes a gene ID by removing the transcript suffix (e.g., '.1').
    Assumes the suffix is a '.' followed by digits at the very end.
    Example: 'mikado.PN10_HapAG1.1' -> 'mikado.PN10_HapAG1'
    """
    parts = gene_id.rsplit('.', 1)
    if len(parts) > 1 and parts[-1].isdigit():
        return parts[0]
    return gene_id

def main():
    parser = argparse.ArgumentParser(
        description="Classify gene IDs based on OrthoFinder orthogroup information."
    )
    parser.add_argument(
        "gene_list_file",
        help="Path to the input file containing a list of gene IDs (one per line)."
    )
    parser.add_argument(
        "orthogroups_file",
        help="Path to the OrthoFinder Orthogroups.tsv file."
    )

    args = parser.parse_args()

    gene_list_path = args.gene_list_file
    orthogroups_path = args.orthogroups_file

    # --- Step 1: Read the Input Gene List ---
    input_gene_ids = set()
    try:
        with open(gene_list_path, 'r') as f:
            for line in f:
                gene_id = line.strip()
                if gene_id:
                    input_gene_ids.add(gene_id)
        print(f"Loaded {len(input_gene_ids)} gene IDs from {gene_list_path}.")
    except FileNotFoundError:
        sys.exit(f"Error: Gene list file not found at '{gene_list_path}'")

    # --- Step 2: Process Orthogroups.tsv ---
    # Maps normalized gene ID to {'orthogroup_id': str, 'hapA_genes': set, 'hapB_genes': set}
    gene_to_orthogroup_context = {}

    try:
        with open(orthogroups_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')

                if not parts or not parts[0]: # Skip empty lines
                    continue

                orthogroup_id = parts[0]

                # Default parsing logic: 
                # Assumes Col 1 is Haplotype A, Col 2 is Haplotype B (indices 1 and 2)
                # Adjust if your OrthoFinder output columns are different
                hapA_genes_raw = parts[1] if len(parts) > 1 else ""
                hapB_genes_raw = parts[2] if len(parts) > 2 else ""

                current_hapA_genes = set()
                if hapA_genes_raw:
                    for gene in hapA_genes_raw.split(','):
                        current_hapA_genes.add(normalize_gene_id(gene.strip()))

                current_hapB_genes = set()
                if hapB_genes_raw:
                    for gene in hapB_genes_raw.split(','):
                        current_hapB_genes.add(normalize_gene_id(gene.strip()))

                # Store context for all genes in this orthogroup
                all_genes_in_orthogroup = current_hapA_genes.union(current_hapB_genes)
                for gene_id in all_genes_in_orthogroup:
                    gene_to_orthogroup_context[gene_id] = {
                        'orthogroup_id': orthogroup_id,
                        'hapA_genes': current_hapA_genes,
                        'hapB_genes': current_hapB_genes
                    }
        print(f"Processed {line_num} lines from orthogroups file.")
    except FileNotFoundError:
        sys.exit(f"Error: Orthogroups file not found at '{orthogroups_path}'")

    # --- Step 3: Classify Each Gene ---
    total_genes_in_input_list = len(input_gene_ids)
    genes_found_in_orthofinder = 0
    unmatched_genes_count = 0

    ortholog_other_hap_yes = 0
    ortholog_other_hap_no = 0
    paralog_same_hap_yes = 0
    paralog_same_hap_no = 0

    for gene_id in input_gene_ids:
        # Normalize input gene ID to match dictionary keys
        norm_id = normalize_gene_id(gene_id)
        
        if norm_id not in gene_to_orthogroup_context:
            unmatched_genes_count += 1
            continue

        genes_found_in_orthofinder += 1
        context = gene_to_orthogroup_context[norm_id]

        # Determine gene's own haplotype based on naming convention
        # Checks for "HapA" or "HapB" in the ID
        is_hapA_gene = "HapA" in gene_id or "HapAG" in gene_id
        is_hapB_gene = "HapB" in gene_id or "HapBG" in gene_id

        gene_haplotype = None
        if is_hapA_gene and not is_hapB_gene:
            gene_haplotype = 'HapA'
        elif is_hapB_gene and not is_hapA_gene:
            gene_haplotype = 'HapB'
        else:
            # Ambiguous or undetectable haplotype
            unmatched_genes_count += 1
            genes_found_in_orthofinder -= 1
            continue

        own_haplotype_genes = context['hapA_genes'] if gene_haplotype == 'HapA' else context['hapB_genes']
        other_haplotype_genes = context['hapB_genes'] if gene_haplotype == 'HapA' else context['hapA_genes']

        # Question 1: Ortholog in other haplotype?
        if len(other_haplotype_genes) > 0:
            ortholog_other_hap_yes += 1
        else:
            ortholog_other_hap_no += 1

        # Question 2: Paralog in same haplotype?
        # Count > 1 because the gene itself counts as 1
        if len(own_haplotype_genes) > 1:
            paralog_same_hap_yes += 1
        else:
            paralog_same_hap_no += 1

    # --- Step 4: Output Report ---
    print("\n" + "="*50)
    print(" ORTHOFINDER ANALYSIS SUMMARY")
    print("="*50)

    print(f"Total genes in input list: {total_genes_in_input_list:,}")
    print(f"Genes classified:          {genes_found_in_orthofinder:,}")
    print(f"Unmatched/Ambiguous genes: {unmatched_genes_count:,}")

    def print_table(title, yes_count, no_count, total):
        if total == 0: return
        yes_pct = (yes_count / total) * 100
        no_pct = (no_count / total) * 100
        print(f"\n{title}")
        print("-" * 65)
        print(f"{'Answer':<10} | {'Count':>10} | {'Percentage':>10}")
        print("-" * 65)
        print(f"{'Yes':<10} | {yes_count:>10,} | {yes_pct:>9.2f}%")
        print(f"{'No':<10} | {no_count:>10,} | {no_pct:>9.2f}%")
        print("-" * 65)

    if genes_found_in_orthofinder > 0:
        print_table("Has Ortholog in OTHER Haplotype?", ortholog_other_hap_yes, ortholog_other_hap_no, genes_found_in_orthofinder)
        print_table("Has Paralog in SAME Haplotype?", paralog_same_hap_yes, paralog_same_hap_no, genes_found_in_orthofinder)

if __name__ == "__main__":
    main()
