#!/bin/bash

# Script Name: bin_methylation_diploid.sh
# Description: Calculates methylation bins (Low/Medium/High) using diploid methylation files.
#              Segregates data by haplotype and genomic feature (Genes, TEs, Genome).
# Input Requirement: 
#   - Methylation files named: Diploid_CG.bed, Diploid_CHG.bed, Diploid_CHH.bed
#   - Feature files named: Hap{A/B}_genes.bed, Hap{A/B}_te.bed
#   - "Genome" mode assumes chromosome names contain "HapA" or "HapB" string.

# Exit on error
set -euo pipefail

# Check for dependencies
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools could not be found."
    exit 1
fi

# --- Arguments ---
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_csv>"
    echo "Example: $0 ./data methylation_summary.csv"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_CSV="$2"

echo "--- Starting Methylation Binning Analysis ---"
echo "Input Directory: $INPUT_DIR"
echo "Output File:     $OUTPUT_CSV"

# Write the header
echo "Haplotype,Region,Context,Low,Medium,High" > "$OUTPUT_CSV"

# Define arrays to loop through
haplotypes=("A" "B")
contexts=("CG" "CHG" "CHH")
regions=("Genes" "TE" "Genome")

# Awk logic for binning
# Low <= 30, Medium 30-70, High > 70
awk_binner='
BEGIN { low=0; med=0; high=0 }
{
    meth_val = $4;
    if (meth_val <= 30) { low++ }
    else if (meth_val > 30 && meth_val <= 70) { med++ }
    else { high++ }
}
END {
    total = low + med + high;
    if (total > 0) {
        low_perc = (low / total) * 100;
        med_perc = (med / total) * 100;
        high_perc = (high / total) * 100;
        printf "%.2f,%.2f,%.2f\n", low_perc, med_perc, high_perc;
    } else {
        printf "0.00,0.00,0.00\n";
    }
}'

# Main loop
for hap in "${haplotypes[@]}"; do
    for ctx in "${contexts[@]}"; do
        
        # Define Input File Path
        diploid_methylation_file="${INPUT_DIR}/Diploid_${ctx}.bed"

        if [ ! -f "$diploid_methylation_file" ]; then
            echo "Warning: File not found: $diploid_methylation_file. Skipping..."
            continue
        fi

        for region in "${regions[@]}"; do
            echo "Processing: Haplotype ${hap} | Context ${ctx} | Region ${region}..."

            # Logic for Genome vs Specific Features
            if [ "$region" == "Genome" ]; then
                # Filter diploid file for lines containing "HapA" or "HapB" in chromosome name
                percentages=$(grep "Hap${hap}" "$diploid_methylation_file" | awk -F'\t' "$awk_binner")
            else
                # Define Feature File Path (e.g., HapA_genes.bed)
                # ${region,,} converts "Genes" to "genes" (lowercase)
                region_file="${INPUT_DIR}/Hap${hap}_${region,,}.bed"
                
                if [ ! -f "$region_file" ]; then
                    echo "  Warning: Feature file not found: $region_file. Outputting 0s."
                    percentages="0.00,0.00,0.00"
                else
                    # Intersect and calculate stats
                    percentages=$(bedtools intersect -a "$diploid_methylation_file" -b "$region_file" -wa | awk -F'\t' "$awk_binner")
                fi
            fi

            # Append to CSV
            echo "Hap${hap},${region},${ctx},${percentages}" >> "$OUTPUT_CSV"
        done
    done
done

echo "--- Analysis Complete. Results saved to ${OUTPUT_CSV} ---"
