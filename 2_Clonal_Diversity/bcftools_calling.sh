#!/bin/bash
#==============================================================================
# SCRIPT: run_bcftools_calling.sh
# DESCRIPTION: Runs bcftools mpileup and call on a specific genomic region.
#              Optimized for ONT data using the '-X ont' configuration.
#
# USAGE: sbatch run_bcftools_calling.sh <REGION> <BAM_LIST> <REF_FASTA> <OUT_DIR>
# EXAMPLE: sbatch run_bcftools_calling.sh Chr1:1-1000000 bam_list.txt ref.fa ./vcf_out
#==============================================================================

#SBATCH --job-name=bcf_SNP
#SBATCH --partition=standard
#SBATCH --cpus-per-task=32
#SBATCH --threads-per-core=2
#SBATCH --time=24:00:00
#SBATCH --output=logs/bcf_%x_%j.out
#SBATCH --error=logs/bcf_%x_%j.err

# Exit on error
set -euo pipefail

# --- MODULES ---
module load java
module load bcftools

# --- ARGUMENTS ---
REGION=${1}        # e.g., "Chr1" or "Chr1:1000-2000"
BAM_LIST=${2}      # Path to a text file containing the list of BAM files
REFERENCE=${3}     # Path to reference FASTA
OUT_DIR=${4}       # Output directory

# --- VALIDATION ---
if [[ -z "$REGION" || -z "$BAM_LIST" || -z "$REFERENCE" || -z "$OUT_DIR" ]]; then
    echo "ERROR: Missing arguments."
    echo "Usage: sbatch $0 <REGION> <BAM_LIST> <REF_FASTA> <OUT_DIR>"
    exit 1
fi

if [[ ! -f "$BAM_LIST" ]]; then
    echo "ERROR: BAM list file not found: $BAM_LIST"
    exit 1
fi

# --- SETUP ---
mkdir -p "$OUT_DIR"
THREADS=${SLURM_CPUS_PER_TASK:-32}

echo "======================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Region: $REGION"
echo "BAM List: $BAM_LIST"
echo "Reference: $REFERENCE"
echo "======================================================"

# --- VARIANT CALLING ---
# Flags explanation:
# -f: Reference genome
# -b: List of input BAM files
# -Ou: Output uncompressed BCF (efficient piping)
# -X ont: Apply Oxford Nanopore specific bias corrections
# -a: Annotate output (AD=Allelic Depth, DP=Depth, etc.)
# -r: Specific region to process
# call -m: Multiallelic calling model
# call -v: Output variant sites only

echo "Starting mpileup | call..."

bcftools mpileup \
    -f "$REFERENCE" \
    -b "$BAM_LIST" \
    --threads "$THREADS" \
    -Ou \
    -X ont \
    -a AD,DP,SP,SCR \
    -r "$REGION" | \
bcftools call \
    -mv \
    -Oz \
    -o "${OUT_DIR}/${REGION}.vcf.gz" \
    --threads "$THREADS"

# Index the output for downstream processing
echo "Indexing VCF..."
bcftools index "${OUT_DIR}/${REGION}.vcf.gz"

echo "======================================================"
echo "Finished. Output: ${OUT_DIR}/${REGION}.vcf.gz"
echo "======================================================"
