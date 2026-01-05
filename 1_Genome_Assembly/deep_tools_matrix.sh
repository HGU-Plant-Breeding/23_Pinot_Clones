#!/bin/bash
#==============================================================================
# SCRIPT: run_deeptools_matrix.sh
# DESCRIPTION: SLURM submission script to run deepTools `computeMatrix`.
#              It calculates coverage/methylation over genomic features 
#              (Genes or TEs) separated by genomic partition (Het, Hom, Hemi).
#
# USAGE: sbatch run_deeptools_matrix.sh <CONTEXT> <FEATURE_TYPE>
# EXAMPLE: sbatch run_deeptools_matrix.sh CG genes
#==============================================================================

#SBATCH --job-name=deeptools_matrix
#SBATCH --output=logs/deeptools_%x_%j.out   # Assumes a 'logs' dir exists
#SBATCH --error=logs/deeptools_%x_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Exit immediately if a command exits with a non-zero status
set -e

#==============================================================================
# CONFIGURATION
#==============================================================================

# Input directory (Current directory by default)
DATA_DIR="."

# Environment Setup
# ACTION: Adjust these lines to match your HPC environment
echo "Setting up environment..."
# module load anaconda3
# source activate deeptools_env
# Or if using module directly:
# module load deeptools

# Check dependency
if ! command -v computeMatrix &> /dev/null; then
    echo "Error: 'computeMatrix' (deepTools) not found. Please check your environment." >&2
    exit 1
fi

#==============================================================================
# ARGUMENT PARSING
#==============================================================================

# Get arguments
CON=${1}          # e.g., CG, CHG, CHH
FEATURE_TYPE=${2} # e.g., genes, tes

# Validation
if [[ -z "$CON" || -z "$FEATURE_TYPE" ]]; then
    echo "ERROR: Missing arguments." >&2
    echo "Usage: sbatch $0 <CONTEXT> <FEATURE_TYPE>" >&2
    exit 1
fi

# Validating feature type
if [[ "$FEATURE_TYPE" != "genes" && "$FEATURE_TYPE" != "tes" ]]; then
    echo "ERROR: Invalid feature type '${FEATURE_TYPE}'. Must be 'genes' or 'tes'." >&2
    exit 1
fi

# Update Job Name dynamically (if running via SLURM)
if [[ -n "$SLURM_JOB_ID" ]]; then
    scontrol update JobId=${SLURM_JOB_ID} JobName=dt_${FEATURE_TYPE}_${CON}
fi

echo "========================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "Context: ${CON}"
echo "Feature: ${FEATURE_TYPE}"
echo "========================================================"

#==============================================================================
# FILE SETUP
#==============================================================================

# Define Input/Output filenames
# Note: This assumes specific naming conventions in your data folder.
BIGWIG_FILE="${DATA_DIR}/Diploid_${CON}.bw"
OUTPUT_MATRIX="${DATA_DIR}/${FEATURE_TYPE}_metaplot_${CON}_matrix.gz"
OUTPUT_TAB="${DATA_DIR}/${FEATURE_TYPE}_metaplot_${CON}_data.tsv"

# Select Region Files based on Feature Type
if [[ "$FEATURE_TYPE" == "genes" ]]; then
    HET_REGIONS="${DATA_DIR}/diploid_heterozygous_genes.bed"
    HOMO_REGIONS="${DATA_DIR}/diploid_homozygous_genes.bed"
    HEMI_REGIONS="${DATA_DIR}/diploid_hemizygous_genes.bed"
elif [[ "$FEATURE_TYPE" == "tes" ]]; then
    HET_REGIONS="${DATA_DIR}/diploid_heterozygous_tes.bed"
    HOMO_REGIONS="${DATA_DIR}/diploid_homozygous_tes.bed"
    HEMI_REGIONS="${DATA_DIR}/diploid_hemizygous_tes.bed"
fi

# Verify input files exist
for f in "$BIGWIG_FILE" "$HET_REGIONS" "$HOMO_REGIONS" "$HEMI_REGIONS"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Input file not found: $f" >&2
        exit 1
    fi
done

#==============================================================================
# RUN COMPUTE MATRIX
#==============================================================================

echo "Running computeMatrix scale-regions..."

computeMatrix scale-regions \
    -S "${BIGWIG_FILE}" \
    -R "${HET_REGIONS}" "${HOMO_REGIONS}" "${HEMI_REGIONS}" \
    --beforeRegionStartLength 2000 \
    --regionBodyLength 3000 \
    --afterRegionStartLength 2000 \
    --skipZeros \
    --missingDataAsZero \
    --smartLabels \
    --samplesLabel "${CON}_Methylation" \
    -o "${OUTPUT_MATRIX}" \
    --outFileNameMatrix "${OUTPUT_TAB}" \
    -p "${SLURM_CPUS_PER_TASK}"

echo "========================================================"
echo "Analysis Complete."
echo "Matrix saved to: ${OUTPUT_MATRIX}"
echo "Raw data saved to: ${OUTPUT_TAB}"
echo "End Time: $(date)"
echo "========================================================"
