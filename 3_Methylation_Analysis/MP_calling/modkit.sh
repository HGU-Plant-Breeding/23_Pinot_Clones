#!/bin/bash
#==============================================================================
# SCRIPT: run_modkit_pileup.sh
# DESCRIPTION: Extracts 5mC methylation calls from aligned BAM files using modkit.
#              1. Generates a genome-wide pileup of all cytosines.
#              2. Splits the output into separate BED files for CG, CHG, and CHH contexts.
#
# USAGE: sbatch run_modkit_pileup.sh <SAMPLE_ID> <INPUT_DIR> <REF_FASTA> <OUT_DIR>
#==============================================================================

#SBATCH --job-name=modkit_pileup
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=2
#SBATCH --mem=32G
#SBATCH --output=logs/modkit_%x_%j.out
#SBATCH --error=logs/modkit_%x_%j.err

# Exit on error
set -euo pipefail

# --- MODULES ---
module load modkit
module load samtools
module load python

# --- ARGUMENTS ---
SAMPLE_ID=${1}
IN_DIR=${2}
REFERENCE=${3}
OUT_DIR=${4}

# --- VALIDATION ---
if [[ -z "$SAMPLE_ID" || -z "$IN_DIR" || -z "$REFERENCE" || -z "$OUT_DIR" ]]; then
    echo "Usage: sbatch $0 <SAMPLE_ID> <INPUT_DIR> <REF_FASTA> <OUT_DIR>"
    exit 1
fi

BAM_FILE="${IN_DIR}/${SAMPLE_ID}.bam"

if [[ ! -f "$BAM_FILE" ]]; then
    echo "Error: Input BAM not found: $BAM_FILE"
    exit 1
fi

# --- SETUP ---
mkdir -p "$OUT_DIR"
THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "======================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Sample: $SAMPLE_ID"
echo "Input BAM: $BAM_FILE"
echo "Output Dir: $OUT_DIR"
echo "======================================================"

# --- STEP 1: MODKIT PILEUP ---
echo "[Step 1] Running modkit pileup..."

# Explanation of flags:
# --motif: Specifies search motifs (CG, CHG, CHH) and offset (0).
# --ignore h: Ignores read haplotype tags (aggregates both haplotypes if unphased).
# --ref: Reference genome required for motif context.

modkit pileup "$BAM_FILE" "${OUT_DIR}/${SAMPLE_ID}_all_contexts.bed" \
    --motif CG 0 \
    --motif CHG 0 \
    --motif CHH 0 \
    --ignore h \
    --ref "$REFERENCE" \
    --threads "$THREADS"

# --- STEP 2: SPLIT CONTEXTS ---
echo "[Step 2] Splitting into context-specific files..."

# Grep is efficient here. We use the pileup output to create 3 separate files.
# The 'modkit' output puts the motif context in a specific column (usually 4 or 10 depending on format, 
# but simply grepping the context string works reliably for standard bedmethyl).

cd "$OUT_DIR"

grep "CG"  "${SAMPLE_ID}_all_contexts.bed" > "${SAMPLE_ID}_CG.bed" &
grep "CHG" "${SAMPLE_ID}_all_contexts.bed" > "${SAMPLE_ID}_CHG.bed" &
grep "CHH" "${SAMPLE_ID}_all_contexts.bed" > "${SAMPLE_ID}_CHH.bed" &

# Wait for parallel grep jobs to finish
wait

# Optional: Compress to save space
# gzip "${SAMPLE_ID}_CG.bed" "${SAMPLE_ID}_CHG.bed" "${SAMPLE_ID}_CHH.bed"

echo "======================================================"
echo "Processing Complete."
echo "Generated:"
echo "  - ${SAMPLE_ID}_CG.bed"
echo "  - ${SAMPLE_ID}_CHG.bed"
echo "  - ${SAMPLE_ID}_CHH.bed"
echo "======================================================"

# --- STEP 2: PROCESS CONTEXTS ---
echo "[Step 2] Processing context files..."

# Ensure the python script is executable or called with python
PYTHON_SCRIPT="/path/to/process_methylation_bed.py"

# Process CG (Merge symmetrical)
python "$PYTHON_SCRIPT" "${SAMPLE_ID}_all_contexts.bed" "${SAMPLE_ID}_CG_merged.bed" --context CG &

# Process CHG (Merge symmetrical)
python "$PYTHON_SCRIPT" "${SAMPLE_ID}_all_contexts.bed" "${SAMPLE_ID}_CHG_merged.bed" --context CHG &

# Process CHH (Filter only)
python "$PYTHON_SCRIPT" "${SAMPLE_ID}_all_contexts.bed" "${SAMPLE_ID}_CHH_filtered.bed" --context CHH &

wait

echo "Processing Complete."
