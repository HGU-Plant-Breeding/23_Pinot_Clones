#!/bin/bash
#==============================================================================
# SCRIPT: run_sniffles_single.sh
# DESCRIPTION: Step 1 of SV Calling. Runs Sniffles on a single BAM to produce 
#              an .snf binary file and a single-sample VCF.
#
# USAGE: sbatch run_sniffles_single.sh <SAMPLE_ID> <BAM_DIR> <REF_FASTA> <OUT_DIR>
#==============================================================================

#SBATCH --job-name=Sniffles_Single
#SBATCH --partition=all
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=2
#SBATCH --time=1:00:00
#SBATCH --output=logs/sniff_%x_%j.out
#SBATCH --error=logs/sniff_%x_%j.err

set -e

# --- MODULES ---
module load sniffles

# --- ARGUMENTS ---
SAMPLE_ID=${1}
IN_DIR=${2}
REFERENCE=${3}
OUT_DIR=${4}

# --- DEFAULTS ---
MIN_SUPPORT=2  # User specific parameter for clonal analysis

# --- VALIDATION ---
if [[ -z "$SAMPLE_ID" || -z "$IN_DIR" || -z "$REFERENCE" || -z "$OUT_DIR" ]]; then
    echo "Usage: sbatch $0 <SAMPLE_ID> <BAM_DIR> <REF_FASTA> <OUT_DIR>"
    exit 1
fi

mkdir -p "$OUT_DIR"
THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "======================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Sample: $SAMPLE_ID"
echo "Reference: $REFERENCE"
echo "Min Support: $MIN_SUPPORT"
echo "======================================================"

echo "Running Sniffles..."

# Note: --snf output is required for the merging step
sniffles \
    --input "${IN_DIR}/${SAMPLE_ID}.bam" \
    --snf "${OUT_DIR}/${SAMPLE_ID}.snf" \
    --vcf "${OUT_DIR}/${SAMPLE_ID}.vcf.gz" \
    --minsupport "$MIN_SUPPORT" \
    --reference "$REFERENCE" \
    --threads "$THREADS" \
    --sample-id "$SAMPLE_ID"

echo "Finished. Created ${OUT_DIR}/${SAMPLE_ID}.snf"
