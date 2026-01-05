#!/bin/bash
#==============================================================================
# SCRIPT: run_minimap_ont.sh
# DESCRIPTION: Maps ONT reads to a reference genome using minimap2.
#              Input FASTQ should be generated using `convert_modbam_to_fastq.sh`
#              if methylation analysis is required later.
#
# USAGE: sbatch run_minimap_ont.sh <SAMPLE_ID> <INPUT_DIR> <REF_FASTA> <OUT_DIR>
#==============================================================================

#SBATCH --job-name=minimap
#SBATCH --partition=all
#SBATCH --cpus-per-task=16
#SBATCH --threads-per-core=2
#SBATCH --mem=32G
#SBATCH --output=logs/minimap_%x_%j.out
#SBATCH --error=logs/minimap_%x_%j.err

set -euo pipefail

module load samtools
module load minimap2

SAMPLE_ID=${1}
IN_DIR=${2}
REFERENCE=${3}
OUT_DIR=${4}

# Validation
if [[ -z "$SAMPLE_ID" || -z "$IN_DIR" || -z "$REFERENCE" || -z "$OUT_DIR" ]]; then
    echo "Usage: sbatch $0 <SAMPLE_ID> <INPUT_DIR> <REF_FASTA> <OUT_DIR>"
    exit 1
fi

# Determine input (handle .fastq or .fq.gz)
if [[ -f "${IN_DIR}/${SAMPLE_ID}.fastq" ]]; then
    INPUT_FILE="${IN_DIR}/${SAMPLE_ID}.fastq"
elif [[ -f "${IN_DIR}/${SAMPLE_ID}.fq.gz" ]]; then
    INPUT_FILE="${IN_DIR}/${SAMPLE_ID}.fq.gz"
else
    echo "ERROR: Input FASTQ not found for sample: ${SAMPLE_ID}"
    exit 1
fi

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "Mapping ${SAMPLE_ID} to ${REFERENCE}..."

# Minimap2 mapping -> Samtools view (bam) -> Samtools sort
minimap2 -t "$THREADS" -L -y -Y -ax map-ont "$REFERENCE" "$INPUT_FILE" | \
samtools view -b -@ "$THREADS" - | \
samtools sort -@ "$THREADS" -o "${SAMPLE_ID}.bam" -

samtools index -@ "$THREADS" "${SAMPLE_ID}.bam"

echo "Finished mapping: ${OUT_DIR}/${SAMPLE_ID}.bam"
