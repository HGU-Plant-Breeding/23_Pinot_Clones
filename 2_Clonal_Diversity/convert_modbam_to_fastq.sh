#!/bin/bash
#==============================================================================
# SCRIPT: convert_modbam_to_fastq.sh
# DESCRIPTION: Converts an unaligned BAM (e.g., from Dorado/Guppy) to FASTQ 
#              while PRESERVING methylation tags (MM, ML).
#
# NOTE: This step is crucial for downstream methylation analysis. 
#       Standard bam2fastq might drop these tags.
#
# USAGE: sbatch convert_modbam_to_fastq.sh <INPUT_BAM> <OUTPUT_FASTQ>
#==============================================================================

#SBATCH --job-name=bam2fq_mods
#SBATCH --partition=all
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/bam2fq_%j.out
#SBATCH --error=logs/bam2fq_%j.err

set -e

# --- ARGUMENTS ---
INPUT_BAM=$1
OUTPUT_FASTQ=$2

if [[ -z "$INPUT_BAM" || -z "$OUTPUT_FASTQ" ]]; then
    echo "Usage: sbatch $0 <input.bam> <output.fastq>"
    exit 1
fi

module load samtools

echo "Starting conversion: $INPUT_BAM -> $OUTPUT_FASTQ"
echo "Preserving tags: MM, ML"

# -T MM,ML : Copies the methylation tags from BAM to the FASTQ header
# -@ : Threads
samtools fastq -T MM,ML -@ "${SLURM_CPUS_PER_TASK:-16}" "$INPUT_BAM" > "$OUTPUT_FASTQ"

echo "Conversion complete."
