#!/bin/bash
#==============================================================================
# SCRIPT: filter_and_merge_vcfs.sh
# DESCRIPTION: Applies hard filters to raw VCF files and merges them.
#              Filters cover Quality, Depth, Strand Bias, and Paralog artifacts.
#
# FILTERS APPLIED:
#   - Low Quality (QUAL <= 10, MQ < 20)
#   - Depth outliers (Default: < 115 or > 690)
#   - Strand/Position Biases (MQBZ, BQBZ, RPBZ, SCBZ)
#   - Fixed Heterozygosity (Removes sites Het in ALL samples)
#
# USAGE: sbatch filter_and_merge_vcfs.sh <IN_DIR> <OUT_DIR> [MIN_DP] [MAX_DP]
# EXAMPLE: sbatch filter_and_merge_vcfs.sh ./raw_vcf ./filtered_vcf 115 690
#==============================================================================

#SBATCH --job-name=bcf_filter
#SBATCH --partition=standard
#SBATCH --cpus-per-task=48
#SBATCH --threads-per-core=2
#SBATCH --time=24:00:00
#SBATCH --output=logs/filter_%x_%j.out
#SBATCH --error=logs/filter_%x_%j.err

set -e

# --- MODULES ---
module load java
module load bcftools

# --- ARGUMENTS ---
IN_DIR=${1}
OUT_DIR=${2}
MIN_DP=${3:-115}  # Default to 115 if not specified
MAX_DP=${4:-690}  # Default to 690 if not specified

# --- VALIDATION ---
if [[ -z "$IN_DIR" || -z "$OUT_DIR" ]]; then
    echo "Usage: sbatch $0 <IN_DIR> <OUT_DIR> [MIN_DP] [MAX_DP]"
    exit 1
fi

mkdir -p "$OUT_DIR"
THREADS=${SLURM_CPUS_PER_TASK:-48}

echo "======================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Input Dir: $IN_DIR"
echo "Output Dir: $OUT_DIR"
echo "Depth Filter: $MIN_DP < DP < $MAX_DP"
echo "======================================================"

# --- STEP 1: FILTERING ---
echo "Starting filtering..."

# Iterate over all VCFs in input directory
# We assume inputs are .vcf.gz
for vcf_file in "${IN_DIR}"/*.vcf.gz; do
    filename=$(basename "$vcf_file")
    base="${filename%.vcf.gz}"
    
    echo "Filtering $filename..."
    
    # Explanation of Filters:
    # QUAL/MQ: Basic quality filtering
    # INFO/DP: Depth filtering (dataset specific)
    # *BZ: Z-score tests for bias (Mapping Quality, Base Quality, Read Position)
    # SCBZ: Soft-Clip Bias Z-score
    # COUNT(GT="het")=23: Removes sites that are Het in ALL 23 clones (Paralog artifact)
    
    bcftools view \
        --threads "$THREADS" \
        -Oz \
        -e "QUAL <= 10 || MQ < 20 || INFO/DP < ${MIN_DP} || INFO/DP > ${MAX_DP} || MQBZ < -2.5 || BQBZ < -2.5 || RPBZ < -2.5 || RPBZ > 2.5 || SCBZ > 5 || COUNT(GT='het')=23" \
        -v snps \
        -M 2 \
        "$vcf_file" \
        -o "${OUT_DIR}/${base}_filtered.vcf.gz"

    # Index individual filtered files immediately
    bcftools index -t "${OUT_DIR}/${base}_filtered.vcf.gz"
done

# --- STEP 2: MERGING ---
echo "Starting merge..."

# Generate list of files to merge to avoid command line length limits
find "$OUT_DIR" -name "*_filtered.vcf.gz" > "${OUT_DIR}/merge_list.txt"

# Check if list is empty
if [[ ! -s "${OUT_DIR}/merge_list.txt" ]]; then
    echo "Error: No filtered VCFs found to merge."
    exit 1
fi

# Merge
bcftools merge \
    --threads "$THREADS" \
    -Oz \
    --file-list "${OUT_DIR}/merge_list.txt" \
    -o "${OUT_DIR}/merged_final.vcf.gz"

# Index final
bcftools index -t "${OUT_DIR}/merged_final.vcf.gz"

echo "======================================================"
echo "Pipeline Finished."
echo "Final Output: ${OUT_DIR}/merged_final.vcf.gz"
echo "======================================================"
