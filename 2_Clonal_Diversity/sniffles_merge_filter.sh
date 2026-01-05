#!/bin/bash
#==============================================================================
# SCRIPT: run_sniffles_merge.sh
# DESCRIPTION: Step 2 of SV Calling. 
#              1. Combines multiple .snf files into a raw multi-sample VCF.
#              2. Applies post-processing filters using bcftools.
#
# FILTERS APPLIED:
#   - Must be PASS
#   - Must be PRECISE (PRECISE=1)
#   - Quality >= 30
#   - Position/Length Deviation <= 50
#   - Length <= 50,000 bp
#   - Exclude Breakends (BND)
#
# USAGE: sbatch run_sniffles_merge.sh <SNF_LIST_FILE> <OUTPUT_PREFIX>
# EXAMPLE: sbatch run_sniffles_merge.sh snf_list.tsv my_cohort
#==============================================================================

#SBATCH --job-name=Sniffles_Merge
#SBATCH --partition=all
#SBATCH --cpus-per-task=16
#SBATCH --threads-per-core=2
#SBATCH --time=2:00:00
#SBATCH --output=logs/sniff_merge_%j.out
#SBATCH --error=logs/sniff_merge_%j.err

set -e

# --- MODULES ---
module load sniffles
module load bcftools

# --- ARGUMENTS ---
SNF_LIST=${1}
OUT_PREFIX=${2} # e.g., "Pinot_Clones"

if [[ -z "$SNF_LIST" || -z "$OUT_PREFIX" ]]; then
    echo "Usage: sbatch $0 <SNF_LIST_FILE> <OUTPUT_PREFIX>"
    exit 1
fi

THREADS=${SLURM_CPUS_PER_TASK:-16}
RAW_VCF="${OUT_PREFIX}_raw.vcf"
FINAL_VCF="${OUT_PREFIX}_filtered.vcf.gz"

echo "======================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Merging files listed in: $SNF_LIST"
echo "Output Prefix: $OUT_PREFIX"
echo "======================================================"

# --- STEP 1: Sniffles Merge ---
echo "[Step 1] Running Sniffles Merge..."

sniffles \
    --input "$SNF_LIST" \
    --vcf "$RAW_VCF" \
    --threads "$THREADS" \
    --combine-low-confidence-abs 1

# --- STEP 2: Filtering ---
echo "[Step 2] Applying Filters..."

# Filter Logic:
# 1. -f PASS: Keep only PASS variants
# 2. -e '...': Exclude variants matching these criteria (Imprecise, Low Qual, High Deviation, Huge Size)
# 3. -e 'SVTYPE="BND"': Exclude Breakends (translocations/complex)

bcftools view \
    -f PASS \
    -e 'PRECISE=0 || QUAL < 30 || STDEV_POS > 50 || STDEV_LEN > 50 || ABS(SVLEN) > 50000' \
    "$RAW_VCF" \
    | bcftools view -e 'SVTYPE="BND"' -Oz -o "$FINAL_VCF"

# Index the final file
bcftools index "$FINAL_VCF"

echo "======================================================"
echo "Pipeline Complete."
echo "Raw VCF:   $RAW_VCF"
echo "Final VCF: $FINAL_VCF"
echo "======================================================"
