#!/bin/bash
#SBATCH --job-name=haplo_mask_covmapq
#SBATCH --output=/
#SBATCH --error=/
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-22

# =============================================================
# 02_coverage_mapq.slurm — v3
# Uses samtools view + samtools depth directly
# No piping to samtools coverage (doesn't accept stdin)
# No temp files needed
# =============================================================

WORK_DIR="/dir"
REG_DIR="${WORK_DIR}/results/regions"
OUT_DIR="${WORK_DIR}/results/covmapq"
BAM_MASKED="${WORK_DIR}/Mapped_BAM_Masked"
BAM_UNMASKED="${WORK_DIR}/Mapped_BAM_Umasked"

mkdir -p "${OUT_DIR}"

CLONES=(
    "1-32" "1-47" "20-13-28-8" "20-16" "20-19-19" "20-20"
    "20-24-3" "20-25-2" "20-27-3" "2-2-13" "2-6-41" "2-9-35"
    "2-9-40" "37-1" "4000-1" "4000-3" "4047-1" "404"
    "Kastenholz" "Probstei" "RC" "S1-1" "S1-2"
)

CLONE="${CLONES[$SLURM_ARRAY_TASK_ID]}"
echo "Processing clone: ${CLONE}"

# -------------------------------------------------------
# Function: mean depth via samtools depth
# samtools depth -b BED outputs: chrom pos depth
# We average across all positions
# -------------------------------------------------------
get_mean_depth() {
    local BAM=$1
    local BED=$2

    samtools depth \
        -@ "${SLURM_CPUS_PER_TASK}" \
        -b "${BED}" \
        -a \
        "${BAM}" | \
    awk '{sum+=$3; n++} END{if(n>0) print sum/n; else print 0}'
}

# -------------------------------------------------------
# Function: mean MAPQ and MAPQ=0 fraction via samtools view
# Excludes supplementary reads (flag 2048)
# -------------------------------------------------------
get_mean_mapq() {
    local BAM=$1
    local BED=$2

    samtools view \
        -@ "${SLURM_CPUS_PER_TASK}" \
        -F 2048 \
        -L "${BED}" \
        "${BAM}" | \
    awk '{sum+=$5; if($5==0) z++; n++}
         END{if(n>0) printf "%.4f\t%.4f\n", sum/n, z/n*100;
             else printf "0\t0\n"}'
}

# -------------------------------------------------------
# Function: run all stats for one BAM x one region
# Output: clone strategy region mean_depth mean_mapq pct_mapq0
# -------------------------------------------------------
run_stats() {
    local BAM=$1
    local BED=$2
    local REGION_TYPE=$3
    local STRATEGY=$4
    local CLONE=$5

    DEPTH=$(get_mean_depth "${BAM}" "${BED}")
    MAPQ_STATS=$(get_mean_mapq "${BAM}" "${BED}")
    MEAN_MAPQ=$(echo "${MAPQ_STATS}" | cut -f1)
    PCT_MAPQ0=$(echo "${MAPQ_STATS}" | cut -f2)

    echo -e "${CLONE}\t${STRATEGY}\t${REGION_TYPE}\t${DEPTH}\t${MEAN_MAPQ}\t${PCT_MAPQ0}"
}

# -------------------------------------------------------
# Function: per-read MAPQ distribution (for histogram)
# -------------------------------------------------------
run_mapq_dist() {
    local BAM=$1
    local BED=$2
    local STRATEGY=$3
    local CLONE=$4
    local OUTFILE=$5

    samtools view \
        -@ "${SLURM_CPUS_PER_TASK}" \
        -F 2048 \
        -L "${BED}" \
        "${BAM}" | \
    awk -v clone="${CLONE}" \
        -v strategy="${STRATEGY}" \
        'BEGIN{OFS="\t"} {print clone, strategy, $5}' \
    >> "${OUTFILE}"
}

# -------------------------------------------------------
# MASKED BAM
# -------------------------------------------------------
echo "  Processing masked BAM..."
BAM_M="${BAM_MASKED}/${CLONE}.bam"

if [ ! -f "${BAM_M}" ]; then
    echo "ERROR: ${BAM_M} not found"; exit 1
fi

{
run_stats "${BAM_M}" "${REG_DIR}/RoH.bed"        "RoH"    "masked" "${CLONE}"
run_stats "${BAM_M}" "${REG_DIR}/border_5kb.bed"  "border" "masked" "${CLONE}"
run_stats "${BAM_M}" "${REG_DIR}/nonRoH.bed"      "nonRoH" "masked" "${CLONE}"
} > "${OUT_DIR}/${CLONE}_masked_coverage.tsv"

echo "  Running masked MAPQ distribution in RoH..."
run_mapq_dist \
    "${BAM_M}" "${REG_DIR}/RoH.bed" \
    "masked" "${CLONE}" \
    "${OUT_DIR}/${CLONE}_masked_mapq_RoH.tsv"

echo "  Running masked MAPQ distribution in nonRoH..."
run_mapq_dist \
    "${BAM_M}" "${REG_DIR}/nonRoH.bed" \
    "masked" "${CLONE}" \
    "${OUT_DIR}/${CLONE}_masked_mapq_nonRoH.tsv"

# -------------------------------------------------------
# UNMASKED BAM
# -------------------------------------------------------
echo "  Processing unmasked BAM..."
BAM_U="${BAM_UNMASKED}/${CLONE}.bam"

if [ ! -f "${BAM_U}" ]; then
    echo "ERROR: ${BAM_U} not found"; exit 1
fi

{
run_stats "${BAM_U}" "${REG_DIR}/RoH_RagTag.bed"        "RoH"    "unmasked" "${CLONE}"
run_stats "${BAM_U}" "${REG_DIR}/border_5kb_RagTag.bed"  "border" "unmasked" "${CLONE}"
run_stats "${BAM_U}" "${REG_DIR}/nonRoH_RagTag.bed"      "nonRoH" "unmasked" "${CLONE}"
} > "${OUT_DIR}/${CLONE}_unmasked_coverage.tsv"

echo "  Running unmasked MAPQ distribution in RoH..."
run_mapq_dist \
    "${BAM_U}" "${REG_DIR}/RoH_RagTag.bed" \
    "unmasked" "${CLONE}" \
    "${OUT_DIR}/${CLONE}_unmasked_mapq_RoH.tsv"

echo "  Running unmasked MAPQ distribution in nonRoH..."
run_mapq_dist \
    "${BAM_U}" "${REG_DIR}/nonRoH_RagTag.bed" \
    "unmasked" "${CLONE}" \
    "${OUT_DIR}/${CLONE}_unmasked_mapq_nonRoH.tsv"

echo "  Clone ${CLONE} done."
