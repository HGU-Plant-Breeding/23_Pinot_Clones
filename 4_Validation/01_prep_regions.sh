#!/bin/bash
# =============================================================
# 01_prep_regions.sh
# Generate all BED files needed for haplotype masking validation
# Callipo et al. 2026
# =============================================================

set -euo pipefail

WORK_DIR="/path"
OUT_DIR="/path"
FAI=".fai"
ROH_RAW="$.bed"

mkdir -p "${OUT_DIR}"
cd "${WORK_DIR}"

echo "=== Step 1: Preparing region BED files ==="

# -------------------------------------------------------
# 1a. Generate genome file for bedtools (masked naming)
#     FAI uses PN1_Hap_A → convert to PN1_HapA
# -------------------------------------------------------
echo "Generating genome size file..."
awk '{print $1"\t"$2}' "${FAI}" | \
    sed 's/Hap_A/HapA/g; s/Hap_B/HapB/g' \
    > "${OUT_DIR}/genome.sizes"
echo "  -> genome.sizes: $(wc -l < ${OUT_DIR}/genome.sizes) chromosomes"

# -------------------------------------------------------
# 1b. Clean RoH BED
#     - Already uses HapA/HapB naming — no rename needed
#     - Take only first 3 columns
#     - Remove any Chr0 if present
# -------------------------------------------------------
echo "Cleaning RoH BED..."
grep -v "Chr0" "${ROH_RAW}" | \
    cut -f1,2,3 | \
    sort -k1,1 -k2,2n \
    > "${OUT_DIR}/RoH.bed"
echo "  -> RoH.bed: $(wc -l < ${OUT_DIR}/RoH.bed) regions"

# -------------------------------------------------------
# 1c. Generate 5kb border regions flanking RoH
#     slop adds 5kb on each side, then subtract RoH itself
# -------------------------------------------------------
echo "Generating 5kb border regions..."
bedtools slop \
    -i "${OUT_DIR}/RoH.bed" \
    -g "${OUT_DIR}/genome.sizes" \
    -b 5000 | \
    bedtools subtract \
        -a - \
        -b "${OUT_DIR}/RoH.bed" | \
    sort -k1,1 -k2,2n | \
    bedtools merge \
    > "${OUT_DIR}/border_5kb.bed"
echo "  -> border_5kb.bed: $(wc -l < ${OUT_DIR}/border_5kb.bed) regions"

# -------------------------------------------------------
# 1d. Generate non-RoH, non-border regions
#     (rest of genome — used for background comparison)
# -------------------------------------------------------
echo "Generating non-RoH background regions..."
# First make full genome BED
awk '{print $1"\t0\t"$2}' "${OUT_DIR}/genome.sizes" \
    > "${OUT_DIR}/genome_full.bed"

# Subtract RoH and border
cat "${OUT_DIR}/RoH.bed" "${OUT_DIR}/border_5kb.bed" | \
    sort -k1,1 -k2,2n | \
    bedtools merge | \
    bedtools subtract \
        -a "${OUT_DIR}/genome_full.bed" \
        -b - | \
    sort -k1,1 -k2,2n \
    > "${OUT_DIR}/nonRoH.bed"
echo "  -> nonRoH.bed: $(wc -l < ${OUT_DIR}/nonRoH.bed) regions"

# -------------------------------------------------------
# 1e. Generate RoH BED in unmasked naming convention
#     For intersecting with unmasked BAMs/VCFs
#     HapA → RagTag, HapB → 2_RagTag
# -------------------------------------------------------
echo "Generating RoH BED in unmasked (RagTag) naming..."
sed 's/_HapA/_RagTag/g; s/_HapB/_2_RagTag/g' \
    "${OUT_DIR}/RoH.bed" \
    > "${OUT_DIR}/RoH_RagTag.bed"

sed 's/_HapA/_RagTag/g; s/_HapB/_2_RagTag/g' \
    "${OUT_DIR}/border_5kb.bed" \
    > "${OUT_DIR}/border_5kb_RagTag.bed"

sed 's/_HapA/_RagTag/g; s/_HapB/_2_RagTag/g' \
    "${OUT_DIR}/nonRoH.bed" \
    > "${OUT_DIR}/nonRoH_RagTag.bed"

echo "  -> RoH_RagTag.bed, border_5kb_RagTag.bed, nonRoH_RagTag.bed"

# -------------------------------------------------------
# Summary
# -------------------------------------------------------
echo ""
echo "=== Region summary ==="
echo "RoH total bp:     $(awk '{sum+=($3-$2)} END{printf "%d\n", sum}' ${OUT_DIR}/RoH.bed)"
echo "Border total bp:  $(awk '{sum+=($3-$2)} END{printf "%d\n", sum}' ${OUT_DIR}/border_5kb.bed)"
echo "NonRoH total bp:  $(awk '{sum+=($3-$2)} END{printf "%d\n", sum}' ${OUT_DIR}/nonRoH.bed)"
echo "Genome total bp:  $(awk '{sum+=($3-$2)} END{printf "%d\n", sum}' ${OUT_DIR}/genome_full.bed)"
echo ""
echo "Step 1 complete."
