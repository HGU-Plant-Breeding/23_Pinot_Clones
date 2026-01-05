#!/bin/bash
#==============================================================================
# SCRIPT: run_snpeff_annotation.sh
# DESCRIPTION: Annotates a VCF file using SnpEff. 
#              - Automatically builds the SnpEff database if it doesn't exist.
#              - Requires Reference FASTA and GFF3 for the build step.
#
# USAGE: sbatch run_snpeff_annotation.sh <VCF_IN> <GENOME_NAME> <SNPEFF_DIR> [REF_FASTA] [REF_GFF]
# EXAMPLE: sbatch run_snpeff_annotation.sh variants.vcf.gz PinotNoir /path/to/snpeff ref.fa annot.gff
#==============================================================================

#SBATCH --job-name=snpeff_ann
#SBATCH --output=logs/snpeff_%j.out
#SBATCH --error=logs/snpeff_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -e

# --- MODULES ---
# Adjust based on your HPC. If using a local jar, java is sufficient.
module load java

# --- ARGUMENTS ---
VCF_IN=${1}
GENOME_NAME=${2}  # e.g. "PinotNoir"
SNPEFF_DIR=${3}   # Directory containing snpEff.jar and snpEff.config
REF_FASTA=${4:-""} # Optional if DB exists
REF_GFF=${5:-""}   # Optional if DB exists

# --- VALIDATION ---
if [[ -z "$VCF_IN" || -z "$GENOME_NAME" || -z "$SNPEFF_DIR" ]]; then
    echo "Usage: sbatch $0 <VCF_IN> <GENOME_NAME> <SNPEFF_DIR> [REF_FASTA] [REF_GFF]"
    exit 1
fi

# Define paths
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
SNPEFF_CONFIG="${SNPEFF_DIR}/snpEff.config"
DB_DIR="${SNPEFF_DIR}/data/${GENOME_NAME}"
OUTPUT_VCF="${VCF_IN%.vcf.gz}_annotated.vcf.gz"
STATS_HTML="${VCF_IN%.vcf.gz}_snpeff_stats.html"
STATS_CSV="${VCF_IN%.vcf.gz}_snpeff_stats.csv"

# --- STEP 1: BUILD DATABASE (If Missing) ---

if [ -d "$DB_DIR" ]; then
    echo "[Info] SnpEff database '${GENOME_NAME}' found. Skipping build."
else
    echo "[Info] SnpEff database '${GENOME_NAME}' not found. Attempting to build..."
    
    if [[ -z "$REF_FASTA" || -z "$REF_GFF" ]]; then
        echo "Error: Database missing. You must provide REF_FASTA and REF_GFF arguments to build it."
        exit 1
    fi

    # 1. Create directory structure
    mkdir -p "$DB_DIR"

    # 2. Copy reference files (SnpEff requires specific names)
    echo "Copying reference files..."
    cp "$REF_FASTA" "${DB_DIR}/sequences.fa"
    cp "$REF_GFF" "${DB_DIR}/genes.gff"

    # 3. Update config file (Check if entry exists first)
    if ! grep -q "${GENOME_NAME}.genome" "$SNPEFF_CONFIG"; then
        echo "Updating snpEff.config..."
        echo "" >> "$SNPEFF_CONFIG"
        echo "# Custom Genome: ${GENOME_NAME} (Added by pipeline)" >> "$SNPEFF_CONFIG"
        echo "${GENOME_NAME}.genome : ${GENOME_NAME}" >> "$SNPEFF_CONFIG"
    fi

    # 4. Build
    echo "Building database..."
    java -Xmx8g -jar "$SNPEFF_JAR" build -gff3 -v "$GENOME_NAME"
    echo "Database build complete."
fi

# --- STEP 2: ANNOTATION ---

echo "Starting Annotation..."
echo "Input: $VCF_IN"
echo "Genome: $GENOME_NAME"

java -Xmx16g -jar "$SNPEFF_JAR" ann \
    -v \
    -stats "$STATS_HTML" \
    -csvStats "$STATS_CSV" \
    "$GENOME_NAME" \
    "$VCF_IN" \
    | gzip > "$OUTPUT_VCF"

echo "========================================================"
echo "SnpEff Pipeline Complete."
echo "Annotated VCF: $OUTPUT_VCF"
echo "Stats: $STATS_HTML"
echo "========================================================"
