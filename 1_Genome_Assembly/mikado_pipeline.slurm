#!/bin/bash
#=========================================================================================
# SCRIPT NAME: run_mikado_pipeline.sh
# DESCRIPTION: SLURM submission script for the Mikado Annotation Pipeline.
#              Consolidates transcript assemblies, aligns proteins (BLAST), 
#              predicts ORFs (TransDecoder), and picks the best loci.
#
# USAGE: sbatch run_mikado_pipeline.sh
#=========================================================================================

#SBATCH --job-name=Mikado_Annotation
#SBATCH --partition=standard         # ACTION: Change to your HPC partition
#SBATCH --cpus-per-task=32           # Adjust based on availability
#SBATCH --mem=160G                   # Adjust based on genome size
#SBATCH --output=mikado.%j.out
#SBATCH --error=mikado.%j.err

# Exit on error
set -e

# --- USER CONFIGURATION ---
# Set the directory where your input files are located
PROJECT_DIR="." 

# Input Filenames
REFERENCE=""
PROTEIN_DB=""
INPUT_LIST="list.tsv"                 # Tab-separated list of transcript assemblies
CONFIG_FILE="configuration.yaml"
SCORING_FILE="plant.yaml"

# Load modules (Adjust names based on your HPC environment)
module load mikado
module load blast
module load seqkit
module load transdecoder

# Move to project directory
cd "$PROJECT_DIR"

echo "======================================================"
echo "Starting Mikado Job: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo "======================================================"

# --- Step 1: Configure ---
# Generates the configuration file for Mikado
echo "[Step 1] Configuring Mikado..."
if [ ! -f "$CONFIG_FILE" ]; then
    mikado configure \
      --list "$INPUT_LIST" \
      --reference "$REFERENCE" \
      --mode permissive \
      --scoring "$SCORING_FILE" \
      --copy-scoring "$SCORING_FILE" \
      --blast_targets "$PROTEIN_DB" \
      "$CONFIG_FILE"
else
    echo "Configuration file $CONFIG_FILE already exists. Skipping generation."
fi

# --- Step 2: Prepare ---
# Cleans and prepares the transcript assemblies
echo "[Step 2] Preparing assemblies..."
mikado prepare --json-conf "$CONFIG_FILE"

# --- Step 3: Homology Search (BLAST) ---
# Aligns the prepared transcripts against the protein database
echo "[Step 3] Running BLAST..."

# Create BLAST DB if it doesn't exist
if [ ! -f "${PROTEIN_DB}.pin" ]; then
    echo "Creating BLAST database..."
    makeblastdb -in "$PROTEIN_DB" -dbtype prot -parse_seqids
fi

# Run BLASTX
# Note: This step can be time-consuming.
if [ ! -f "mikado_prepared.blast.tsv" ]; then
    blastx -query mikado_prepared.fasta -db "$PROTEIN_DB" \
      -evalue 1e-10 -num_threads "$SLURM_CPUS_PER_TASK" -max_target_seqs 5 \
      -out mikado_prepared.blast.tsv \
      -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop"
else
    echo "BLAST output found. Skipping alignment."
fi

# --- Step 4: ORF Prediction (TransDecoder) ---
# Predicts open reading frames in the transcripts
echo "[Step 4] Predicting ORFs with TransDecoder..."
TransDecoder.LongOrfs -t mikado_prepared.fasta
TransDecoder.Predict -t mikado_prepared.fasta --single_best_only

# --- Step 5: Serialise ---
# Consolidates all data (structure, blast, orfs) into a database
echo "[Step 5] Serialising data..."
mikado serialise --json-conf "$CONFIG_FILE" \
  --orfs mikado_prepared.fasta.transdecoder.gff3 \
  --xml  mikado_prepared.blast.tsv

# --- Step 6: Pick ---
# Selects the best gene models
echo "[Step 6] Picking best loci..."
mikado pick --configuration "$CONFIG_FILE" --subloci-out mikado.subloci.gff3

echo "======================================================"
echo "Mikado Pipeline Finished Successfully"
echo "End Time: $(date)"
echo "======================================================"
