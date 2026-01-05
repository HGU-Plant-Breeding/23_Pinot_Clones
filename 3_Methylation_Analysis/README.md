# Methylation Analysis Pipeline

This directory contains the computational workflow used to characterize the methylation landscape of the 23 Pinot Noir clones. 
The pipeline processes Nanopore methylation data to identify **Differentially Methylated Cytosines (DMCs)**, analyze **Gene Body Methylation (gbM)** dynamics, and generate the data matrices required for clonal lineage reconstruction.

## Dependencies

*   **Bioinformatics Tools:**
    *   `modkit` (Nanopore modification calling)
    *   `samtools`
*   **Python Libraries:**
    *   `pandas`, `numpy`, `scipy`
    *   `intervaltree`, `argparse`

---

## 1. Methylation Calling & Extraction

### `run_modkit_pileup.sh`
**Description:**  
A SLURM submission script that runs `modkit` on aligned BAM files (containing MM/ML tags) to extract 5mC probabilities and prepare the bedMethyl files for each context for downstream analysis.
*   **Output:** Generates a genome-wide pileup bedMethyl file and automatically splits it into context-specific files (`CG`, `CHG`, `CHH`) for downstream processing.
*   **Reference:** Requires the diploid reference FASTA.

**Usage:**
```bash
sbatch run_modkit_pileup.sh SampleID ./bam_dir reference.fasta ./meth_out
```
## 2. Data Standardization & Binarization
### `process_methylation_bed.py`
**Description:**
Parses raw modkit BED files to standardize calls.
CG/CHG Contexts: Merges symmetrical sites from the forward (+) and reverse (-) strands into a single high-confidence call per site using weighted averages. Unpaired sites are discarded to ensure robustness.
CHH Context: Extracts relevant data without merging (as CHH is asymmetric).
**Usage:**
```Bash
python process_methylation_bed.py raw_CG.bed merged_CG.bed --context CG
```
### `binarize_methylation.py`
**Description:**
Converts continuous methylation percentages into discrete binary states to create an "Epigenetic Genotype" for phylogenetic analysis.
0 (Unmethylated): Methylation level ≤30%
1 (Methylated): Methylation level ≥70%
. (Missing): Sites with low coverage (< 4x) or intermediate methylation (30-70%) because considered ambigous.
**Usage:**
```Bash
python binarize_methylation.py merged_CG.bed binary_CG.tsv --min-cov 4 --low 30 --high 70
```
## 3. Population Matrix Generation
###  `merge_methylation_matrix.py`
**Description:**
Aggregates the individual binarized files from all 23 clones into a single Population Matrix (Rows = CpG Sites, Columns = Clones).
Method: Performs an Outer Join on site IDs. Missing data (where a site is covered in some clones but not others) is filled with . values.
**Usage:**
```Bash
python merge_methylation_matrix.py population_matrix_CG.tsv *_CG_binary.tsv
```
### `filter_methylation_matrix.py`
**Description:**
Filters the raw population matrix to retain only informative Differentially Methylated Cytosines (DMCs).
Missingness Filter: Removes sites with >10% missing data across the population.
MAF Filter: Removes invariant sites and rare variants (Minor Allele Frequency < 0.05).
Output: This filtered matrix is the direct input for the Clonal Lineage Reconstruction (Folder 2).
**Usage:**
```Bash
python filter_methylation_matrix.py population_matrix_CG.tsv filtered_DMCs.tsv --missing 0.1 --maf 0.05
```
## 4. Genomic Distribution & Statistics
### `calc_maf_spectrum.py`
**Description:**
Calculates the Minor Allele Frequency (MAF) for every site in the matrix. The output is used to generate Site Frequency Spectrum (SFS) plots (Figure 5A) to compare the evolutionary stability of CG vs. non-CG methylation.
**Usage:**
```Bash
python calc_maf_spectrum.py filtered_DMCs.tsv maf_results.tsv
```
### `calc_dmc_density.py`
**Description:**
Calculates the genomic density of DMCs (events per 100kb) across Exons, Introns, and Intergenic regions. It uses interval trees to map DMCs to features, prioritizing coding regions (Exon > Intron > Intergenic).
**Usage:**

```Bash
python calc_dmc_density.py \
    --dmcs filtered_DMCs.tsv \
    --exons exons.bed --introns introns.bed --intergenic intergenic.bed \
    --fai genome.fasta.fai \
    --output dmc_density_stats.tsv
```
## 5. Gene Body Methylation (gbM) Analysis
###  `GBM_calc.py`
**Description:**
Performs the statistical classification of Gene Body Methylation (gbM) states for the entire cohort.
Calculates a clone-specific global methylation background (_pCG_).
Performs a Binomial Test for every gene to determine if it is significantly methylated (BM) or undermethylated (UM) compared to the background.
Generates summary matrices (Class Matrix and P-value Matrix) for all clones.
**Usage:**
```Bash
python GBM_calc.py \
    genes.bed \
    --clone "Clone_20-13=path/to/20-13_CG.bed" \
    --clone "Clone_Ab48=path/to/Ab48_CG.bed" \
    --outdir gbm_results \
    --format percent
```
###  `analyze_gbm_stability.py`
**Description:**
Analyzes the stability of epigenetic states across the population using the output from GBM_calc.py.
gbM Shift: Identifies genes that switch states (e.g., BM in one clone, UM in another).
Stability Metrics: Calculates the fraction of clones holding a specific state per gene.
**Usage:**
```Bash
python analyze_gbm_stability.py gbm_results/gbm_class_matrix.tsv --output gbm_stability_report.tsv
```
