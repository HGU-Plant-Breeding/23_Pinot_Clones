# Clonal Diversity and Somatic Evolution Analysis

This directory contains the computational pipeline used to characterize somatic variation across the 23 Pinot noir clones. 

The analysis covers **Read Mapping**, **Variant Calling (SNPs & SVs)**, **Functional Annotation**, **Mechanistic Classification**, and **Phylogenetic Reconstruction**.

## Dependencies

*   **Bioinformatics Tools:** `samtools`, `minimap2`, `bcftools`, `sniffles`, `SnpEff`, `BLAST+`
*   **Python:** `pandas`, `intervaltree`, `argparse`
*   **R:** `ape`, `vegan`, `pheatmap`, `dendextend`

---

## 1. Read Processing and Diploid Alignment

To minimize reference bias and correctly handle Runs of Homozygosity (RoH), reads were aligned to a **haplotype-masked diploid reference** (see Manuscript Methods).

### `convert_modbam_to_fastq.sh`
**Description:**  
Converts unaligned BAM files (output by Dorado basecalling) to FASTQ format. Crucially, it preserves the **MM (Modified Base)** and **ML (Probability)** tags, ensuring that methylation information is retained for downstream epigenetic analysis.
**Usage:**
```bash
sbatch convert_modbam_to_fastq.sh input.bam output.fastq
```
### `minimap_ont.sh`
**Description:**
Aligns long reads to the reference genome using minimap2. It sorts and indexes the output to generate analysis-ready BAM files.
**Usage:**

```Bash
sbatch minimap_ont.sh SampleID ./fastq_dir reference_masked.fasta ./bam_output
```
## 2. Somatic SNP Calling
### bcftools_calling.sh
**Description:**
Performs joint variant calling using bcftools mpileup.
Optimization: Uses the -X ont parameter to apply Nanopore-specific error models.
Parallelization: Designed to run on specific genomic regions (e.g., per chromosome) to expedite processing.
**Usage:**

```Bash
sbatch bcftools_calling.sh "Chr1" bam_list.txt reference.fasta ./vcf_out
```
### filter_and_merge_vcfs.sh
**Description:**
Applies the stringent Hard Filtering pipeline described in the manuscript to remove sequencing artifacts and paralogous mappings.
Filters: Quality (QUAL>10), Depth (115<DP<690), Strand Bias, and Fixed Heterozygosity (COUNT(het)=23).
Merging: Combines regional VCFs into the final cohort VCF.
**Usage:**
```bash
sbatch filter_and_merge_vcfs.sh ./raw_vcf_dir ./final_vcf_dir 115 690
```
### snpeff_annotation.sh
**Description:**
Annotates the filtered somatic SNPs using SnpEff to predict functional impacts (e.g., missense, nonsense, intergenic). It automatically builds a custom SnpEff database for the Pinot noir diploid assembly if one does not exist.
**Usage:**
```Bash
sbatch snpeff_annotation.sh filtered.vcf.gz PinotNoir /path/to/snpeff
```
## 3. Structural Variant (SV) Analysis
### sniffles_single_sample.sh
**Description:**
Step 1 of SV discovery. Runs Sniffles2 on individual clones to generate binary .snf files. Uses --minsupport 2 to maximize sensitivity for rare somatic events.
**Usage:**
```Bash
sbatch sniffles_single_sample.sh SampleID ./bam_dir reference.fasta ./snf_out
```
### sniffles_merge_filter.sh
**Description:**
Step 2 of SV discovery. Merges population .snf files into a unified multi-sample VCF.
Refinement: Automatically filters the output using bcftools to remove imprecise calls, extremely large events (>50kb), and complex breakends (BNDs), retaining only high-confidence insertions and deletions.
**Usage:**
```Bash
sbatch sniffles_merge_filter.sh snf_list.tsv Output_Prefix
```
## 4. Variant Characterization & Mechanism
### parse_blast_TE_hits.py
**Description:**
Determines if somatic Insertions are Transposable Elements (TEs). It parses BLAST results of insertion sequences against the EDTA TE library and assigns identity based on the best bit-score hit.
**Usage:**
```Bash
# First run BLAST: blastn -query insertions.fa -db TE_lib -outfmt "6 ..."
python parse_blast_TE_hits.py blast_output.tsv --output te_families.tsv
```
### classify_TE_types.py
**Description:**
Classifies Structural Variants into mechanistic categories:
TE-derived: (Matched via parse_blast_TE_hits.py)
Centromeric/Satellite: (Length matches multiples of 107/79/135/187 bp)
Unknown/Other
**Usage:**
```Bash
python classify_TE_types.py sv_list.tsv --te-ins te_families.tsv ...
```
### calc_variant_density.py
**Description:**
Calculates the genomic density of somatic variants (SNPs or SVs) across Exons, Introns, and Intergenic regions. It uses interval trees to correctly handle SVs that span multiple features, assigning priority to coding regions.
**Usage:**
```Bash
python calc_variant_density.py --variants variants.bed --exons exons.bed ...
```
## 5. Clonal Lineage Reconstruction
### clonal_lineage_analysis.R
**Description:**
The unified R script used to generate Figure 6 and Figure S11.
Methodology: Calculates a Haploid Genomic Relationship Matrix (GRM) tailored for somatic variation.
Phylogeny: Reconstructs clonal trees using UPGMA.
Comparison: Performs a Mantel test and generates Tanglegrams to visualize the concordance between the Genetic Tree (SNPs) and the Epigenetic Tree (CG-Methylation).
**Usage:**
```Bash
Rscript clonal_lineage_analysis.R
```
