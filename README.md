# A Dual Genomic-Epigenomic Map of Clonal Evolution in Grapevine

**Author:** Paolo Callipo  
**Department of Plant Breeding, Hochschule Geisenheim University**

---

## 🍇 Abstract

This repository contains the custom code and analysis pipelines developed for the study **"A dual genomic-epigenomic map of clonal evolution in grapevine"**.

---

## 📂 Repository Structure

The analysis is partitioned into three main directories, following the logical flow of the manuscript. Each directory contains a detailed `README` explaining the specific scripts.

### 1. [Genome Assembly & Evaluation](./1_Genome_Assembly)
Scripts for characterizing the diploid reference genome structure.
*   **Structural Variation:** `SyRI` classification and `PlotSR` visualization.
*   **Composition:** Partitioning the genome into Hemizygous, Heterozygous, and Homozygous regions.
*   **Evolutionary History:** TE divergence/age estimation and Orthology analysis.

### 2. [Clonal Diversity Analysis](./2_Clonal_Diversity)
The pipeline for characterizing somatic variation across the 23 clones.
*   **Mapping:** Diploid-aware "Haplotype-masked" alignment strategy.
*   **Variant Calling:** Somatic SNP (`bcftools`) and SV (`Sniffles2`) identification.
*   **Phylogeny:** Reconstruction of clonal lineages using a Haploid Genomic Relationship Matrix (GRM) in R.

### 3. [Methylation Analysis](./3_Methylation_Analysis)
The custom pipeline to process the raw nanopore 5mC data and downstream analyses.
*   **Processing:** `modkit` pileup, context splitting, and standardization.
*   **Population Epigenetics:** Identification of Differentially Methylated Cytosines (DMCs).
*   **Gene Regulation:** Statistical classification of Gene Body Methylation (gbM) states and detection of epigenetic shifting.

---

## 🚀 Getting Started

### Prerequisites
The pipelines rely on standard bioinformatics tools and Python/R environments.

**Core Tools:**
*   `samtools`, `bcftools`, `bedtools`
*   `minimap2`
*   `modkit` (Nanopore methylation)
*   `Sniffles2` (Structural Variants)

**Python Dependencies:**
```bash
pip install pandas numpy scipy intervaltree pybedtools argparse
```
**R Dependencies:**
```R
install.packages(c("ape", "vegan", "pheatmap", "dendextend", "RColorBrewer"))
```

---

## 💾 Data Availability
*  **Raw Sequencing Data:** Available at the European Nucleotide Archive (ENA) under project PRJEB106155.
*  **Assembly & Annotations:** The phased diploid assembly (PN_1 and PN_2) and GFF3 files are available at https://zenodo.org/records/18154549.
*  **Modified PlotSR:** The custom visualization tool used for Figure 2 is included in 1_Assembly_Evaluation/plotsr_modified.

---

## 📞 Contact
If you have questions regarding the code or the analysis workflow, please open an Issue in this repository or contact: paolo.callipo@hs-gm.de
