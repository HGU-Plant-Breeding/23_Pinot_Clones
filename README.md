# A Dual Genome-Methylome Map of Clonal Evolution in Grapevine

**Authors:** Paolo Callipo, Hannah Robinson, Maximilian Schmidt, Kai P. Voss-Fels  
**Department of Plant Breeding, Hochschule Geisenheim University**

---

## Abstract

This repository contains all custom code and analysis pipelines developed for the study *"A dual genome-methylome map of clonal evolution in grapevine"*. The study generates a phased diploid reference genome for Pinot noir and integrates it with Oxford Nanopore sequencing of 23 distinct clones to build a genome-wide map of clonal genetic and DNA methylation variation in *Vitis vinifera*.

---

## Repository Structure

The repository is organised into five directories following the logical flow of the manuscript.

### [1. Genome Assembly](./1_Genome_Assembly)
Scripts for characterising the diploid reference genome structure, annotation, and epigenomic landscape.

### [2. Clonal Diversity](./2_Clonal_Diversity)
The full pipeline for mapping, somatic variant calling (SNPs and SVs), functional annotation, mechanistic classification, and clonal lineage reconstruction across the 23 clones.

### [3. Methylation Analysis](./3_Methylation_Analysis)
The complete methylation analysis pipeline, split into two parts:
- **MP_Calling/** — processing raw Nanopore methylation data into methylation polymorphisms (MPs)
- **GBM/** — gene body methylation classification, shifting gene detection, and GO enrichment

### [4. Validation](./4_Validation)
Validation analyses added during peer review: ONT vs PacBio HiFi methylome concordance, and quantitative assessment of the haplotype-masking strategy.

### [5. Genetic-Methylation Interplay](./5_Genetic_Methylation_Interplay)
Analyses of the relationship between genetic and epigenetic variation: MP-variant cis-associations and 5mC deamination contribution to somatic mutagenesis.

---

## Data Availability

- **Raw sequencing data:** European Nucleotide Archive (ENA), project accession [PRJEB106155](https://www.ebi.ac.uk/ena/browser/view/PRJEB106155)
- **Genome assembly, annotations, and methylation data:** [Zenodo record 18154549](https://zenodo.org/records/18154549)

---

## Dependencies

**Core bioinformatics tools:**
`samtools`, `bcftools`, `bedtools`, `minimap2`, `modkit`, `Sniffles2`, `SnpEff`, `BLAST+`, `deepTools`, `SyRI`, `Helixer`, `Liftoff`, `Mikado`, `EDTA`, `OrthoFinder`, `Tandem Repeats Finder`

**Python (v3.11+):**
```bash
pip install pandas numpy scipy intervaltree pybedtools matplotlib seaborn statsmodels
```

**R (v4.2+):**
```R
install.packages(c("ape", "vegan", "pheatmap", "dendextend", "RColorBrewer", "ggplot2"))
```

---

## Contact

For questions about the code or analysis workflow, please open an Issue or contact: paolo.callipo@hs-gm.de
