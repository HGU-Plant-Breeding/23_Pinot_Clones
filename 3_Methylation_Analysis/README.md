# 3. Methylation Analysis

The complete methylation analysis pipeline, from raw Nanopore data to methylation polymorphisms (MPs) and gene body methylation (GBM) classification. Organised into two subdirectories.

---

## MP_Calling/

Pipeline for identifying methylation polymorphisms (MPs) across the 23-clone panel using a binned, binarise-first approach. Run scripts in numbered order.

### Dependencies
**Tools:** `modkit`  
**Python 3.8+:** `pandas`, `numpy`, `intervaltree`

---

### `00_methylation_distributions.py`
Samples up to 2 million sites per clone (post coverage filter ≥4×) and plots genome-wide methylation density curves for CG, CHG, and CHH contexts. Used to derive the context-specific binarization thresholds empirically (Figure S18).

```bash
python 00_methylation_distributions.py \
    --cg-dir ./pileup/CG \
    --chg-dir ./pileup/CHG \
    --chh-dir ./pileup/CHH \
    --output methylation_distributions.pdf
```

---

### `01_parse_modkit.py`
Parses raw modkit bedMethyl files. For each context (CG, CHG, CHH), extracts Nmod and Nvalid counts. Merges symmetrical strand pairs for CG and CHG (1bp offset for CHG). CHH sites are written directly without merging.

```bash
python 01_parse_modkit.py input.bedMethyl output_prefix
# Produces: output_prefix.CG.bed, output_prefix.CHG.bed, output_prefix.CHH.bed
```

---

### `02_bin_methylation.py`
Aggregates per-site methylation counts from `01_parse_modkit.py` into fixed-size non-overlapping genomic bins. For each bin, sums Nmod and Nvalid across all cytosines and computes mean methylation percentage.

```bash
python 02_bin_methylation.py \
    sample.CG.bed \
    genome.chrom.sizes \
    sample.CG.bins.bed \
    --bin_size 200 \
    --min_sites 3
```

---

### `03_build_matrix.py`
Core MP calling script. Loads binned methylation files from all clones, applies coverage and missingness filters, binarizes each bin using context-specific thresholds (CG: <30%/≥70%; CHG: <25%/≥50%; CHH: <5%/≥15%), and identifies MPs as bins exhibiting both methylated and unmethylated states across the panel.

```bash
python 03_build_matrix.py \
    --binned_dir ./02_binned \
    --context CG \
    --output_prefix ./03_matrix/CG
```

---

### `04_calc_mp_maf.py`
Calculates the Minor Epiallele Frequency (MEF) for every MP in the binary matrix output of `03_build_matrix.py`. Used to generate the MEF site frequency spectra (Figure 5A).

```bash
python 04_calc_mp_maf.py CG.vmr.binary.tsv CG.maf.tsv
```

---

### `05_calc_mp_density.py`
Calculates MP rate (MP bins / callable bins) per genomic feature (Exon, Intron, Intergenic) for each chromosome. Uses majority overlap (>50%) for feature assignment with priority: Exon > Intron > Intergenic.

```bash
python 05_calc_mp_density.py \
    --continuous CG.continuous.tsv \
    --vmrs CG.vmr.binary.tsv \
    --exons exons.bed \
    --introns introns.bed \
    --intergenic intergenic.bed \
    --fai genome.fasta.fai \
    --output CG.mp_rate.tsv
```

---

## GBM/

Gene body methylation classification, shifting gene detection, allele-specific methylation analysis, and GO enrichment.

### Dependencies
**Python 3.8+:** `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `statsmodels`  
**External:** eggNOG-mapper annotations for GO enrichment; `go-basic.obo` (optional, for GO term labels)

---

### `genomic_cg_per_gene.py`
Counts CpG dinucleotides and CG% per gene from the genome FASTA and an exon BED file. Calculates genomic CG content independently of methylation calling, used to determine which genes have sufficient cytosine density for classification.

```bash
python genomic_cg_per_gene.py genome.fasta exons.bed --output cg_per_gene.tsv
```

### `GBM_classify_CG_CHG.py`
Classifies genes as **gbM** (gene-body methylated), **teM** (TE-like methylation), **UM** (unmethylated), or **Unclassified** using CG and CHG contexts only (CHH excluded due to ONT noise). Uses a one-sided binomial test against a per-clone genome-wide CDS background (pCG), with BH-FDR correction.

Classification logic:
- **teM**: CHG significantly above background AND CHG fraction ≥ threshold
- **gbM**: CG significant AND not teM AND sufficient CG/CHG coverage
- **UM**: CG not significant AND not teM AND CG fraction ≤ effect-size guard

```bash
python GBM_classify_CG_CHG.py \
    genes_cds.bed cg.bed chg.bed output_prefix \
    --min-n-cg 15 --min-n-chg 15 --format percent
```

### `GBM_classify_multiclone.py`
Driver script that runs `GBM_classify_CG_CHG.py` independently on each clone from a sample sheet, then merges per-clone results into a single wide table. Each clone uses its own per-clone background (pCG) to avoid bias from global methylation drift.

```bash
python GBM_classify_multiclone.py \
    samples.tsv genes_cds.bed ./gbm_output \
    --min-n-cg 15 --min-n-chg 15 --format percent --threads 4
```

Sample sheet format (tab-separated): `clone_name  cg_file  chg_file`

### `call_switching_genes.py`
Identifies genes that shift methylation state (gbM/teM/UM) across clones using conservative criteria: classified in ≥N clones, minority class supported by ≥2 clones, and effect-size guards on CG and CHG fraction differences between class groups.

```bash
python call_switching_genes.py ./gbm_output \
    --min-minority 2 --delta-cg 0.30 --delta-chg 0.10
```

### `go_enrichment.py`
GO enrichment analysis for shifting genes using Fisher's exact test with BH FDR correction. Reads eggNOG-mapper annotations and optionally a `go-basic.obo` file for GO term labels. Tests all shifting genes and key transition subsets (gbM↔teM, UM↔teM).

```bash
python go_enrichment.py \
    --emapper grapevine.emapper.annotations \
    --classification 20-13_classification.tsv \
    --switching switching_genes.tsv \
    --obo go-basic.obo \
    --outdir ./go_results
```

### `plot_cg_chg_scatter.py`
Generates a CG vs CHG methylation scatter plot with marginal histograms, coloured by classification (gbM/teM/UM/ambiguous). Used for QC and Figure S9.

```bash
python plot_cg_chg_scatter.py per_clone/20-13_classification.tsv \
    --output 20-13_scatter.png --min-n-cg 15 --min-n-chg 15
```

### `ortho_concordance.py`
Analyses methylation state concordance between the two haplotypes for 1-to-1 ortholog pairs. Produces a concordance heatmap and CG/CHG scatter panels (Figure S10).

```bash
python ortho_concordance.py
# Expects: one_to_one_orthologs.tsv, 20-13_classification.tsv
```

### `plot_AMP_bodies.py`
Visualises the gene-body methylation states of the 1,240 asymmetrically methylated promoter pairs (AMPs). Produces CG and CHG scatter plots and a categorical bar chart of gene-body outcomes (Figure S13).

```bash
python plot_AMP_bodies.py
# Expects: 20-13_classification.tsv, significant_AMPs.tsv
```
