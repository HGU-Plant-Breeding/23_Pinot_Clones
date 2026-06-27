# 4. Validation

Validation analyses added during peer review to support key methodological choices: concordance of the ONT-derived methylome with PacBio HiFi, and quantitative assessment of the haplotype-masking strategy for variant calling in runs of homozygosity (RoH).

---

## Dependencies

**Tools:** `samtools`, `bedtools`  
**Python 3.8+:** `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`

---

## ONT vs PacBio HiFi Methylome Validation

### `hifi_ont_comparison.py`
Loads per-site CG methylation from both platforms (ONT via modkit, PacBio HiFi via pb-CpG-tools), joins them on genomic coordinates, and computes Pearson and Spearman correlations at per-site and 200bp-binned resolution. Also performs a misclassification analysis at 30/70% binarization thresholds. Produces hexbin scatter and distribution figures.

Reported in the manuscript: per-site Pearson r = 0.959, 200bp-binned Pearson r = 0.982 (Figure S4).

```bash
python hifi_ont_comparison.py \
    --hifi PacBio_CG.bed \
    --ont ONT_CG.bed \
    --output hifi_ont_comparison \
    --min-cov 4 \
    --bin-size 200
```

---

## Haplotype Masking Validation

A numbered pipeline that quantifies how much the haplotype-masking strategy improves read mapping quality and variant recovery within RoH regions compared to standard unmasked alignment. Run scripts in order.

### `01_prep_regions.sh`
Generates all BED files needed for the analysis: RoH regions, 5kb border zones flanking RoH, and non-RoH background regions. Produces versions in both masked (HapA/HapB) and unmasked (RagTag) chromosome naming conventions.

```bash
bash 01_prep_regions.sh
```

### `02_coverage_mapq.slurm`
SLURM array job (one task per clone) that computes mean depth and mean MAPQ for masked and unmasked BAMs across the three region types (RoH, border, nonRoH). Also records per-read MAPQ distributions for histogram plots.

```bash
sbatch 02_coverage_mapq.slurm
```

### `03_aggregate_covmapq.py`
Aggregates per-clone coverage and MAPQ results from `02_coverage_mapq.slurm` into a single summary table.

```bash
python 03_aggregate_covmapq.py \
    --results_dir ./results/covmapq \
    --output covmapq_summary.tsv
```

### `04_plot_RoH_statistics.py`
Plots MAPQ distributions within RoH for masked vs unmasked strategies across all clones (Figures S14–S15).

```bash
python 04_plot_RoH_statistics.py \
    --summary covmapq_summary.tsv \
    --output roh_mapq_plots
```

### `05_plot.py`
Generates the final validation figures showing depth and MAPQ across all three region types and the improvement in variant recovery (Figures S16–S17).

```bash
python 05_plot.py \
    --summary covmapq_summary.tsv \
    --output masking_validation_figures
```
