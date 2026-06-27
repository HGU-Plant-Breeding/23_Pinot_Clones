# 5. Genetic-Methylation Interplay

Scripts investigating the relationship between somatic genetic variation and methylation polymorphisms: cis-associations between MPs and nearby SNPs/SVs, and the contribution of 5mC deamination to somatic mutagenesis.

---

## Dependencies

**Python 3.8+:** `pandas`, `numpy`, `scipy`, `matplotlib`

---

## Scripts

### `01_mp_genotype_association.py`
Computes phi coefficient associations between MP bins and nearby genetic variants (SNPs, SVs, or combined) within a defined cis window. For each MP bin, finds all variants within ±flank bp, computes phi coefficient and Fisher's exact p-value for each MP-variant pair, and reports the single best association per MP. Runs all three methylation contexts (CG, CHG, CHH) against all three variant types (SNP, SV, SNP+SV), producing nine output files.

Reported in the manuscript: at ±10kb, 63% of CG, 66% of CHG, and 62% of CHH MPs lack any nearby genetic variant; where associations exist, they are consistently weak (median |phi| = 0.13–0.15) (Figure S23).

```bash
python 01_mp_genotype_association.py \
    --vmr_dir ./mp_matrices \
    --snp GT_snp.tsv \
    --sv GT_SV_with_end.tsv \
    --output_dir ./associations \
    --flank 10000 \
    --min_samples 18
```

---

### `02_genome_overview_plot.py`
Generates a multi-page PDF with one page per chromosome pair (PN1–PN19). Each page shows SNP, SV, CG MP, CHG MP, and CHH MP density across both haplotypes in 50kb bins, with globally normalised y-axes per track type. Useful for visualising the genome-wide co-distribution of genetic and epigenetic variation.

```bash
python 02_genome_overview_plot.py \
    --snp GT_snp.tsv \
    --sv GT_SV_with_end.tsv \
    --cg_vmr CG.vmr.binary.tsv \
    --chg_vmr CHG.vmr.binary.tsv \
    --chh_vmr CHH.vmr.binary.tsv \
    --fai genome.fasta.fai \
    --output genome_overview.pdf \
    --bin_size 50000
```

---

### `03_5mc_deamination.py`
Tests whether somatic C>T transitions are enriched at highly methylated CG sites compared to sparsely methylated CG sites in the reference methylome. Computes C>T mutation rates across methylation deciles and performs a Fisher's exact test comparing rates at highly (>70%) vs sparsely (<30%) methylated CG sites. Also reports the full substitution type spectrum.

Reported in the manuscript: 6.3-fold enrichment of C>T transitions at highly methylated CG sites (0.457 vs 0.073 per 1,000 CG sites; Fisher's exact p < 0.0001) (Figure S22).

```bash
python 03_5mc_deamination.py \
    --vcf somatic_SNPs.vcf.gz \
    --methylome 20-13_CG.bed \
    --output 5mc_deamination.tsv
```
