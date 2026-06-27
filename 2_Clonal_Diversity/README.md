# 2. Clonal Diversity

The complete pipeline for characterising somatic variation across the 23 Pinot noir clones: read mapping, SNP and SV calling, functional annotation, mechanistic classification, and clonal lineage reconstruction.

---

## Dependencies

**Tools:** `samtools`, `minimap2`, `bcftools`, `Sniffles2`, `SnpEff`, `BLAST+`  
**Python 3.8+:** `pandas`, `intervaltree`  
**R:** `ape`, `vegan`, `pheatmap`, `dendextend`, `RColorBrewer`

---

## 1. Read Processing & Diploid Alignment

Reads were aligned to a **haplotype-masked diploid reference** to minimise reference bias and ensure accurate variant calling in runs of homozygosity (RoH). See Methods and `4_Validation/` for a quantitative assessment of this strategy.

### `convert_modbam_to_fastq.sh`
Converts unaligned Dorado BAM files to FASTQ, preserving the MM (Modified Base) and ML (Probability) tags required for downstream methylation analysis.

```bash
sbatch convert_modbam_to_fastq.sh input.bam output.fastq.gz
```

### `minimap_ont.sh`
Aligns ONT long reads to the masked diploid reference using Minimap2 (`map-ont` preset). Outputs sorted, indexed BAM files.

```bash
sbatch minimap_ont.sh SampleID ./fastq_dir reference_masked.fasta ./bam_output
```

---

## 2. Somatic SNP Calling

### `bcftools_calling.sh`
Joint variant calling across all 23 clones using `bcftools mpileup` with the ONT-specific error model (`-X ont`). Designed to run per-chromosome for parallelisation.

```bash
sbatch bcftools_calling.sh "Chr1" bam_list.txt reference.fasta ./vcf_out
```

### `filter_and_merge_vcfs.sh`
Applies the stringent hard-filtering pipeline to remove sequencing artefacts. Filters: QUAL > 10, mapping quality > 20, depth 115–690×, strand bias tests (MQBZ, BQBZ, RPBZ, SCBZ), and exclusion of loci heterozygous in all 23 samples (fixed germline variants). Merges per-chromosome VCFs into the final cohort VCF.

```bash
sbatch filter_and_merge_vcfs.sh ./raw_vcf_dir ./final_vcf_dir 115 690
```

### `snpeff_annotation.sh`
Annotates filtered somatic SNPs using SnpEff. Automatically builds a custom database from the diploid assembly if one does not exist. Classifies variants into impact categories (high, moderate, low, modifier).

```bash
sbatch snpeff_annotation.sh filtered.vcf.gz PinotNoir /path/to/snpeff ref.fasta annot.gff3
```

---

## 3. Structural Variant Calling

### `sniffles_single_sample.sh`
Step 1: Runs Sniffles2 on individual clones to generate `.snf` binary files. Uses `--minsupport 2` to maximise sensitivity for rare somatic events.

```bash
sbatch sniffles_single_sample.sh SampleID ./bam_dir reference.fasta ./snf_out
```

### `sniffles_merge_filter.sh`
Step 2: Merges population `.snf` files into a unified multi-sample VCF and applies filtering: PASS, PRECISE, QUAL > 30, STDEV_POS/LEN ≤ 50, size 50bp–50kb, excluding complex breakends (BNDs).

```bash
sbatch sniffles_merge_filter.sh snf_list.tsv Output_Prefix
```

---

## 4. SV Mechanistic Classification

Classification follows a hierarchical strategy. Run scripts in order:

### `parse_blast_TE_hits.py`
Parses BLASTN results of insertion sequences against the EDTA TE library. Assigns TE family to insertions with ≥80% identity over ≥80% query length.

```bash
# First run BLASTN:
# blastn -query insertions.fa -db TElib.fa -outfmt "6 qseqid sseqid pident length qlen"
python parse_blast_TE_hits.py blast_output.tsv --output te_insertion_families.tsv
```

### `classify_TE_types.py`
First-pass SV classification into: TE-derived (from BLAST results), Centromeric/Satellite (length matches multiples of grapevine satellite repeat units: 79, 107, 135, 187 bp ±1bp), or Unknown.

```bash
python classify_TE_types.py sv_list.tsv \
    --te-ins te_insertion_families.tsv \
    --te-del te_deletion_families.tsv \
    --output sv_classification.tsv
```

### `refine_unknown_svs.py`
Second-pass classification of SVs in the Unknown category. Decomposes them into: Tandem repeats (TRF, non-centromeric), Segmental duplications (BLASTN against diploid reference, ≥90% identity, ≥50% coverage), Organellar insertions/NUMTs/NUPTs (BLASTN against mitochondrial/chloroplast genomes), and Complex/Unclassified.

```bash
python refine_unknown_svs.py \
    --svs all_svs_with_lengths.tsv \
    --insertions insertions.fasta \
    --deletions deletions.bed \
    --te-ins te_insertion_families.tsv \
    --te-del te_deletion_families.tsv \
    --reference genome.fasta \
    --edta-lib TElib.fa
```

### `calc_variant_density.py`
Calculates SNP and SV density (events per 100kb) across exons, introns, and intergenic regions for all 38 pseudo-chromosomes. Uses interval trees to correctly handle SVs spanning multiple features.

```bash
python calc_variant_density.py \
    --variants variants.bed \
    --exons exons.bed \
    --introns introns.bed \
    --intergenic intergenic.bed \
    --fai genome.fasta.fai \
    --output variant_density.tsv
```

---

## 5. Clonal Lineage Reconstruction

### `clonal_lineage_analysis.R`
Reconstructs clonal phylogenies from all molecular layers (SNPs, SVs, CG MPs, CHG MPs, CHH MPs) using a Haploid Genomic Relationship Matrix (GRM) adapted for somatic/binary markers. Produces UPGMA dendrograms, tanglegrams comparing SNP vs CG and SNP vs CHG topologies (Figure 6), pairwise distance heatmaps, and Mantel tests for concordance between all layer pairs.

```bash
Rscript clonal_lineage_analysis.R
```

> Input files are expected in a `Distance_Matrix/` subdirectory. See script header for the full file list.
