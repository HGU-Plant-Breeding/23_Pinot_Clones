# Enhanced PlotSR with Variant Density Heatmaps

This repository contains a modified version of **[PlotSR](https://github.com/schneebergerlab/plotsr)**, extended to transform chromosome visualizations from schematic representations into data-rich quantitative plots. 

The core enhancement replaces solid-colored chromosome bars with **variant density heatmaps**. This allows for the simultaneous visualization of large-scale structural variants (SVs) via ribbons and the fine-scale landscape of nucleotide diversity (SNPs/indels) via the chromosome tracks.

## Key Features

### 1. Variant Density Heatmaps
Instead of abstract colored bars, chromosomes are rendered as heatmaps representing the density of small variants (SNPs, indels) in user-defined genomic windows.
*   **Logic:** The script calculates the raw count of variants within genomic windows (e.g., 100kb) and normalizes them against a global maximum (default: 98th percentile) to prevent hotspots from washing out the signal.
*   **Visuals:** High-density regions appear "hotter" (darker/more intense color), while conserved regions appear lighter.

### 2. Aesthetic Refinements
*   **Controllable Bar Thickness:** Users can adjust the thickness of the heatmap bars to suit the figure's aspect ratio.
*   **Smart Legends:** Redundant "Genomes" legends are automatically hidden when SNP density is active to reduce clutter.
*   **Integrated Color Bar:** A horizontal color bar indicating variant density is automatically placed in the top-left corner.

---

## Installation & Dependencies

This version relies on the standard PlotSR dependencies. Ensure you have the following installed:

*   Python >= 3.8
*   pandas
*   numpy
*   matplotlib

---

## Usage Guide

### Step 1: Format Variant Data
PlotSR requires a specific simplified format for variant data (`GenomeID`, `Chromosome`, `Position`). Use the provided helper script to convert your data (e.g., from a parsed VCF or SyRI output).

```bash
python parse_variants_for_plotsr.py \
    --input raw_variants.txt \
    --output plotsr_density_input.txt \
    --genome-a-id Reference \
    --genome-b-id Query \
    --variant-types SNP INS DEL 
```

### Step 2: Run PlotSR
Run plotsr with the new SNP arguments.
```bash
./bin/plotsr \
    --sr syri.out \
    --genomes genomes.txt \
    -o output_plot.png \
    --snp-density plotsr_density_input.txt \
    --snp-window-size 100000 \
    --snp-colormap YlGn \
    --snp-bar-thickness 0.12
```

### New Command-Line Arguments

| Argument                | Description                                                                 | Default    |
| ----------------------- | --------------------------------------------------------------------------- | ---------- |
| `--snp-density`         | Path to the simplified variant file created in Step 1.                      | `None`     |
| `--snp-window-size`     | Window size (bp) for calculating variant density.                           | `100000`   |
| `--snp-colormap`        | Matplotlib colormap name (e.g., `YlGn`, `viridis`, `Reds`).                  | `viridis`  |
| `--snp-bar-thickness`   | Vertical thickness of the chromosome heatmap bars (e.g., 0.08).             | `0.08`     |
| `--snp-norm-max-perc`   | Percentile cutoff for color normalization (0-100).                          | `98.0`     |
| `--snp-zero-color`      | A distinct color for windows with zero variants. Use 'none' for transparent.| `lightgray`|

### Modified File List
*   `plotsr/scripts/plotsr.py`: Main execution logic and argument parsing.
*   `plotsr/scripts/func.py`: Core plotting functions, including `calculate_snp_density` and `draw_heatmap_bar`.
*   `parse_variants_for_plotsr.py`: New standalone utility for data preparation.

## License

This project is a fork/modification of **[PlotSR](https://github.com/schneebergerlab/plotsr)**.

*   **Original Code:** Copyright (c) 2021 Manish Goel (MIT License)
*   **Modifications (Variant Density Heatmaps):** Copyright (c) 2025 Paolo Callipo

This software is distributed under the MIT License. See the `LICENSE` file for full details.

## Citation

If you use this modified version of PlotSR in your research, please cite **both** the original tool and our paper describing the enhancements:

1.  **Original PlotSR:**
    > Goel, M., & Schneeberger, K. (2022). plotSR: visualizing structural rearrangements between genomes. *Bioinformatics*, 38(10), 2922-2924. [DOI: 10.1093/bioinformatics/btac196](https://doi.org/10.1093/bioinformatics/btac196)
