#!/usr/bin/env python3
"""
Script: 02_genome_overview_plot.py
Description:
    Generates a multi-page PDF with one page per chromosome pair (PN1..PN19).
    Each page shows:
      - 5 rows: SNP | SV | CG VMR | CHG VMR | CHH VMR density
      - 2 columns: HapA (left) | HapB (right)
    Density computed in fixed-size bins (default 50kb).
    Y-axis normalised globally per track type across all chromosomes.

Usage:
    python 02_genome_overview_plot.py \
        --snp GT_snp.tsv \
        --sv GT_SV_with_end.tsv \
        --cg_vmr CG.vmr.binary.tsv \
        --chg_vmr CHG.vmr.binary.tsv \
        --chh_vmr CHH.vmr.binary.tsv \
        --fai genome.fai \
        --output genome_overview.pdf \
        --bin_size 50000

Author: Paolo Callipo
Date: 2026
Dependencies: python3, matplotlib, numpy
"""

import argparse
import os
import sys
import math
import re
from collections import defaultdict

import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.patches as mpatches
except ImportError:
    sys.exit("Error: 'matplotlib' required. Install via: pip install matplotlib")


# ─────────────────────────────────────────────
# Colour palette
# ─────────────────────────────────────────────
COLORS = {
    'SNP':  '#4C72B0',
    'SV':   '#DD8452',
    'CG':   '#55A868',
    'CHG':  '#8172B2',
    'CHH':  '#C44E52',
}

TRACK_KEYS   = ['SNP', 'SV', 'CG', 'CHG', 'CHH']
TRACK_LABELS = ['SNP', 'SV', 'CG VMR', 'CHG VMR', 'CHH VMR']


# ─────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────

def load_fai(filepath):
    sizes = {}
    with open(filepath) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            sizes[parts[0]] = int(parts[1])
    return sizes


def load_point_variants(filepath):
    positions = defaultdict(list)
    total = 0
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            try:
                positions[parts[0]].append(int(parts[1]))
                total += 1
            except ValueError:
                continue
    print(f"  -> {total:,} SNPs loaded")
    return positions


def load_interval_variants(filepath):
    intervals = defaultdict(list)
    total = 0
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            try:
                intervals[parts[0]].append((int(parts[1]), int(parts[2])))
                total += 1
            except ValueError:
                continue
    print(f"  -> {total:,} SVs loaded")
    return intervals


def load_vmr(filepath):
    bins = defaultdict(list)
    total = 0
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('chrom'):
                continue
            parts = line.strip().split('\t')
            try:
                bins[parts[0]].append((int(parts[1]), int(parts[2])))
                total += 1
            except (ValueError, IndexError):
                continue
    print(f"  -> {total:,} bins loaded from {os.path.basename(filepath)}")
    return bins


# ─────────────────────────────────────────────
# Density computation
# ─────────────────────────────────────────────

def point_density(positions, chrom, chrom_len, bin_size):
    n_bins = max(1, math.ceil(chrom_len / bin_size))
    counts = np.zeros(n_bins)
    for pos in positions.get(chrom, []):
        b = min(int(pos // bin_size), n_bins - 1)
        counts[b] += 1
    return counts


def interval_density(intervals, chrom, chrom_len, bin_size):
    n_bins = max(1, math.ceil(chrom_len / bin_size))
    counts = np.zeros(n_bins)
    for start, end in intervals.get(chrom, []):
        b0 = max(0, int(start // bin_size))
        b1 = min(n_bins - 1, int(end // bin_size))
        counts[b0:b1 + 1] += 1
    return counts


def get_density(key, snp_pos, sv_ivs, cg_vmr, chg_vmr, chh_vmr,
                chrom, chrom_len, bin_size):
    if key == 'SNP':
        return point_density(snp_pos, chrom, chrom_len, bin_size)
    elif key == 'SV':
        return interval_density(sv_ivs, chrom, chrom_len, bin_size)
    elif key == 'CG':
        return interval_density(cg_vmr, chrom, chrom_len, bin_size)
    elif key == 'CHG':
        return interval_density(chg_vmr, chrom, chrom_len, bin_size)
    elif key == 'CHH':
        return interval_density(chh_vmr, chrom, chrom_len, bin_size)


# ─────────────────────────────────────────────
# Chromosome pairing
# ─────────────────────────────────────────────

def get_chr_pairs(chrom_sizes):
    """
    Extract chromosome number and pair HapA/HapB.
    Returns list of (chr_num, hapa_name, hapb_name) sorted by chr_num.
    """
    hapa = {}
    hapb = {}
    for chrom in chrom_sizes:
        m = re.match(r'PN(\d+)_HapA', chrom)
        if m:
            hapa[int(m.group(1))] = chrom
        m = re.match(r'PN(\d+)_HapB', chrom)
        if m:
            hapb[int(m.group(1))] = chrom

    pairs = []
    for num in sorted(set(list(hapa.keys()) + list(hapb.keys()))):
        pairs.append((num, hapa.get(num), hapb.get(num)))
    return pairs


# ─────────────────────────────────────────────
# Compute global maxima for normalisation
# ─────────────────────────────────────────────

def compute_global_maxima(chrom_sizes, snp_pos, sv_ivs, cg_vmr, chg_vmr, chh_vmr,
                           bin_size):
    print("  -> Computing global density maxima for normalisation...")
    global_max = {k: 1.0 for k in TRACK_KEYS}
    for chrom, clen in chrom_sizes.items():
        for key in TRACK_KEYS:
            d = get_density(key, snp_pos, sv_ivs, cg_vmr, chg_vmr, chh_vmr,
                            chrom, clen, bin_size)
            mx = d.max()
            if mx > global_max[key]:
                global_max[key] = mx
    return global_max


# ─────────────────────────────────────────────
# Plot one page
# ─────────────────────────────────────────────

def plot_page(pdf, chr_num, hapa, hapb, chrom_sizes,
              snp_pos, sv_ivs, cg_vmr, chg_vmr, chh_vmr,
              bin_size, global_max):

    n_rows = len(TRACK_KEYS)
    n_cols = 2

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(14, 8),
        gridspec_kw={'wspace': 0.06, 'hspace': 0.35}
    )

    haps = [hapa, hapb]
    hap_labels = ['HapA', 'HapB']

    # Column headers
    for col, (hap, label) in enumerate(zip(haps, hap_labels)):
        chrom_label = hap if hap else 'N/A'
        axes[0][col].set_title(
            f'{label}  —  {chrom_label}',
            fontsize=8, fontweight='bold', color='#333333', pad=5
        )

    for row, (key, track_label) in enumerate(zip(TRACK_KEYS, TRACK_LABELS)):
        color = COLORS[key]
        gmax = global_max[key]

        for col, hap in enumerate(haps):
            ax = axes[row][col]

            if hap is None or hap not in chrom_sizes:
                ax.set_visible(False)
                continue

            clen = chrom_sizes[hap]
            d = get_density(key, snp_pos, sv_ivs, cg_vmr, chg_vmr, chh_vmr,
                            hap, clen, bin_size)
            n_bins = len(d)
            x = np.arange(n_bins) * bin_size / 1e6  # x in Mb

            d_norm = d / gmax

            ax.fill_between(x, 0, d_norm, color=color, alpha=0.80, linewidth=0)
            ax.plot(x, d_norm, color=color, linewidth=0.4, alpha=0.6)

            ax.set_xlim(0, clen / 1e6)
            ax.set_ylim(0, 1.05)

            # Y label on left column only
            if col == 0:
                ax.set_ylabel(track_label, fontsize=7, color=color,
                              fontweight='bold', rotation=90, labelpad=4)
            else:
                ax.set_ylabel('')

            # X axis on bottom row only
            if row == n_rows - 1:
                ax.set_xlabel('Position (Mb)', fontsize=6.5)
                ax.tick_params(axis='x', labelsize=6)
            else:
                ax.set_xticklabels([])
                ax.tick_params(axis='x', length=0)

            ax.set_yticks([])
            ax.tick_params(axis='y', length=0)

            # Subtle grid
            ax.set_facecolor('#fafafa')
            ax.grid(axis='x', color='#dddddd', linewidth=0.4, linestyle='-')

            for spine in ax.spines.values():
                spine.set_linewidth(0.4)
                spine.set_color('#cccccc')

    fig.suptitle(
        f'Chromosome PN{chr_num}  |  bin = {bin_size // 1000}kb',
        fontsize=10, fontweight='bold', y=0.98, color='#222222'
    )

    # Legend
    patches = [mpatches.Patch(color=COLORS[k], label=l)
               for k, l in zip(TRACK_KEYS, TRACK_LABELS)]
    fig.legend(handles=patches, loc='lower center', ncol=5,
               fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.01))

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Per-chromosome genome overview PDF (5 tracks x 2 haplotypes).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--snp", required=True)
    parser.add_argument("--sv", required=True)
    parser.add_argument("--cg_vmr", required=True)
    parser.add_argument("--chg_vmr", required=True)
    parser.add_argument("--chh_vmr", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--output", default="genome_overview.pdf")
    parser.add_argument("--bin_size", type=int, default=50000)
    args = parser.parse_args()

    print("--- Genome Overview Plot ---\n")
    print(f"  Bin size: {args.bin_size:,} bp ({args.bin_size // 1000}kb)\n")

    print("--- Loading data ---")
    chrom_sizes = load_fai(args.fai)
    print(f"  -> {len(chrom_sizes)} chromosomes in FAI")
    snp_pos = load_point_variants(args.snp)
    sv_ivs  = load_interval_variants(args.sv)
    cg_vmr  = load_vmr(args.cg_vmr)
    chg_vmr = load_vmr(args.chg_vmr)
    chh_vmr = load_vmr(args.chh_vmr)

    global_max = compute_global_maxima(
        chrom_sizes, snp_pos, sv_ivs, cg_vmr, chg_vmr, chh_vmr, args.bin_size
    )

    pairs = get_chr_pairs(chrom_sizes)
    print(f"\n--- Generating {len(pairs)} pages ---")

    with PdfPages(args.output) as pdf:
        for chr_num, hapa, hapb in pairs:
            print(f"  -> PN{chr_num}  ({hapa or 'missing'} / {hapb or 'missing'})")
            plot_page(
                pdf, chr_num, hapa, hapb, chrom_sizes,
                snp_pos, sv_ivs, cg_vmr, chg_vmr, chh_vmr,
                args.bin_size, global_max
            )

    print(f"\n  -> PDF saved to {args.output}")


if __name__ == "__main__":
    main()
