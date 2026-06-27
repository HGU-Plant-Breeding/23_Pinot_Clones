#!/usr/bin/env python3
"""
Script: 03_5mc_deamination.py
Description:
    Tests whether somatic C->T SNPs are enriched at highly methylated CG sites
    compared to lowly methylated CG sites in the reference genome.

    Key comparison:
      For all CG sites in the reference methylome, ask:
        - Are sites with high methylation (>70%) more likely to carry a C->T SNP
          than sites with low methylation (<30%)?

    This is the correct test for 5mC deamination: not C->T vs other SNPs at CG
    sites, but C->T rate at high-methylation CG sites vs low-methylation CG sites.

    Also reports:
      - Substitution type spectrum
      - Methylation distribution at C->T vs other SNP positions (original comparison)
      - Binned C->T mutation rate across methylation deciles

Usage:
    python 03_5mc_deamination.py \
        --vcf SNP.vcf \
        --methylome 20-13_CG.bed \
        --output 5mc_deamination.tsv

Author: Paolo Callipo
Date: 2026
Dependencies: python3, scipy
"""

import argparse
import os
import sys
from collections import defaultdict

try:
    from scipy.stats import fisher_exact, mannwhitneyu
except ImportError:
    sys.exit("Error: 'scipy' required. Install via: pip install scipy")


COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def complement(base):
    return COMP.get(base.upper(), 'N')


def classify_snp(ref, alt):
    """Classify SNP to pyrimidine strand convention."""
    ref = ref.upper()
    alt = alt.upper()
    if len(ref) != 1 or len(alt) != 1:
        return 'other'
    if ref in ('C', 'T'):
        return f"{ref}>{alt}"
    elif ref in ('G', 'A'):
        return f"{complement(ref)}>{complement(alt)}"
    return 'other'


def load_methylome(filepath):
    """Load reference CG methylome. Returns {(chrom, pos): methylation%}."""
    print(f"  -> Loading reference methylome from {os.path.basename(filepath)}...")
    meth = {}
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            try:
                meth[(parts[0], int(parts[1]))] = float(parts[3])
            except ValueError:
                continue
    print(f"     {len(meth):,} CG sites loaded")
    return meth


def load_vcf(filepath):
    """Parse VCF. Returns {(chrom, pos): subst_type}."""
    print(f"  -> Loading SNPs from {os.path.basename(filepath)}...")
    snps = {}
    subst_counts = defaultdict(int)
    with open(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            ref, alt = parts[3], parts[4]
            if len(ref) != 1 or len(alt) != 1 or alt == '.':
                continue
            try:
                key = (parts[0], int(parts[1]))
            except ValueError:
                continue
            subst = classify_snp(ref, alt)
            snps[key] = subst
            subst_counts[subst] += 1
    total = len(snps)
    print(f"     {total:,} SNPs loaded")
    return snps, subst_counts, total


def print_histogram(values, label, n_bins=10, lo=0, hi=100):
    bin_size = (hi - lo) / n_bins
    counts = defaultdict(int)
    for v in values:
        bucket = min(int((v - lo) / bin_size), n_bins - 1)
        counts[bucket] += 1
    total = len(values)
    max_count = max(counts.values()) if counts else 1
    bar_width = 35
    print(f"\n  {label} (n={total:,}):")
    print(f"  {'Methylation':>14} {'Count':>8}  {'%':>6}  Bar")
    print(f"  {'-'*58}")
    for i in range(n_bins):
        l = lo + i * bin_size
        h = lo + (i + 1) * bin_size
        count = counts[i]
        pct = 100 * count / total if total > 0 else 0
        bar = '█' * int(bar_width * count / max_count) if max_count > 0 else ''
        print(f"  {l:>5.0f}% -{h:>5.0f}%  {count:>8,}  {pct:>5.1f}%  {bar}")


def main():
    parser = argparse.ArgumentParser(
        description="5mC deamination: C->T rate at high vs low methylated CG sites.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--methylome", required=True)
    parser.add_argument("--output", default="5mc_deamination.tsv")
    parser.add_argument("--high_meth", type=float, default=70.0,
                        help="Threshold for highly methylated CG sites (%%)")
    parser.add_argument("--low_meth", type=float, default=30.0,
                        help="Threshold for lowly methylated CG sites (%%)")
    args = parser.parse_args()

    print("--- 5mC Deamination Analysis ---\n")

    methylome = load_methylome(args.methylome)
    snps, subst_counts, total_snps = load_vcf(args.vcf)

    # ── Substitution spectrum ──
    print(f"\n{'='*58}")
    print(f"  SUBSTITUTION TYPE BREAKDOWN")
    print(f"{'='*58}")
    print(f"  {'Type':<10} {'Count':>8}  {'%':>7}")
    print(f"  {'-'*30}")
    for subst in sorted(subst_counts.keys()):
        pct = 100 * subst_counts[subst] / total_snps
        marker = '  <-- putative 5mC deamination' if subst == 'C>T' else ''
        print(f"  {subst:<10} {subst_counts[subst]:>8,}  {pct:>6.1f}%{marker}")

    # ── KEY ANALYSIS: C->T rate across methylation levels ──
    # For every CG site in the reference, record:
    #   - its methylation level
    #   - whether it carries a C->T SNP
    print(f"\n  -> Computing C->T mutation rate across methylation levels...")

    high_with_ct = 0    # high-meth CG sites with C->T SNP
    high_without_ct = 0  # high-meth CG sites without any SNP
    low_with_ct = 0     # low-meth CG sites with C->T SNP
    low_without_ct = 0   # low-meth CG sites without any SNP

    # For decile analysis
    n_deciles = 10
    decile_size = 100 / n_deciles
    decile_ct = defaultdict(int)      # C->T SNPs per decile
    decile_total = defaultdict(int)   # total CG sites per decile

    # For methylation histograms (original comparison)
    ct_meth_vals = []
    other_meth_vals = []

    total_cg = len(methylome)
    processed = 0

    for (chrom, pos), meth_val in methylome.items():
        decile = min(int(meth_val / decile_size), n_deciles - 1)
        decile_total[decile] += 1

        snp_type = snps.get((chrom, pos))
        is_ct = (snp_type == 'C>T')

        if is_ct:
            decile_ct[decile] += 1
            ct_meth_vals.append(meth_val)
        elif snp_type is not None:
            other_meth_vals.append(meth_val)

        if meth_val > args.high_meth:
            if is_ct:
                high_with_ct += 1
            else:
                high_without_ct += 1
        elif meth_val < args.low_meth:
            if is_ct:
                low_with_ct += 1
            else:
                low_without_ct += 1

    # ── Fisher's exact on high vs low meth ──
    high_total = high_with_ct + high_without_ct
    low_total = low_with_ct + low_without_ct
    rate_high = 1000 * high_with_ct / high_total if high_total > 0 else 0
    rate_low = 1000 * low_with_ct / low_total if low_total > 0 else 0

    odds, p_fisher = fisher_exact(
        [[high_with_ct, high_without_ct],
         [low_with_ct,  low_without_ct]],
        alternative='greater'
    )

    print(f"\n{'='*58}")
    print(f"  C->T MUTATION RATE: HIGH vs LOW METHYLATION CG SITES")
    print(f"{'='*58}")
    print(f"  High methylation (>{args.high_meth:.0f}%):")
    print(f"    CG sites:       {high_total:>12,}")
    print(f"    With C->T SNP:  {high_with_ct:>12,}")
    print(f"    Mutation rate:  {rate_high:>11.4f} per 1000 sites")
    print(f"\n  Low methylation (<{args.low_meth:.0f}%):")
    print(f"    CG sites:       {low_total:>12,}")
    print(f"    With C->T SNP:  {low_with_ct:>12,}")
    print(f"    Mutation rate:  {rate_low:>11.4f} per 1000 sites")
    print(f"\n  Fold enrichment (high/low): {rate_high/rate_low:.2f}x" if rate_low > 0 else "")
    print(f"  Fisher's exact p-value:     {p_fisher:.4e}")
    print(f"  Odds ratio:                 {odds:.3f}")
    print(f"{'='*58}")

    # ── Decile analysis ──
    print(f"\n  C->T MUTATION RATE PER METHYLATION DECILE")
    print(f"  {'Methylation':>16} {'CG sites':>12} {'C->T SNPs':>10} {'Rate/1000':>10}")
    print(f"  {'-'*52}")
    for i in range(n_deciles):
        lo = i * decile_size
        hi = (i + 1) * decile_size
        total = decile_total[i]
        ct = decile_ct[i]
        rate = 1000 * ct / total if total > 0 else 0
        print(f"  {lo:>5.0f}% - {hi:>5.0f}%  {total:>12,} {ct:>10,} {rate:>10.4f}")

    # ── Methylation histograms ──
    print_histogram(ct_meth_vals, "CG methylation at C->T SNP positions")
    print_histogram(other_meth_vals, "CG methylation at other SNP positions")

    # ── Write output ──
    print(f"\n  -> Writing output to {args.output}...")
    with open(args.output, 'w') as fout:
        fout.write("chrom\tpos\tmethylation\thas_ct_snp\tsnp_type\n")
        for (chrom, pos), meth_val in methylome.items():
            snp_type = snps.get((chrom, pos), 'none')
            is_ct = 1 if snp_type == 'C>T' else 0
            fout.write(f"{chrom}\t{pos}\t{meth_val:.2f}\t{is_ct}\t{snp_type}\n")
    print(f"  Done.")


if __name__ == "__main__":
    main()
