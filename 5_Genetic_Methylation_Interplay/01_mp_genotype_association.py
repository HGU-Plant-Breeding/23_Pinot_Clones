#!/usr/bin/env python3
"""
Script: 01_mp_genotype_association.py
Description:
    Computes phi coefficient associations between VMR bins and genetic variants
    (SNPs, SVs, or combined) within a cis window.

    For each VMR bin:
      - SNPs: variant point position falls within [bin_start - flank, bin_end + flank]
      - SVs:  SV interval [sv_start, sv_end] overlaps [bin_start - flank, bin_end + flank]
      - For each VMR-variant pair, computes phi coefficient and Fisher's exact p-value
        using pairwise exclusion of missing data (requires >= min_samples valid pairs)
      - Reports the single best association (highest |phi|) per VMR

    Genotype recoding:
      - 0/0 -> 0
      - 0/1 or 1/1 -> 1
      - ./. -> missing

    VMR recoding:
      - 0 -> 0, 1 -> 1, . -> missing

    SNP file format:  chrom, pos, gt x23  (no header)
    SV file format:   chrom, start, end, svtype, gt x23  (no header, from bcftools query)

    Output columns:
      chrom, bin_start, bin_end, vmr_maf, n_variants_tested,
      best_variant_start, best_variant_end, variant_type, phi, fisher_p, n_samples

    Runs all three contexts (CG, CHG, CHH) x three variant types (SNP, SV, SNP+SV)
    producing 9 output files.

Usage:
    python 01_vmr_genotype_association.py \
        --vmr_dir ./ \
        --snp GT_snp.tsv \
        --sv GT_SV_with_end.tsv \
        --output_dir ./associations \
        --flank 1000 \
        --min_samples 18

Author: Paolo Callipo
Date: 2026
Dependencies: python3, scipy
"""

import argparse
import os
import sys
import math
from collections import defaultdict

try:
    from scipy.stats import fisher_exact
except ImportError:
    sys.exit("Error: 'scipy' required. Install via: pip install scipy")


# ─────────────────────────────────────────────
# Genotype / VMR parsing
# ─────────────────────────────────────────────

def recode_genotype(gt):
    """Recode VCF-style genotype to 0/1/None."""
    gt = gt.strip()
    if gt in ('./.', '.', ''):
        return None
    if gt == '0/0':
        return 0
    if gt in ('0/1', '1/0', '1/1'):
        return 1
    return None


def recode_vmr(val):
    """Recode VMR binary value to 0/1/None."""
    val = val.strip()
    if val == '.':
        return None
    if val == '0':
        return 0
    if val == '1':
        return 1
    return None


# ─────────────────────────────────────────────
# Statistics
# ─────────────────────────────────────────────

def compute_contingency(vmr_vec, gt_vec):
    """
    Build 2x2 contingency table using pairwise exclusion.
    Returns (n11, n10, n01, n00, n) or None if no valid pairs.
    """
    n11 = n10 = n01 = n00 = 0
    for v, g in zip(vmr_vec, gt_vec):
        if v is None or g is None:
            continue
        if v == 1 and g == 1:
            n11 += 1
        elif v == 1 and g == 0:
            n10 += 1
        elif v == 0 and g == 1:
            n01 += 1
        else:
            n00 += 1
    n = n11 + n10 + n01 + n00
    return (n11, n10, n01, n00, n) if n > 0 else None


def phi_from_table(n11, n10, n01, n00):
    """Compute phi coefficient from 2x2 contingency table."""
    r1 = n11 + n10  # VMR=1
    r0 = n01 + n00  # VMR=0
    c1 = n11 + n01  # GT=1
    c0 = n10 + n00  # GT=0
    denom = math.sqrt(r1 * r0 * c1 * c0)
    if denom == 0:
        return 0.0
    return (n11 * n00 - n10 * n01) / denom


# ─────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────

def load_snps(filepath):
    """
    Load SNP file: chrom, pos, gt x23
    Returns {chrom: [(start, end, type, gt_vec)]} sorted by start.
    SNPs: start == end (point variants).
    """
    print(f"  -> Loading SNPs from {os.path.basename(filepath)}...")
    variants = defaultdict(list)
    total = 0
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            gt_vec = [recode_genotype(g) for g in parts[2:]]
            variants[chrom].append((pos, pos, 'SNP', gt_vec))
            total += 1
    for chrom in variants:
        variants[chrom].sort(key=lambda x: x[0])
    print(f"     {total:,} SNPs loaded")
    return variants


def load_svs(filepath):
    """
    Load SV file: chrom, start, end, svtype, gt x23
    Returns {chrom: [(start, end, svtype, gt_vec)]} sorted by start.
    """
    print(f"  -> Loading SVs from {os.path.basename(filepath)}...")
    variants = defaultdict(list)
    total = 0
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            svtype = parts[3]
            gt_vec = [recode_genotype(g) for g in parts[4:]]
            variants[chrom].append((start, end, svtype, gt_vec))
            total += 1
    for chrom in variants:
        variants[chrom].sort(key=lambda x: x[0])
    print(f"     {total:,} SVs loaded")
    return variants


def load_vmr(filepath):
    """Load VMR binary matrix."""
    print(f"  -> Loading VMRs from {os.path.basename(filepath)}...")
    vmrs = []
    with open(filepath) as f:
        for line in f:
            if not line.strip() or line.startswith('chrom'):
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            vmr_vec = [recode_vmr(v) for v in parts[3:]]
            zeros = vmr_vec.count(0)
            ones = vmr_vec.count(1)
            total = zeros + ones
            maf = min(ones, zeros) / total if total > 0 else None
            vmrs.append((chrom, start, end, maf, vmr_vec))
    print(f"     {len(vmrs):,} VMR bins loaded")
    return vmrs


def get_variants_in_window(variants, chrom, win_start, win_end):
    """
    Return all variants whose interval overlaps [win_start, win_end].
    Overlap condition: var_start <= win_end AND var_end >= win_start.
    Works for both SNPs (start==end) and SVs (start < end).
    """
    chrom_vars = variants.get(chrom, [])
    result = []
    for var_start, var_end, vtype, gt_vec in chrom_vars:
        if var_start <= win_end and var_end >= win_start:
            result.append((var_start, var_end, vtype, gt_vec))
        elif var_start > win_end:
            break
    return result


def merge_variant_dicts(snp_variants, sv_variants):
    """Merge SNP and SV variant dicts, sorted by start position."""
    merged = defaultdict(list)
    for chrom in set(list(snp_variants.keys()) + list(sv_variants.keys())):
        combined = snp_variants.get(chrom, []) + sv_variants.get(chrom, [])
        merged[chrom] = sorted(combined, key=lambda x: x[0])
    return merged


# ─────────────────────────────────────────────
# Main association logic
# ─────────────────────────────────────────────

def run_association(vmrs, variants, min_samples, flank, output_path, label):
    """Run VMR-variant association for one context/variant combination."""
    print(f"  -> Running association ({label})...")

    total_vmrs = len(vmrs)
    vmrs_tested = 0
    vmrs_no_variants = 0

    with open(output_path, 'w') as fout:
        fout.write("chrom\tbin_start\tbin_end\tvmr_maf\t"
                   "n_variants_tested\tbest_variant_start\tbest_variant_end\t"
                   "variant_type\tphi\tfisher_p\tn_samples\n")

        for chrom, start, end, maf, vmr_vec in vmrs:
            win_start = start - flank
            win_end = end + flank
            maf_str = f"{maf:.4f}" if maf is not None else "NA"

            candidates = get_variants_in_window(variants, chrom, win_start, win_end)

            if not candidates:
                vmrs_no_variants += 1
                fout.write(f"{chrom}\t{start}\t{end}\t{maf_str}\t"
                           f"0\tNA\tNA\tNA\tNA\tNA\tNA\n")
                continue

            vmrs_tested += 1
            best_phi = None
            best_p = None
            best_start = None
            best_end = None
            best_type = None
            best_n = None
            n_tested = 0

            for var_start, var_end, vtype, gt_vec in candidates:
                result = compute_contingency(vmr_vec, gt_vec)
                if result is None:
                    continue
                n11, n10, n01, n00, n = result
                if n < min_samples:
                    continue

                n_tested += 1
                phi = phi_from_table(n11, n10, n01, n00)
                _, p = fisher_exact([[n11, n10], [n01, n00]])

                if best_phi is None or abs(phi) > abs(best_phi):
                    best_phi = phi
                    best_p = p
                    best_start = var_start
                    best_end = var_end
                    best_type = vtype
                    best_n = n

            if best_phi is not None:
                fout.write(f"{chrom}\t{start}\t{end}\t{maf_str}\t"
                           f"{n_tested}\t{best_start}\t{best_end}\t"
                           f"{best_type}\t{best_phi:.4f}\t{best_p:.6e}\t{best_n}\n")
            else:
                fout.write(f"{chrom}\t{start}\t{end}\t{maf_str}\t"
                           f"{n_tested}\tNA\tNA\tNA\tNA\tNA\tNA\n")

    print(f"     VMRs total:          {total_vmrs:,}")
    print(f"     VMRs with variants:  {vmrs_tested:,}")
    print(f"     VMRs no variants:    {vmrs_no_variants:,}")
    print(f"     Output: {output_path}")


# ─────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="VMR-genotype phi coefficient association analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vmr_dir", required=True,
                        help="Directory containing CG/CHG/CHH.vmr.binary.tsv files")
    parser.add_argument("--snp", required=True,
                        help="SNP genotype file (chrom, pos, gt x23)")
    parser.add_argument("--sv", required=True,
                        help="SV file with end coords (chrom, start, end, svtype, gt x23)")
    parser.add_argument("--output_dir", required=True, help="Output directory")
    parser.add_argument("--flank", type=int, default=1000,
                        help="Flanking window around VMR bin edges (bp)")
    parser.add_argument("--min_samples", type=int, default=18,
                        help="Minimum valid samples per VMR-variant pair")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("  VMR-Genotype Association Analysis")
    print(f"  Flank:       ±{args.flank}bp from bin edges")
    print(f"  Min samples: {args.min_samples}")
    print(f"  SV lookup:   interval overlap")
    print(f"  SNP lookup:  point position")
    print("=" * 60)

    print("\n--- Loading variant files ---")
    snp_variants = load_snps(args.snp)
    sv_variants = load_svs(args.sv)

    print("  -> Merging SNP + SV...")
    merged_variants = merge_variant_dicts(snp_variants, sv_variants)
    print(f"     {sum(len(v) for v in merged_variants.values()):,} total variants")

    variant_sets = [
        ("snp",    snp_variants),
        ("sv",     sv_variants),
        ("snp_sv", merged_variants),
    ]

    contexts = ["CG", "CHG", "CHH"]

    for ctx in contexts:
        vmr_file = os.path.join(args.vmr_dir, f"{ctx}.vmr.binary.tsv")
        if not os.path.exists(vmr_file):
            print(f"\n  WARNING: {vmr_file} not found, skipping {ctx}")
            continue

        print(f"\n--- Context: {ctx} ---")
        vmrs = load_vmr(vmr_file)

        for vtype, variants in variant_sets:
            label = f"{ctx}.{vtype}"
            output_path = os.path.join(args.output_dir, f"{label}.associations.tsv")
            run_association(vmrs, variants, args.min_samples,
                            args.flank, output_path, label)

    print("\n" + "=" * 60)
    print("  Done!")
    print(f"  Output files in: {args.output_dir}/")
    print("=" * 60)


if __name__ == "__main__":
    main()
