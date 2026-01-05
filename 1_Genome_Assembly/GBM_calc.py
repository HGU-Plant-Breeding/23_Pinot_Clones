#!/usr/bin/env python3
"""
Script Name: GBM_calc.py
Description: Classifies genes into Body Methylated (BM), Intermediate Methylated (IM), 
             or Unmethylated (UM) states based on CG methylation data.
             
Methodology:
    1. Calculates the global methylation rate (pCG) across all provided genes.
    2. Performs a one-sided Binomial Test for each gene comparison to pCG.
    3. Classifies based on p-value thresholds (default alpha=0.05).
       - BM: p-value < alpha (Significantly methylated)
       - UM: p-value > 1-alpha (Significantly unmethylated)
       - IM: Between thresholds

Author: Paolo Callipo
Date: 2025
Dependencies: pandas, scipy, numpy
"""

import sys
import argparse
import os
from collections import defaultdict
from bisect import bisect_left, bisect_right

try:
    import pandas as pd
    from scipy.stats import binomtest
except ImportError:
    sys.exit("Error: This script requires 'pandas' and 'scipy'. Please install them (pip install pandas scipy).")

# ---------------- Defaults ----------------
DEFAULT_MIN_COV = 1          # Per-site min coverage
DEFAULT_MIN_N_CG = 20        # Min number of CG sites per gene to classify
DEFAULT_GENE_ALPHA = 0.05    # P-value threshold
DEFAULT_MIN_COVERAGE_FRACTION = 0.60  # Fraction of total cytosines covered

# ---------------- Loaders ----------------

def load_genes_bed(path):
    """
    Loads gene coordinates.
    Input BED: chrom, start, end, gene_id
    Returns: dict {gene_id: [(chrom, start, end), ...]}
    """
    print(f"Loading genes from {os.path.basename(path)}...")
    genes = defaultdict(list)
    try:
        with open(path, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith('#'): continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4: continue
                
                chrom, start, end, gid = parts[0], int(parts[1]), int(parts[2]), parts[3]
                genes[gid].append((chrom, start, end))
    except FileNotFoundError:
        sys.exit(f"Error: Gene file '{path}' not found.")
    return genes

def load_cg_sites_percent(path, min_cov):
    """
    Loads methylation data in Percent format.
    Format: chrom, pos, coverage, percent_methylated (0-100)
    Returns: dicts for pos, methylated_reads, coverage
    """
    print(f"Loading CG sites (Percent Format) from {os.path.basename(path)}...")
    pos = defaultdict(list)
    mreads = defaultdict(list)
    cov = defaultdict(list)
    
    try:
        with open(path, 'r') as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4: continue
                
                chrom = parts[0]
                p = int(parts[1])
                c = float(parts[2]) # Coverage
                pct = float(parts[3]) # Percent
                
                if c < min_cov: continue
                
                # Convert percent back to estimated read counts
                m = int(round(c * (pct/100.0)))
                
                pos[chrom].append(p)
                mreads[chrom].append(m)
                cov[chrom].append(int(c))
    except FileNotFoundError:
        sys.exit(f"Error: CG file '{path}' not found.")

    # Sort arrays by position for bisect usage
    for ch in list(pos.keys()):
        idx = sorted(range(len(pos[ch])), key=lambda i: pos[ch][i])
        pos[ch]    = [pos[ch][i]    for i in idx]
        mreads[ch] = [mreads[ch][i] for i in idx]
        cov[ch]    = [cov[ch][i]    for i in idx]
    return pos, mreads, cov

def load_cg_sites_counts(path, min_cov):
    """
    Loads methylation data in Count format.
    Format: chrom, pos, methylated_reads, coverage
    """
    print(f"Loading CG sites (Count Format) from {os.path.basename(path)}...")
    pos = defaultdict(list)
    mreads = defaultdict(list)
    cov = defaultdict(list)
    
    try:
        with open(path, 'r') as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4: continue
                
                chrom = parts[0]
                p = int(parts[1])
                m = int(float(parts[2])) # Methylated reads
                c = int(float(parts[3])) # Total reads (coverage)
                
                if c < min_cov: continue
                
                pos[chrom].append(p)
                mreads[chrom].append(m)
                cov[chrom].append(c)
    except FileNotFoundError:
        sys.exit(f"Error: CG file '{path}' not found.")

    for ch in list(pos.keys()):
        idx = sorted(range(len(pos[ch])), key=lambda i: pos[ch][i])
        pos[ch]    = [pos[ch][i]    for i in idx]
        mreads[ch] = [mreads[ch][i] for i in idx]
        cov[ch]    = [cov[ch][i]    for i in idx]
    return pos, mreads, cov

def load_cg_catalog(path):
    """
    Loads a catalog of ALL theoretical CG positions (genome-wide).
    Used to check if a gene has enough sequencing data coverage.
    """
    print(f"Loading CG catalog from {os.path.basename(path)}...")
    cat = defaultdict(list)
    try:
        with open(path, 'r') as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2: continue # Need at least Chrom, Start
                chrom, start = parts[0], int(parts[1])
                cat[chrom].append(start)
    except FileNotFoundError:
        sys.exit(f"Error: Catalog file '{path}' not found.")

    for ch in cat:
        cat[ch].sort()
    return cat

# ---------------- Counting Logic ----------------

def sum_counts_in_interval(pos_list, mreads_list, cov_list, start, end):
    """
    Sums methylated reads and coverage for all sites within [start, end).
    Uses binary search (bisect) for speed.
    """
    i = bisect_left(pos_list, start)
    j = bisect_right(pos_list, end-1)
    
    m_sum = 0
    n_sum = 0
    n_sites = j - i
    
    # Iterate only the slice relevant to the gene
    for k in range(i, j):
        m_sum += mreads_list[k]
        n_sum += cov_list[k]
        
    return m_sum, n_sum, n_sites

def count_catalog_in_interval(cat_pos_list, start, end):
    """Counts how many theoretical CGs exist in the interval."""
    i = bisect_left(cat_pos_list, start)
    j = bisect_right(cat_pos_list, end-1)
    return max(0, j - i)

# ---------------- Main ----------------

def main():
    parser = argparse.ArgumentParser(
        description="Classify genes (BM/IM/UM) using a binomial test on CG methylation data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("genes_bed", help="BED file with gene coordinates (chrom, start, end, gene_id).")
    parser.add_argument("cg", help="CG sites file (format depends on --format argument).")
    parser.add_argument("output", help="Output TSV filename.")

    parser.add_argument("--format", choices=["percent","counts"], default="percent",
                    help="Format of CG input: 'percent' (Chrom, Pos, Cov, %) or 'counts' (Chrom, Pos, MethCounts, Cov).")
    parser.add_argument("--min-cov", type=int, default=DEFAULT_MIN_COV,
                    help="Minimum read coverage per site to include.")
    parser.add_argument("--min-n-cg", type=int, default=DEFAULT_MIN_N_CG,
                    help="Minimum number of covered CG sites required per gene.")
    parser.add_argument("--gene-alpha", type=float, default=DEFAULT_GENE_ALPHA,
                    help="Significance threshold (alpha) for the binomial test.")
    parser.add_argument("--catalog-cg", help="Optional BED file of all theoretical CG positions in the genome. Used to enforce coverage completeness.")
    parser.add_argument("--min-coverage-fraction", type=float, default=DEFAULT_MIN_COVERAGE_FRACTION,
                    help="If catalog provided: Gene must have data for this fraction of theoretical CGs (0.0-1.0).")
    
    args = parser.parse_args()

    # 1. Load Data
    genes = load_genes_bed(args.genes_bed)
    
    if args.format == "counts":
        cg_pos, cg_m, cg_cov = load_cg_sites_counts(args.cg, args.min_cov)
    else:
        cg_pos, cg_m, cg_cov = load_cg_sites_percent(args.cg, args.min_cov)

    cat = None
    if args.catalog_cg:
        cat = load_cg_catalog(args.catalog_cg)
    else:
        print("[INFO] No CG catalog provided. Skipping 'coverage fraction' filter.")

    # 2. Aggregation Pass
    # Calculate total K (methylated) and N (total) for the background calculation
    print("Aggregating counts per gene...")
    rows = []
    total_K = 0
    total_N = 0

    for gid, intervals in genes.items():
        K = 0   # Methylated reads sum
        N = 0   # Total reads sum
        nCG = 0 # Number of cytosines covered
        
        cov60_ok = True
        
        for chrom, start, end in intervals:
            if chrom not in cg_pos:
                continue
            
            # Sum experimental data
            k, n, s = sum_counts_in_interval(cg_pos[chrom], cg_m[chrom], cg_cov[chrom], start, end)
            K += k
            N += n
            nCG += s
            
            # Check theoretical coverage if catalog is present
            if cat is not None and chrom in cat:
                t = count_catalog_in_interval(cat[chrom], start, end)
                if t > 0 and (s / t) < args.min_coverage_fraction:
                    cov60_ok = False

        rows.append({"Gene_ID": gid, "K": K, "N": N, "nCG": nCG, "cov60_ok": cov60_ok})
        total_K += K
        total_N += N

    # 3. Calculate Background Probability (pCG)
    # This represents the average methylation level of the provided gene set
    pCG = (total_K / total_N) if total_N > 0 else 0.0
    print(f"Global Genic Background (pCG): {pCG:.5f}")

    # 4. Statistical Testing & Classification
    print("Performing binomial tests...")
    
    final_rows = []
    
    for r in rows:
        gene_id = r["Gene_ID"]
        K, N, nCG = r["K"], r["N"], r["nCG"]
        
        # Filter: Low Data
        if nCG < args.min_n_cg:
            final_rows.append([gene_id, "LowCov", f"nCG<{args.min_n_cg}", nCG, K, N, 1.0])
            continue
        
        # Filter: Low Completeness
        if args.catalog_cg and not r["cov60_ok"]:
            final_rows.append([gene_id, "LowCov", "coverage_fraction<threshold", nCG, K, N, 1.0])
            continue

        # Binomial Test
        # Null Hypothesis: Gene methylation rate == Global Background (pCG)
        # Alternative: Gene methylation rate > Global Background
        if N > 0:
            result = binomtest(k=int(K), n=int(N), p=pCG, alternative='greater')
            pval = result.pvalue
        else:
            pval = 1.0

        # Classification Logic
        # If p < alpha: Methylation is significantly HIGHER than background -> BM
        # If p > 1-alpha: Methylation is significantly LOWER than background -> UM
        # Else: IM
        
        classification = "IM"
        reason = "Intermediate"
        
        if pval < args.gene_alpha:
            classification = "BM"
            reason = f"p<{args.gene_alpha}"
        elif pval > (1.0 - args.gene_alpha):
            classification = "UM"
            reason = f"p>1-{args.gene_alpha}"
        
        final_rows.append([gene_id, classification, reason, nCG, K, N, pval])

    # 5. Output
    df = pd.DataFrame(final_rows, columns=["Gene_ID", "Class", "Reason", "nCG", "K", "N", "PCG"])
    
    # Save to TSV
    df.to_csv(args.output, sep="\t", index=False, float_format="%.5g")
    
    print("\n--- Classification Summary ---")
    print(df["Class"].value_counts())
    print(f"Results saved to {args.output}")

if __name__ == "__main__":
    main()
