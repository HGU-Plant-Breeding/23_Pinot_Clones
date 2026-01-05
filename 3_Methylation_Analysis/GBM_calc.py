#!/usr/bin/env python3
"""
Script Name: GBM_calc.py
Description: Performs Gene Body Methylation (gbM) classification for multiple clones.
             
Methodology:
    1. Loads gene coordinates and methylation data for multiple samples.
    2. Calculates the global genic methylation background (pCG) per clone.
    3. Performs a Binomial Test for every gene in every clone.
    4. Classifies genes as BM (Body Methylated), UM (Unmethylated), or IM (Intermediate).
    5. Aggregates results into matrices (Rows=Genes, Cols=Clones) for downstream plotting.

Input Format:
    --clone Name=Path/To/CG.bed (Repeat for each clone)

Author: Paolo Callipo
Date: 2025
Dependencies: python3, pandas, scipy
"""

import argparse
import sys
import os
from collections import defaultdict, OrderedDict
from bisect import bisect_left, bisect_right

try:
    import pandas as pd
    from scipy.stats import binomtest
except ImportError:
    sys.exit("Error: 'pandas' and 'scipy' are required. Install via: pip install pandas scipy")

# ---------------- Defaults ----------------
DEFAULT_MIN_COV = 2          # Per-site min coverage
DEFAULT_MIN_N_CG = 20        # Min number of CG sites per gene
DEFAULT_ALPHA = 0.05         # P-value threshold
DEFAULT_MIN_COVERAGE_FRACTION = 0.60  # Optional completeness rule

# ---------------- Loaders ----------------

def load_genes_bed(path):
    """
    Reads gene definitions.
    Expected: chrom, start, end, gene_id
    """
    print(f"Loading genes from {os.path.basename(path)}...")
    genes = defaultdict(list)
    try:
        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith('#'): continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 4: continue
                chrom, start, end, gid = p[0], int(p[1]), int(p[2]), p[3]
                genes[gid].append((chrom, start, end))
    except FileNotFoundError:
        sys.exit(f"Error: Gene file '{path}' not found.")
    return genes

def load_cg_percent(path, min_cov):
    """
    Loads CG sites in 'percent' format: chrom, pos, cov, percent
    Returns sorted arrays for bisect lookup.
    """
    pos = defaultdict(list)
    mreads = defaultdict(list)
    cov = defaultdict(list)
    
    try:
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4: continue
                chrom = p[0]
                position = int(p[1])
                c = float(p[2])
                pct = float(p[3])
                
                if c < min_cov: continue
                
                m = int(round(c * (pct/100.0)))
                pos[chrom].append(position)
                mreads[chrom].append(m)
                cov[chrom].append(int(c))
    except FileNotFoundError:
        sys.exit(f"Error: File '{path}' not found.")

    # Sort for binary search
    for ch in list(pos.keys()):
        idx = sorted(range(len(pos[ch])), key=lambda i: pos[ch][i])
        pos[ch]    = [pos[ch][i]    for i in idx]
        mreads[ch] = [mreads[ch][i] for i in idx]
        cov[ch]    = [cov[ch][i]    for i in idx]
    return pos, mreads, cov

def load_cg_counts(path, min_cov):
    """
    Loads CG sites in 'counts' format: chrom, pos, methylated_reads, total_reads
    """
    pos = defaultdict(list)
    mreads = defaultdict(list)
    cov = defaultdict(list)
    
    try:
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 4: continue
                chrom = p[0]
                position = int(p[1])
                m = int(float(p[2]))
                c = int(float(p[3]))
                
                if c < min_cov: continue
                
                pos[chrom].append(position)
                mreads[chrom].append(m)
                cov[chrom].append(c)
    except FileNotFoundError:
        sys.exit(f"Error: File '{path}' not found.")

    for ch in list(pos.keys()):
        idx = sorted(range(len(pos[ch])), key=lambda i: pos[ch][i])
        pos[ch]    = [pos[ch][i]    for i in idx]
        mreads[ch] = [mreads[ch][i] for i in idx]
        cov[ch]    = [cov[ch][i]    for i in idx]
    return pos, mreads, cov

def load_cg_catalog(path):
    """Loads a catalog of ALL theoretical CG sites for completeness check."""
    cat = defaultdict(list)
    if not path: return None
    print(f"Loading CG catalog from {os.path.basename(path)}...")
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2: continue
            chrom, start = p[0], int(p[1])
            cat[chrom].append(start)
    for ch in cat:
        cat[ch].sort()
    return cat

# ---------------- Core Logic ----------------

def sum_counts_in_interval(pos, mreads, cov, start, end):
    """Sums reads within [start, end) using binary search."""
    i = bisect_left(pos, start)
    j = bisect_right(pos, end-1)
    K = 0
    N = 0
    n_sites = j - i
    for k in range(i, j):
        K += mreads[k]
        N += cov[k]
    return K, N, n_sites

def count_catalog_in_interval(cat_pos, start, end):
    """Counts theoretical sites in interval."""
    i = bisect_left(cat_pos, start)
    j = bisect_right(cat_pos, end-1)
    return max(0, j - i)

def classify_one_clone(clone_name, cg_file, genes, fmt, min_cov, min_n_cg, alpha, catalog=None, min_cov_frac=0.60):
    print(f"Processing clone: {clone_name}...")
    
    # Load data
    if fmt == "counts":
        cg_pos, cg_m, cg_cov = load_cg_counts(cg_file, min_cov)
    else:
        cg_pos, cg_m, cg_cov = load_cg_percent(cg_file, min_cov)

    rows = []
    total_K = 0
    total_N = 0

    # Pass 1: Aggregate counts per gene
    for gid, intervals in genes.items():
        K = N = nCG = 0
        cov60_ok = True
        
        for chrom, start, end in intervals:
            if chrom not in cg_pos: continue
            
            k, n, s = sum_counts_in_interval(cg_pos[chrom], cg_m[chrom], cg_cov[chrom], start, end)
            K += k
            N += n
            nCG += s
            
            # Catalog check
            if catalog is not None and chrom in catalog:
                t = count_catalog_in_interval(catalog[chrom], start, end)
                if t > 0 and (s / t) < min_cov_frac:
                    cov60_ok = False

        rows.append({"Gene_ID": gid, "K": K, "N": N, "nCG": nCG, "cov60_ok": cov60_ok})
        total_K += K
        total_N += N

    # Calculate Clone-Specific Background (pCG)
    pCG = (total_K / total_N) if total_N > 0 else 0.0
    
    # Pass 2: Statistical Test & Classification
    PCs = []
    classes = []
    reasons = []
    
    for _, r in pd.DataFrame(rows).iterrows():
        # Low Coverage Filters
        if r['nCG'] < min_n_cg:
            PCs.append(1.0)
            classes.append("LowCov")
            reasons.append(f"nCG<{min_n_cg}")
            continue
        
        if (catalog is not None) and (not r['cov60_ok']):
            PCs.append(1.0)
            classes.append("LowCov")
            reasons.append("coverage_fraction<60%")
            continue
        
        # Binomial Test
        # Null: Gene Methylation == pCG
        if r['N'] > 0:
            pval = binomtest(int(r['K']), int(r['N']), pCG, alternative='greater').pvalue
        else:
            pval = 1.0
            
        PCs.append(pval)
        
        # Classify
        if pval < alpha:
            classes.append("BM")
            reasons.append(f"p<{alpha}")
        elif pval > (1.0 - alpha):
            classes.append("UM")
            reasons.append(f"p>1-{alpha}")
        else:
            classes.append("IM")
            reasons.append("Intermediate")

    # Construct DataFrame
    df = pd.DataFrame(rows)
    df["PCG"] = PCs
    df["Class"] = classes
    df["Reason"] = reasons
    df["Clone"] = clone_name
    
    return df, pCG

# ---------------- Main CLI ----------------

def parse_clones(arglist):
    """Parses 'Name=Path' arguments into a dict."""
    clones = OrderedDict()
    for item in arglist:
        if "=" not in item:
            raise ValueError(f"Invalid format for --clone: '{item}'. Expected 'Name=Path'.")
        name, path = item.split("=", 1)
        clones[name] = path
    return clones

def main():
    parser = argparse.ArgumentParser(
        description="Calculate Gene Body Methylation (gbM) classes and matrices for multiple clones.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("genes_bed", help="BED file with gene coordinates (CDS start to stop).")
    
    # Multi-clone input
    parser.add_argument("--clone", action="append", required=True,
                        help="Input clone data. Format: 'ID=Path/To/File.bed'. Can be used multiple times.")
    
    parser.add_argument("--outdir", default="gbm_output", help="Directory to save output tables.")
    parser.add_argument("--format", choices=["percent","counts"], default="percent",
                        help="Format of input methylation files.")
    
    # Thresholds
    parser.add_argument("--min-cov", type=int, default=DEFAULT_MIN_COV, help="Min per-site coverage.")
    parser.add_argument("--min-n-cg", type=int, default=DEFAULT_MIN_N_CG, help="Min CG sites per gene.")
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA, help="Significance threshold (alpha).")
    parser.add_argument("--catalog-cg", help="Optional BED of all theoretical CGs for coverage check.")
    
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    
    try:
        clones = parse_clones(args.clone)
    except ValueError as e:
        sys.exit(f"Error: {e}")

    genes = load_genes_bed(args.genes_bed)
    catalog = load_cg_catalog(args.catalog_cg)

    # --- Processing ---
    all_tables = []
    backgrounds = []
    
    for name, path in clones.items():
        try:
            df_clone, pCG_val = classify_one_clone(
                name, path, genes, args.format,
                args.min_cov, args.min_n_cg, args.alpha, catalog
            )
            
            # Save individual result
            outfile = os.path.join(args.outdir, f"{name}_classification.tsv")
            df_clone.to_csv(outfile, sep="\t", index=False)
            
            all_tables.append(df_clone)
            backgrounds.append({"Clone": name, "pCG_Background": pCG_val})
            
        except Exception as e:
            print(f"Error processing clone {name}: {e}", file=sys.stderr)

    if not all_tables:
        sys.exit("No clones were processed successfully.")

    print("\nGenerating Aggregate Matrices...")
    big_df = pd.concat(all_tables, ignore_index=True)

    # 1. Class Matrix (BM/UM/IM)
    class_mat = big_df.pivot(index="Gene_ID", columns="Clone", values="Class").fillna("NA")
    class_mat.to_csv(os.path.join(args.outdir, "gbm_class_matrix.tsv"), sep="\t")

    # 2. P-value Matrix
    pcg_mat = big_df.pivot(index="Gene_ID", columns="Clone", values="PCG")
    pcg_mat.to_csv(os.path.join(args.outdir, "pcg_pvalue_matrix.tsv"), sep="\t", float_format="%.5g")

    # 3. Background Stats
    pd.DataFrame(backgrounds).to_csv(
        os.path.join(args.outdir, "background_stats.tsv"), sep="\t", index=False
    )

    print(f"Analysis complete. Results saved in: {args.outdir}/")

if __name__ == "__main__":
    main()
