#!/usr/bin/env python3
"""
Script Name: GBM_classify_CG_CHG.py
Description: Classifies genes as gbM / teM / UM / Unclassified from ONT methylation
             data using only CG and CHG contexts. Adapted for site-averaged
             aggregation of ONT per-site calls to avoid pseudoreplication.

    Per-gene test (for CG and CHG):
        For each cytosine site with cv > 0, compute the site-level methylation
        fraction (mr / cv). Across all sites of that context in the gene:
            avg_frac = unweighted mean of site fractions
            eff_N    = number of sites with cv > 0
            eff_K    = round(avg_frac * eff_N)
        Null p = genome-wide CDS average methylation
        One-sided binomial test (greater) on (eff_K, eff_N) -> P-value per gene
        BH-FDR correction across genes, per context

    Classification (CHH removed for noise reduction):
        teM if (CHG significant AND CHG frac >= min_chg_frac)
        gbM if CG significant AND not teM AND adequate coverage in CG & CHG
        UM  if CG not significant AND not teM AND adequate coverage in CG & CHG
              AND CG_frac <= max_cg_frac_um   (effect-size guard)
        Unclassified otherwise

    The max_cg_frac_um guard prevents genes with intermediate CG methylation
    (below background, but not biologically unmethylated) from being called UM
    in genomes with high background pCG. Without it, the one-sided "greater"
    binomial test lumps together truly unmethylated genes and intermediate
    genes, since both fail significance.

Usage:
    python GBM_classify_CG_CHG.py genes_cds.bed cg.bed chg.bed output_prefix [options]

Dependencies: pandas, scipy (>=1.11), matplotlib, numpy
"""

import sys
import argparse
import os
from collections import defaultdict
from bisect import bisect_left, bisect_right

try:
    import numpy as np
    import pandas as pd
    from scipy.stats import binomtest
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("Error: pip install pandas scipy matplotlib numpy")

# BH correction
try:
    from scipy.stats import false_discovery_control as _fdr
    def bh_correct(pvals):
        pvals = np.asarray(pvals, dtype=float)
        mask = ~np.isnan(pvals)
        out = np.full_like(pvals, np.nan)
        if mask.sum() > 0:
            out[mask] = _fdr(pvals[mask], method="bh")
        return out
except ImportError:
    try:
        from statsmodels.stats.multitest import multipletests
        def bh_correct(pvals):
            pvals = np.asarray(pvals, dtype=float)
            mask = ~np.isnan(pvals)
            out = np.full_like(pvals, np.nan)
            if mask.sum() > 0:
                _, q, _, _ = multipletests(pvals[mask], method="fdr_bh")
                out[mask] = q
            return out
    except ImportError:
        sys.exit("Error: need scipy>=1.11 or statsmodels for BH correction")

# ---------------- Defaults ----------------
DEFAULT_MIN_COV         = 3
DEFAULT_MIN_N_CG        = 20
DEFAULT_MIN_N_CHG       = 20
DEFAULT_ALPHA           = 0.05
DEFAULT_MIN_CHG_FRAC    = 0.10
DEFAULT_MAX_CG_FRAC_UM  = 0.10   # effect-size guard for UM call

# ---------------- Loaders ----------------

def load_genes_bed(path):
    print(f"Loading genes from {os.path.basename(path)}...")
    genes = defaultdict(list)
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, gid = parts[0], int(parts[1]), int(parts[2]), parts[3]
            genes[gid].append((chrom, start, end))
    print(f"  {len(genes):,} genes loaded.")
    return genes


def load_methylation(path, min_cov, fmt):
    print(f"  Loading ({fmt}): {os.path.basename(path)}...")
    pos = defaultdict(list)
    mr  = defaultdict(list)
    cv  = defaultdict(list)
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            p = int(parts[1])
            if fmt == "percent":
                c = float(parts[2])
                pct = float(parts[3])
                if c < min_cov:
                    continue
                c = int(round(c))
                m = int(round(c * (pct / 100.0)))
            else:
                m = int(float(parts[2]))
                c = int(float(parts[3]))
                if c < min_cov:
                    continue
            pos[chrom].append(p)
            mr[chrom].append(m)
            cv[chrom].append(c)

    for ch in list(pos.keys()):
        idx = sorted(range(len(pos[ch])), key=lambda i: pos[ch][i])
        pos[ch] = [pos[ch][i] for i in idx]
        mr[ch]  =[mr[ch][i]  for i in idx]
        cv[ch]  = [cv[ch][i]  for i in idx]
    total = sum(len(v) for v in pos.values())
    print(f"    {total:,} sites kept (min_cov={min_cov}).")
    return pos, mr, cv

# ---------------- Interval site collection ----------------

def sum_interval(pos_list, mr_list, cv_list, start, end):
    i = bisect_left(pos_list, start)
    j = bisect_right(pos_list, end - 1)
    site_fracs =[]
    raw_K = 0
    raw_N = 0
    for k in range(i, j):
        c = cv_list[k]
        if c <= 0:
            continue
        m = mr_list[k]
        site_fracs.append(m / c)
        raw_K += m
        raw_N += c
    return site_fracs, raw_K, raw_N

# ---------------- Per-gene aggregation ----------------

def aggregate_gene(intervals, pos_data, mr_data, cv_data):
    site_fracs =[]
    raw_K = 0
    raw_N = 0
    for chrom, start, end in intervals:
        if chrom not in pos_data:
            continue
        sf, k, n = sum_interval(
            pos_data[chrom], mr_data[chrom], cv_data[chrom], start, end
        )
        site_fracs.extend(sf)
        raw_K += k
        raw_N += n

    eff_N = len(site_fracs)
    if eff_N > 0:
        avg_frac = float(np.mean(site_fracs))
        eff_K = int(round(avg_frac * eff_N))
    else:
        avg_frac = 0.0
        eff_K = 0

    return {
        "K":       eff_K,
        "N":       eff_N,
        "n_sites": eff_N,
        "frac":    avg_frac,
        "raw_K":   raw_K,
        "raw_N":   raw_N,
    }


def aggregate_all(genes, pos_data, mr_data, cv_data):
    return {gid: aggregate_gene(ivs, pos_data, mr_data, cv_data)
            for gid, ivs in genes.items()}

# ---------------- Background ----------------

def compute_background(gene_agg):
    total_K = sum(a["K"] for a in gene_agg.values())
    total_N = sum(a["N"] for a in gene_agg.values())
    return (total_K / total_N) if total_N > 0 else 0.0

# ---------------- Binomial + BH ----------------

def run_binomial_with_bh(gene_agg, background_p, min_n):
    gids = list(gene_agg.keys())
    pvals = np.full(len(gids), np.nan)

    for i, gid in enumerate(gids):
        a = gene_agg[gid]
        if a["n_sites"] < min_n or a["N"] == 0:
            continue
        pvals[i] = binomtest(
            int(a["K"]), int(a["N"]), background_p, alternative="greater"
        ).pvalue

    qvals = bh_correct(pvals)

    out = {}
    for i, gid in enumerate(gids):
        a = gene_agg[gid]
        out[gid] = {
            "pval":    None if np.isnan(pvals[i]) else float(pvals[i]),
            "qval":    None if np.isnan(qvals[i]) else float(qvals[i]),
            "frac":    a["frac"],
            "n_sites": a["n_sites"],
            "K":       a["K"],
            "N":       a["N"],
        }
    return out

# ---------------- Classification ----------------

def classify(cg_res, chg_res, alpha, min_n_cg, min_n_chg, min_chg_frac,
             max_cg_frac_um):
    def sig(q):
        return q is not None and q < alpha

    def ok_n(res, min_n):
        return res["n_sites"] >= min_n

    empty = {"pval": None, "qval": None, "frac": None, "n_sites": 0, "K": 0, "N": 0}
    all_genes = set(cg_res) | set(chg_res)
    rows =[]

    for gid in all_genes:
        cg  = cg_res.get(gid, empty)
        chg = chg_res.get(gid, empty)

        is_tem = sig(chg["qval"]) and chg["frac"] is not None and chg["frac"] >= min_chg_frac
        has_all_cov = (ok_n(cg, min_n_cg) and ok_n(chg, min_n_chg))

        if has_all_cov and is_tem:
            cl = "teM"
        elif has_all_cov and sig(cg["qval"]):
            cl = "gbM"
        elif (has_all_cov and not sig(cg["qval"])
              and cg["frac"] is not None and cg["frac"] <= max_cg_frac_um):
            cl = "UM"
        else:
            cl = "Unclassified"

        def fmt(x, sp="{:.5g}"):
            return sp.format(x) if x is not None else "NA"

        rows.append({
            "Gene_ID":   gid,
            "Classification": cl,
            "CG_nSites":  cg["n_sites"],  "CG_K":  cg["K"],  "CG_N":  cg["N"],
            "CG_frac":    fmt(cg["frac"],  "{:.4f}"),
            "CG_pval":    fmt(cg["pval"]),
            "CG_qval":    fmt(cg["qval"]),
            "CHG_nSites": chg["n_sites"], "CHG_K": chg["K"], "CHG_N": chg["N"],
            "CHG_frac":   fmt(chg["frac"], "{:.4f}"),
            "CHG_pval":   fmt(chg["pval"]),
            "CHG_qval":   fmt(chg["qval"])
        })

    return pd.DataFrame(rows)

# ---------------- Plotting ----------------

def _class_fracs(df, frac_col):
    vals = pd.to_numeric(df[frac_col], errors="coerce")
    out = {}
    for cl in["gbM", "teM", "UM"]:
        mask = (df["Classification"] == cl) & vals.notna()
        out[cl] = vals[mask].values
    return out


def plot_distribution(df, frac_col, context, threshold, outpath, zoom_max=None, show_threshold=True):
    colors = {"gbM": "#4A90D9", "teM": "#E57373", "UM": "#BDBDBD"}
    data = _class_fracs(df, frac_col)
    counts = {cl: len(v) for cl, v in data.items()}

    if zoom_max is None:
        all_nonzero = np.concatenate([v for v in data.values() if len(v) > 0])
        zoom_max = min(0.30, np.quantile(all_nonzero, 0.99))

    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    fig.suptitle(f"{context} Methylation Fraction Distribution by Gene Class",
                 fontsize=13, fontweight="bold")

    for ax, xmax, title in [(axes[0], 1.0, "Full range"),
                             (axes[1], zoom_max, f"Zoomed: 0 \u2013 {zoom_max:.2f}")]:
        bins = np.linspace(0, xmax, 120)
        for cl in["UM", "gbM", "teM"]:
            if len(data[cl]) == 0:
                continue
            ax.hist(data[cl], bins=bins, alpha=0.55, color=colors[cl],
                    density=True, label=f"{cl} (n={counts[cl]:,})")
        if show_threshold and threshold is not None:
            ax.axvline(threshold, color="black", linestyle="--", linewidth=1.2,
                       label=f"Threshold = {threshold}")
        ax.set_xlabel(f"{context} methylation fraction")
        ax.set_ylabel("Density")
        ax.set_xlim(0, xmax)
        ax.set_title(title)
        ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def plot_cg_sanity(df, outpath, threshold=None):
    colors = {"gbM": "#4A90D9", "teM": "#E57373", "UM": "#BDBDBD"}
    data = _class_fracs(df, "CG_frac")
    counts = {cl: len(v) for cl, v in data.items()}

    fig, ax = plt.subplots(figsize=(10, 5.5))
    bins = np.linspace(0, 1, 120)
    for cl in["UM", "gbM", "teM"]:
        if len(data[cl]) == 0:
            continue
        ax.hist(data[cl], bins=bins, alpha=0.55, color=colors[cl],
                density=True, label=f"{cl} (n={counts[cl]:,})")
    if threshold is not None:
        ax.axvline(threshold, color="black", linestyle="--", linewidth=1.2,
                   label=f"UM cap = {threshold}")
    ax.set_xlabel("CG methylation fraction")
    ax.set_ylabel("Density")
    ax.set_xlim(0, 1)
    ax.set_title("CG Methylation Fraction by Gene Class (sanity check)",
                 fontsize=12, fontweight="bold")
    ax.legend(loc="upper center", fontsize=10)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def plot_classification_bar(df, outpath):
    categories = ["gbM", "teM", "UM", "Unclassified"]
    colors_cat =["#4A90D9", "#E57373", "#BDBDBD", "#EEEEEE"]
    counts = [(df["Classification"] == c).sum() for c in categories]
    total = len(df)

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(categories, counts, color=colors_cat, edgecolor="white")
    for bar, n in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(counts) * 0.01,
                f"{n:,}\n({100*n/total:.1f}%)",
                ha="center", va="bottom", fontsize=10)
    ax.set_ylabel("Number of genes")
    ax.set_title("Gene classification summary", fontweight="bold")
    ax.set_ylim(0, max(counts) * 1.18)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")

# ---------------- Main ----------------

def main():
    parser = argparse.ArgumentParser(
        description="gbM/teM/UM classification (site-averaged effective counts). CHH removed.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("genes_bed")
    parser.add_argument("cg_file")
    parser.add_argument("chg_file")
    parser.add_argument("output_prefix")

    parser.add_argument("--format",        choices=["percent", "counts"], default="percent")
    parser.add_argument("--min-cov",       type=int,   default=DEFAULT_MIN_COV)
    parser.add_argument("--min-n-cg",      type=int,   default=DEFAULT_MIN_N_CG)
    parser.add_argument("--min-n-chg",     type=int,   default=DEFAULT_MIN_N_CHG)
    parser.add_argument("--alpha",         type=float, default=DEFAULT_ALPHA,
                        help="FDR threshold on BH-corrected q-values")
    parser.add_argument("--min-chg-frac",  type=float, default=DEFAULT_MIN_CHG_FRAC,
                        help="Effect-size guard: CHG frac must be >= this to call teM")
    parser.add_argument("--max-cg-frac-um", type=float, default=DEFAULT_MAX_CG_FRAC_UM,
                        help="Effect-size guard: CG frac must be <= this to call UM. "
                             "Genes that are non-significant but have intermediate "
                             "CG_frac become Unclassified rather than UM. Important "
                             "in genomes with high background pCG, where the one-sided "
                             "binomial test cannot distinguish unmethylated from "
                             "intermediate genes.")

    args = parser.parse_args()

    # ---- Load ----
    genes = load_genes_bed(args.genes_bed)

    print("\nLoading methylation data...")
    cg_pos,  cg_mr,  cg_cv  = load_methylation(args.cg_file,  args.min_cov, args.format)
    chg_pos, chg_mr, chg_cv = load_methylation(args.chg_file, args.min_cov, args.format)

    # ---- Aggregate ----
    print("\nAggregating per gene...")
    cg_agg  = aggregate_all(genes, cg_pos,  cg_mr,  cg_cv)
    chg_agg = aggregate_all(genes, chg_pos, chg_mr, chg_cv)

    # ---- Background ----
    p_cg  = compute_background(cg_agg)
    p_chg = compute_background(chg_agg)
    print(f"\nBackgrounds (CDS site-averaged):")
    print(f"  pCG  = {p_cg:.5f}")
    print(f"  pCHG = {p_chg:.5f}")

    # ---- Binomial + BH per context ----
    print("\nRunning binomial tests + BH correction...")
    cg_res  = run_binomial_with_bh(cg_agg,  p_cg,  args.min_n_cg)
    chg_res = run_binomial_with_bh(chg_agg, p_chg, args.min_n_chg)

    # ---- Classify ----
    print("Classifying...")
    df = classify(cg_res, chg_res,
                  args.alpha,
                  args.min_n_cg, args.min_n_chg,
                  args.min_chg_frac,
                  args.max_cg_frac_um)

    out_tsv = f"{args.output_prefix}_classification.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"\nSaved: {out_tsv}")

    # ---- Summary ----
    print("\n--- Classification summary ---")
    print(df["Classification"].value_counts().to_string())

    # ---- Plots ----
    print("\nGenerating plots...")
    plot_distribution(df, "CHG_frac", "CHG",
                      args.min_chg_frac,
                      f"{args.output_prefix}_dist_CHG.png",
                      zoom_max=0.30)
    plot_cg_sanity(df, f"{args.output_prefix}_dist_CG.png",
                   threshold=args.max_cg_frac_um)
    plot_classification_bar(df, f"{args.output_prefix}_summary_bar.png")

    print("\nDone.")

if __name__ == "__main__":
    main()
