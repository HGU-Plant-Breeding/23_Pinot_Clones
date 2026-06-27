#!/usr/bin/env python3
# =============================================================
# 05_plot.py
# Generate all figures for haplotype masking validation
# Callipo et al. 2026
# =============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

WORK_DIR = "/"
AGG_DIR  = f"{WORK_DIR}/results/aggregated"
VAR_DIR  = f"{WORK_DIR}/results/variants"
OUT_DIR  = f"{WORK_DIR}/results/figures"

import os; os.makedirs(OUT_DIR, exist_ok=True)

# Colors
COL_MASKED   = "#003f5c"
COL_UNMASKED = "#f95d6a"
COL_ROH      = "#ffa600"
COL_BORDER   = "#8fbc8f"
COL_NONROH   = "grey70"

REGION_ORDER   = ["RoH", "border", "nonRoH"]
REGION_LABELS  = {"RoH": "RoH", "border": "5kb border", "nonRoH": "Non-RoH"}
REGION_COLORS  = {"RoH": COL_ROH, "border": COL_BORDER, "nonRoH": COL_NONROH}

# -------------------------------------------------------
# Load data
# -------------------------------------------------------
print("Loading data...")
cov_df      = pd.read_csv(f"{AGG_DIR}/coverage_all_clones.tsv", sep="\t")
cov_summary = pd.read_csv(f"{AGG_DIR}/coverage_summary.tsv",    sep="\t")
mapq_freq   = pd.read_csv(f"{AGG_DIR}/mapq_frequency.tsv",      sep="\t")
mapq_summ   = pd.read_csv(f"{AGG_DIR}/mapq_summary.tsv",        sep="\t")
snp_counts  = pd.read_csv(f"{VAR_DIR}/snp_counts.tsv",          sep="\t")
sv_counts   = pd.read_csv(f"{VAR_DIR}/sv_counts.tsv",           sep="\t")
quality_df  = pd.read_csv(f"{VAR_DIR}/snp_quality.tsv",         sep="\t")
quality_df  = quality_df.dropna()
quality_df["MQ"] = pd.to_numeric(quality_df["MQ"], errors="coerce")
quality_df["AF"] = pd.to_numeric(quality_df["AF"], errors="coerce")
quality_df = quality_df.dropna()

# -------------------------------------------------------
# FIGURE 1: Coverage and MAPQ (2x2 grid)
# Panel A: Mean depth per region (masked vs unmasked)
# Panel B: Mean MAPQ per region (masked vs unmasked)
# Panel C: MAPQ=0 fraction per region (masked vs unmasked)
# Panel D: MAPQ distribution in RoH (masked vs unmasked)
# -------------------------------------------------------
print("Plotting Figure 1: Coverage and MAPQ...")

fig1, axes = plt.subplots(2, 2, figsize=(12, 10))
fig1.subplots_adjust(hspace=0.38, wspace=0.35)

x       = np.arange(len(REGION_ORDER))
width   = 0.35

# --- Panel A: Mean depth ---
ax = axes[0, 0]
for i, (strategy, color, label) in enumerate([
        ("masked",   COL_MASKED,   "Masked"),
        ("unmasked", COL_UNMASKED, "Unmasked")]):
    vals = []
    errs = []
    for region in REGION_ORDER:
        sub = cov_df[(cov_df["strategy"] == strategy) &
                     (cov_df["region"] == region)]["mean_depth"]
        vals.append(sub.mean())
        errs.append(sub.std())
    bars = ax.bar(x + i * width - width/2, vals, width,
                  yerr=errs, capsize=3,
                  color=color, alpha=0.85,
                  label=label, edgecolor="white", linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER], fontsize=10)
ax.set_ylabel("Mean sequencing depth (x)", fontsize=11, fontweight="bold")
ax.set_title("A. Sequencing depth by region", fontsize=12, fontweight="bold", loc="left")
ax.legend(fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

# --- Panel B: Mean MAPQ per region ---
ax = axes[0, 1]
for i, (strategy, color, label) in enumerate([
        ("masked",   COL_MASKED,   "Masked"),
        ("unmasked", COL_UNMASKED, "Unmasked")]):
    vals = []
    errs = []
    for region in REGION_ORDER:
        sub = cov_df[(cov_df["strategy"] == strategy) &
                     (cov_df["region"] == region)]["mean_mapq"]
        vals.append(sub.mean())
        errs.append(sub.std())
    ax.bar(x + i * width - width/2, vals, width,
           yerr=errs, capsize=3,
           color=color, alpha=0.85,
           label=label, edgecolor="white", linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER], fontsize=10)
ax.set_ylabel("Mean MAPQ", fontsize=11, fontweight="bold")
ax.set_title("B. Mean MAPQ by region", fontsize=12, fontweight="bold", loc="left")
ax.legend(fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

# --- Panel C: MAPQ=0 fraction in RoH ---
ax = axes[1, 0]
for i, (strategy, color, label) in enumerate([
        ("masked",   COL_MASKED,   "Masked"),
        ("unmasked", COL_UNMASKED, "Unmasked")]):
    vals = []
    for region in REGION_ORDER:
        sub = mapq_summ[(mapq_summ["strategy"] == strategy) &
                        (mapq_summ["region"] == region)]["pct_mapq0"]
        vals.append(sub.values[0] if len(sub) > 0 else 0)
    ax.bar(x + i * width - width/2, vals, width,
           color=color, alpha=0.85,
           label=label, edgecolor="white", linewidth=0.5)

ax.set_xticks(x)
ax.set_xticklabels([REGION_LABELS[r] for r in REGION_ORDER], fontsize=10)
ax.set_ylabel("Reads with MAPQ=0 (%)", fontsize=11, fontweight="bold")
ax.set_title("C. MAPQ=0 reads by region", fontsize=12, fontweight="bold", loc="left")
ax.legend(fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

# --- Panel D: MAPQ distribution in RoH ---
ax = axes[1, 1]
for strategy, color, label in [
        ("masked",   COL_MASKED,   "Masked RoH"),
        ("unmasked", COL_UNMASKED, "Unmasked RoH")]:
    sub = mapq_freq[(mapq_freq["strategy"] == strategy) &
                    (mapq_freq["region"] == "RoH")]
    if len(sub) > 0:
        ax.plot(sub["mapq"], sub["fraction"] * 100,
                color=color, linewidth=1.5, label=label)
        ax.fill_between(sub["mapq"], sub["fraction"] * 100,
                        alpha=0.2, color=color)

ax.set_xlabel("MAPQ", fontsize=11, fontweight="bold")
ax.set_ylabel("Reads (%)", fontsize=11, fontweight="bold")
ax.set_title("D. MAPQ distribution in RoH", fontsize=12, fontweight="bold", loc="left")
ax.set_xlim(0, 60)
ax.legend(fontsize=9)
ax.spines[["top", "right"]].set_visible(False)

plt.savefig(f"{OUT_DIR}/Fig_coverage_mapq.pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{OUT_DIR}/Fig_coverage_mapq.png", bbox_inches="tight", dpi=300)
print(f"  -> Fig_coverage_mapq saved")
plt.close()

# -------------------------------------------------------
# FIGURE 2: Variant recovery (SNPs + SVs)
# Panel A: SNP counts in RoH — filtered
# Panel B: SNP counts in RoH — unfiltered
# Panel C: SV counts in RoH — filtered
# Panel D: SV counts in RoH — unfiltered
# -------------------------------------------------------
print("Plotting Figure 2: Variant recovery...")

fig2, axes = plt.subplots(1, 4, figsize=(16, 5))
fig2.subplots_adjust(wspace=0.4)

def plot_variant_counts(ax, df, var_type, filtered, title, panel_label):
    sub = df[df["filtered"] == filtered]
    categories  = ["RoH", "nonRoH"]
    cat_labels  = ["RoH", "Non-RoH"]
    x = np.arange(len(categories))

    for i, (strategy, color, label) in enumerate([
            ("masked",   COL_MASKED,   "Masked"),
            ("unmasked", COL_UNMASKED, "Unmasked")]):
        vals = []
        for cat in categories:
            row = sub[(sub["vcf"] == strategy) & (sub["region"] == cat)]
            vals.append(row["count"].values[0] if len(row) > 0 else 0)
        bars = ax.bar(x + i * width - width/2, vals, width,
                      color=color, alpha=0.85, label=label,
                      edgecolor="white", linewidth=0.5)
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + max(vals) * 0.01,
                    f"{val:,}", ha="center", va="bottom",
                    fontsize=7, rotation=45)

    ax.set_xticks(x)
    ax.set_xticklabels(cat_labels, fontsize=10)
    ax.set_ylabel(f"Number of {var_type}s", fontsize=10, fontweight="bold")
    ax.set_title(f"{panel_label}. {var_type} ({filtered})",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

plot_variant_counts(axes[0], snp_counts, "SNP", "filtered",   "SNPs filtered",   "A")
plot_variant_counts(axes[1], snp_counts, "SNP", "unfiltered", "SNPs unfiltered", "B")
plot_variant_counts(axes[2], sv_counts,  "SV",  "filtered",   "SVs filtered",    "C")
plot_variant_counts(axes[3], sv_counts,  "SV",  "unfiltered", "SVs unfiltered",  "D")

plt.savefig(f"{OUT_DIR}/Fig_variant_recovery.pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{OUT_DIR}/Fig_variant_recovery.png", bbox_inches="tight", dpi=300)
print(f"  -> Fig_variant_recovery saved")
plt.close()

# -------------------------------------------------------
# FIGURE 3: Variant quality validation
# Panel A: MQ distribution — masked RoH vs masked nonRoH
# Panel B: AF distribution — masked RoH vs masked nonRoH
# Panel C: MQ distribution — masked RoH vs unmasked RoH
# Panel D: AF distribution — masked RoH vs unmasked RoH
# -------------------------------------------------------
print("Plotting Figure 3: Variant quality validation...")

fig3, axes = plt.subplots(1, 4, figsize=(16, 5))
fig3.subplots_adjust(wspace=0.4)

def plot_quality_dist(ax, df, metric, groups, colors, labels, title, panel, xlabel):
    bins = np.linspace(df[metric].min(), df[metric].max(), 40)
    for group, color, label in zip(groups, colors, labels):
        sub = df[df["label"] == group][metric].dropna()
        ax.hist(sub, bins=bins, density=True, alpha=0.6,
                color=color, label=f"{label}\n(n={len(sub):,})",
                linewidth=0)
        ax.axvline(sub.median(), color=color, linestyle="--",
                   linewidth=1.2, alpha=0.8)
    ax.set_xlabel(xlabel, fontsize=10, fontweight="bold")
    ax.set_ylabel("Density", fontsize=10, fontweight="bold")
    ax.set_title(f"{panel}. {title}", fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

# MQ: masked RoH vs masked nonRoH
plot_quality_dist(
    axes[0], quality_df, "MQ",
    ["masked_RoH", "masked_nonRoH"],
    [COL_ROH, COL_NONROH],
    ["Masked RoH", "Masked non-RoH"],
    "MQ: RoH vs non-RoH\n(masked)", "A", "Mapping Quality (MQ)"
)

# AF: masked RoH vs masked nonRoH
plot_quality_dist(
    axes[1], quality_df, "AF",
    ["masked_RoH", "masked_nonRoH"],
    [COL_ROH, COL_NONROH],
    ["Masked RoH", "Masked non-RoH"],
    "AF: RoH vs non-RoH\n(masked)", "B", "Allele Frequency (AF)"
)

# MQ: masked RoH vs unmasked RoH
plot_quality_dist(
    axes[2], quality_df, "MQ",
    ["masked_RoH", "unmasked_RoH"],
    [COL_MASKED, COL_UNMASKED],
    ["Masked RoH", "Unmasked RoH"],
    "MQ: masked vs unmasked\n(RoH only)", "C", "Mapping Quality (MQ)"
)

# AF: masked RoH vs unmasked RoH
plot_quality_dist(
    axes[3], quality_df, "AF",
    ["masked_RoH", "unmasked_RoH"],
    [COL_MASKED, COL_UNMASKED],
    ["Masked RoH", "Unmasked RoH"],
    "AF: masked vs unmasked\n(RoH only)", "D", "Allele Frequency (AF)"
)

plt.savefig(f"{OUT_DIR}/Fig_variant_quality.pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{OUT_DIR}/Fig_variant_quality.png", bbox_inches="tight", dpi=300)
print(f"  -> Fig_variant_quality saved")
plt.close()

print("\nAll figures saved to:", OUT_DIR)
print("Done.")
