import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# =============================================================
# SNP and SV count comparison: masked vs unmasked in RoH
# Callipo et al. 2026
# =============================================================

WORK_DIR = "/"
OUT_DIR  = f"{WORK_DIR}/results/figures"

import os; os.makedirs(OUT_DIR, exist_ok=True)

COL_MASKED   = "#003f5c"
COL_UNMASKED = "#f95d6a"

# -------------------------------------------------------
# Data — from 04_variant_analysis.sh output
# -------------------------------------------------------
snp_data = {
    "filtered":   {"masked": 2704,  "unmasked": 494},
    "unfiltered": {"masked": 28415, "unmasked": 8814}
}
sv_data = {
    "filtered":   {"masked": 133,  "unmasked": 65},
    "unfiltered": {"masked": 679,  "unmasked": 276}
}

# -------------------------------------------------------
# Plot
# -------------------------------------------------------
fig = plt.figure(figsize=(10, 5))
gs  = gridspec.GridSpec(1, 2, figure=fig, wspace=0.4)

def plot_counts(ax, data, title, panel_label, ylabel):
    categories = ["filtered", "unfiltered"]
    labels     = ["Filtered", "Unfiltered"]
    x     = np.arange(len(categories))
    width = 0.35

    for i, (strategy, color, label) in enumerate([
            ("masked",   COL_MASKED,   "Masked (HapA)"),
            ("unmasked", COL_UNMASKED, "Unmasked (HapA+HapB)")]):
        vals = [data[cat][strategy] for cat in categories]
        bars = ax.bar(x + i * width - width/2, vals, width,
                      color=color, alpha=0.85, label=label,
                      edgecolor="white", linewidth=0.5)
        # Value labels on bars
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + max([data[c][s]
                        for c in categories
                        for s in ["masked","unmasked"]]) * 0.01,
                    f"{val:,}", ha="center", va="bottom",
                    fontsize=8, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10, fontweight="bold")
    ax.set_title(f"{panel_label}. {title}", fontsize=11,
                 fontweight="bold", loc="left")
    ax.legend(fontsize=9, framealpha=0.85)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=9)

plot_counts(fig.add_subplot(gs[0]), snp_data,
            "SNPs in RoH regions", "A", "Number of SNPs")
plot_counts(fig.add_subplot(gs[1]), sv_data,
            "SVs in RoH regions",  "B", "Number of SVs")

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/Fig_variant_counts_RoH.pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{OUT_DIR}/Fig_variant_counts_RoH.png", bbox_inches="tight", dpi=300)
print("-> Saved Fig_variant_counts_RoH")
plt.close()
