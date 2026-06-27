import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob
import warnings
warnings.filterwarnings("ignore")

# =============================================================
# RoH mapping metrics plots
# Callipo et al. 2026
# =============================================================

WORK_DIR = ""
METRICS_DIR = f"{WORK_DIR}/results/metrics"
OUT_DIR     = f"{WORK_DIR}/results/figures"

import os; os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------------------------------------
# 1. Load data
# -------------------------------------------------------
files = glob.glob(f"{METRICS_DIR}/*_RoH_metrics.tsv")
df = pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)
print(f"Loaded {len(df)} rows, {df['clone'].nunique()} clones")

# -------------------------------------------------------
# 2. Strategy display settings
# Only show the 3 most informative strategies:
#   masked_HapA, unmasked_HapA, unmasked_both
# unmasked_HapB is symmetric to HapA so redundant
# -------------------------------------------------------
strategies = ["masked_HapA", "unmasked_HapA", "unmasked_both"]
labels     = {
    "masked_HapA"   : "Masked\n(HapA)",
    "unmasked_HapA" : "Unmasked\n(HapA)",
    "unmasked_both" : "Unmasked\n(HapA+HapB)"
}
colors = {
    "masked_HapA"   : "#003f5c",
    "unmasked_HapA" : "#f95d6a",
    "unmasked_both" : "#ffa600"
}

df_plot = df[df["strategy"].isin(strategies)].copy()

# -------------------------------------------------------
# 3. Plot
# -------------------------------------------------------
fig = plt.figure(figsize=(14, 5))
gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.4)

# --- Panel A: Primary MAPQ=0 % ---
ax1 = fig.add_subplot(gs[0])
data_a = [df_plot[df_plot["strategy"] == s]["primary_mapq_lt10_pct"].values
          for s in strategies]
bp1 = ax1.boxplot(data_a, patch_artist=True, widths=0.5,
                  medianprops=dict(color="white", linewidth=2),
                  whiskerprops=dict(linewidth=1.2),
                  capprops=dict(linewidth=1.2),
                  flierprops=dict(marker="o", markersize=3, alpha=0.5))
for patch, s in zip(bp1["boxes"], strategies):
    patch.set_facecolor(colors[s])
    patch.set_alpha(0.85)
for flier, s in zip(bp1["fliers"], strategies):
    flier.set_markerfacecolor(colors[s])

ax1.set_xticks(range(1, len(strategies) + 1))
ax1.set_xticklabels([labels[s] for s in strategies], fontsize=9)
ax1.set_ylabel("Primary reads with MAPQ<10 (%)", fontsize=10, fontweight="bold")
ax1.set_title("A. Low MAPQ (<10) primary reads\nin RoH regions",
              fontsize=11, fontweight="bold", loc="left")
ax1.spines[["top", "right"]].set_visible(False)
ax1.tick_params(labelsize=9)

# --- Panel B: Mean MAPQ of primary reads ---
ax2 = fig.add_subplot(gs[1])
data_b = [df_plot[df_plot["strategy"] == s]["primary_mean_mapq"].values
          for s in strategies]
bp2 = ax2.boxplot(data_b, patch_artist=True, widths=0.5,
                  medianprops=dict(color="white", linewidth=2),
                  whiskerprops=dict(linewidth=1.2),
                  capprops=dict(linewidth=1.2),
                  flierprops=dict(marker="o", markersize=3, alpha=0.5))
for patch, s in zip(bp2["boxes"], strategies):
    patch.set_facecolor(colors[s])
    patch.set_alpha(0.85)
for flier, s in zip(bp2["fliers"], strategies):
    flier.set_markerfacecolor(colors[s])

ax2.set_xticks(range(1, len(strategies) + 1))
ax2.set_xticklabels([labels[s] for s in strategies], fontsize=9)
ax2.set_ylabel("Mean MAPQ (primary reads)", fontsize=10, fontweight="bold")
ax2.set_title("B. Mean MAPQ of primary reads\nin RoH regions",
              fontsize=11, fontweight="bold", loc="left")
ax2.spines[["top", "right"]].set_visible(False)
ax2.tick_params(labelsize=9)

# --- Panel C: Primary vs Secondary read counts ---
ax3 = fig.add_subplot(gs[2])
x      = np.arange(len(strategies))
width  = 0.35

for i, (read_type, col, label) in enumerate([
        ("primary_reads",   "#2c7fb8", "Primary"),
        ("secondary_reads", "#d7301f", "Secondary")]):
    vals = [df_plot[df_plot["strategy"] == s][read_type].mean() / 1000
            for s in strategies]
    errs = [df_plot[df_plot["strategy"] == s][read_type].std() / 1000
            for s in strategies]
    ax3.bar(x + i * width - width/2, vals, width,
            yerr=errs, capsize=3,
            color=col, alpha=0.85, label=label,
            edgecolor="white", linewidth=0.5)

ax3.set_xticks(x)
ax3.set_xticklabels([labels[s] for s in strategies], fontsize=9)
ax3.set_ylabel("Number of reads (×1,000)", fontsize=10, fontweight="bold")
ax3.set_title("C. Primary vs secondary reads\nin RoH regions",
              fontsize=11, fontweight="bold", loc="left")
ax3.legend(fontsize=9, framealpha=0.85)
ax3.spines[["top", "right"]].set_visible(False)
ax3.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/Fig_RoH_mapping_metrics.pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{OUT_DIR}/Fig_RoH_mapping_metrics.png", bbox_inches="tight", dpi=300)
print(f"-> Saved Fig_RoH_mapping_metrics")
plt.close()
