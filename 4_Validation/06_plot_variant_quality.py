import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

# =============================================================
# Variant quality validation
# Masked RoH vs masked nonRoH (genome-wide baseline)
# Callipo et al. 2026
# =============================================================

WORK_DIR = "/"
OUT_DIR  = f"{WORK_DIR}/results/figures"

import os; os.makedirs(OUT_DIR, exist_ok=True)

# Load
snp = pd.read_csv(f"{WORK_DIR}/results/variant_metrics/snp_metrics.tsv", sep="\t")
snp["MQ"]   = pd.to_numeric(snp["MQ"],   errors="coerce")
snp["DP"]   = pd.to_numeric(snp["DP"],   errors="coerce")
snp["MQ0F"] = pd.to_numeric(snp["MQ0F"], errors="coerce")
snp["AF"]   = pd.to_numeric(snp["AF"],   errors="coerce")

# Focus: masked only
roh    = snp[(snp["strategy"] == "masked") & (snp["region"] == "RoH")].dropna()
nonroh = snp[(snp["strategy"] == "masked") & (snp["region"] == "nonRoH")].dropna()

print(f"Masked RoH SNPs:    {len(roh):,}")
print(f"Masked nonRoH SNPs: {len(nonroh):,}")

# Subsample nonRoH for plotting (too many points otherwise)
np.random.seed(42)
nonroh_plot = nonroh.sample(min(10000, len(nonroh)))

COL_ROH    = "#ffa600"
COL_NONROH = "#003f5c"

# -------------------------------------------------------
# Plot
# -------------------------------------------------------
fig = plt.figure(figsize=(8, 5))
gs  = gridspec.GridSpec(1, 2, figure=fig, wspace=0.4)

# -------------------------------------------------------
# Panel A: MQ distribution
# -------------------------------------------------------
ax1 = fig.add_subplot(gs[0])
bins_mq = np.arange(0, 62, 2)
ax1.hist(nonroh_plot["MQ"], bins=bins_mq, density=True, alpha=0.7,
         color=COL_NONROH, label=f"Non-RoH (n={len(nonroh):,})", linewidth=0)
ax1.hist(roh["MQ"], bins=bins_mq, density=True, alpha=0.7,
         color=COL_ROH, label=f"RoH (n={len(roh):,})", linewidth=0)
ax1.axvline(x=20, color="red", linestyle="--", linewidth=1,
            label="MQ filter (20)", alpha=0.8)
ax1.set_xlabel("Mapping Quality (MQ)", fontsize=10, fontweight="bold")
ax1.set_ylabel("Density", fontsize=10, fontweight="bold")
ax1.set_title("A. MQ distribution", fontsize=11, fontweight="bold", loc="left")
ax1.legend(fontsize=8, framealpha=0.85)
ax1.spines[["top", "right"]].set_visible(False)
ax1.text(0.97, 0.97,
         f"RoH median: {roh['MQ'].median():.1f}\nGenome median: {nonroh['MQ'].median():.1f}",
         transform=ax1.transAxes, fontsize=8,
         ha="right", va="top",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   alpha=0.85, edgecolor="grey"))

# -------------------------------------------------------
# Panel B: AF spectrum
# -------------------------------------------------------
ax2 = fig.add_subplot(gs[1])
bins_af = np.linspace(0, 1, 30)
ax2.hist(nonroh_plot["AF"], bins=bins_af, density=True, alpha=0.7,
         color=COL_NONROH, label="Non-RoH", linewidth=0)
ax2.hist(roh["AF"], bins=bins_af, density=True, alpha=0.7,
         color=COL_ROH, label="RoH", linewidth=0)
ax2.set_xlabel("Allele Frequency (AF)", fontsize=10, fontweight="bold")
ax2.set_ylabel("Density", fontsize=10, fontweight="bold")
ax2.set_title("B. AF spectrum", fontsize=11, fontweight="bold", loc="left")
ax2.legend(fontsize=8, framealpha=0.85)
ax2.spines[["top", "right"]].set_visible(False)
ax2.text(0.97, 0.97,
         f"RoH median: {roh['AF'].median():.3f}\nGenome median: {nonroh['AF'].median():.3f}",
         transform=ax2.transAxes, fontsize=8,
         ha="right", va="top",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                   alpha=0.85, edgecolor="grey"))

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/Fig_variant_quality.pdf", bbox_inches="tight", dpi=300)
plt.savefig(f"{OUT_DIR}/Fig_variant_quality.png", bbox_inches="tight", dpi=300)
print(f"-> Saved Fig_variant_quality")
plt.close()
