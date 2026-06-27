import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

# =======================================================
# HiFi vs ONT CG Methylation Comparison
# Reference clone: 20-13 Gm
# Callipo et al. 2026
# =======================================================

WORKING_DIR = "/"
HIFI_FILE   = WORKING_DIR + "PacBio.combined.bed"
ONT_FILE    = WORKING_DIR + "ONT_Methylation_20-13.tsv"
OUT_PREFIX  = WORKING_DIR + "HiFi_vs_ONT"

MIN_COV       = 4
BIN_SIZE      = 200
METH_THRESH   = 70
UNMETH_THRESH = 30

# -------------------------------------------------------
# 1. Load files
# -------------------------------------------------------
print("Loading HiFi data...")
hifi = pd.read_csv(HIFI_FILE, sep="\t", comment="#", header=None,
                   usecols=[0, 1, 2, 3, 5],
                   names=["chrom", "start", "end", "mod_score", "cov"],
                   compression="infer")
hifi = hifi[hifi["cov"] >= MIN_COV].copy()
hifi["meth_hifi"] = hifi["mod_score"]
print(f"  HiFi sites after coverage filter: {len(hifi):,}")

print("Loading ONT data...")
ont = pd.read_csv(ONT_FILE, sep="\t", comment="#", header=None,
                  usecols=[0, 1, 2, 3, 4],
                  names=["chrom", "start", "end", "cov", "meth_ont"],
                  compression="infer")
ont = ont[ont["cov"] >= MIN_COV].copy()
print(f"  ONT sites after coverage filter: {len(ont):,}")

# -------------------------------------------------------
# 2. Join on chrom + start + end
# -------------------------------------------------------
print("Joining datasets...")
merged = pd.merge(
    hifi[["chrom", "start", "end", "meth_hifi"]],
    ont[["chrom",  "start", "end", "meth_ont"]],
    on=["chrom", "start", "end"],
    how="inner"
)
print(f"  Sites in common: {len(merged):,}")

# -------------------------------------------------------
# 3. Per-site correlation
# -------------------------------------------------------
r_pearson,  p_pearson  = stats.pearsonr(merged["meth_hifi"], merged["meth_ont"])
r_spearman, p_spearman = stats.spearmanr(merged["meth_hifi"], merged["meth_ont"])

print(f"\nPer-site correlation:")
print(f"  Pearson  r = {r_pearson:.4f}  (p = {p_pearson:.2e})")
print(f"  Spearman r = {r_spearman:.4f}  (p = {p_spearman:.2e})")

# -------------------------------------------------------
# 4. Misclassification — per site (70/30)
# -------------------------------------------------------
def binarize(series):
    b = pd.Series(np.nan, index=series.index)
    b[series >= METH_THRESH]   = 1
    b[series <= UNMETH_THRESH] = 0
    return b

merged["bin_hifi"] = binarize(merged["meth_hifi"])
merged["bin_ont"]  = binarize(merged["meth_ont"])

both_called = merged.dropna(subset=["bin_hifi", "bin_ont"])
n_both = len(both_called)

concordant_meth   = ((both_called["bin_ont"] == 1) & (both_called["bin_hifi"] == 1)).sum()
concordant_unmeth = ((both_called["bin_ont"] == 0) & (both_called["bin_hifi"] == 0)).sum()
ont1_hifi0        = ((both_called["bin_ont"] == 1) & (both_called["bin_hifi"] == 0)).sum()
ont0_hifi1        = ((both_called["bin_ont"] == 0) & (both_called["bin_hifi"] == 1)).sum()

print(f"\nMisclassification analysis (70/30 thresholds):")
print(f"  Sites with unambiguous call in both: {n_both:,}")
print(f"  Concordant methylated   (ONT=1, HiFi=1): {concordant_meth:,}  ({concordant_meth/n_both*100:.2f}%)")
print(f"  Concordant unmethylated (ONT=0, HiFi=0): {concordant_unmeth:,}  ({concordant_unmeth/n_both*100:.2f}%)")
print(f"  ONT overcalls  (ONT=1, HiFi=0): {ont1_hifi0:,}  ({ont1_hifi0/n_both*100:.2f}%)")
print(f"  ONT undercalls (ONT=0, HiFi=1): {ont0_hifi1:,}  ({ont0_hifi1/n_both*100:.2f}%)")
print(f"  Total discordant: {ont1_hifi0+ont0_hifi1:,}  ({(ont1_hifi0+ont0_hifi1)/n_both*100:.2f}%)")

# Save discordant sites list
discordant_sites = both_called[
    ((both_called["bin_ont"] == 1) & (both_called["bin_hifi"] == 0)) |
    ((both_called["bin_ont"] == 0) & (both_called["bin_hifi"] == 1))
].copy()
discordant_sites["discordance_type"] = np.where(
    (discordant_sites["bin_ont"] == 1) & (discordant_sites["bin_hifi"] == 0),
    "ONT_overcall", "ONT_undercall"
)
discordant_sites[["chrom", "start", "end", "meth_hifi", "meth_ont",
                   "bin_hifi", "bin_ont", "discordance_type"]].to_csv(
    OUT_PREFIX + "_discordant_sites.tsv", sep="\t", index=False
)
print(f"  -> Discordant sites saved: {OUT_PREFIX}_discordant_sites.tsv")

# -------------------------------------------------------
# 5. Bin into 200bp windows
# -------------------------------------------------------
print(f"\nBinning into {BIN_SIZE}bp windows...")
merged["bin_id"]    = merged["chrom"].astype(str) + "_" + (merged["start"] // BIN_SIZE).astype(str)
merged["bin_start"] = (merged["start"] // BIN_SIZE) * BIN_SIZE
merged["bin_end"]   = merged["bin_start"] + BIN_SIZE

binned = merged.groupby(["chrom", "bin_start", "bin_end", "bin_id"]).agg(
    meth_hifi=("meth_hifi", "mean"),
    meth_ont =("meth_ont",  "mean"),
    n_sites  =("meth_hifi", "count")
).reset_index()

binned = binned[binned["n_sites"] >= 1]
print(f"  Total 200bp bins: {len(binned):,}")

r_bin_p, p_bin_p = stats.pearsonr(binned["meth_hifi"], binned["meth_ont"])
r_bin_s, p_bin_s = stats.spearmanr(binned["meth_hifi"], binned["meth_ont"])
print(f"\n200bp binned correlation:")
print(f"  Pearson  r = {r_bin_p:.4f}  (p = {p_bin_p:.2e})")
print(f"  Spearman r = {r_bin_s:.4f}  (p = {p_bin_s:.2e})")

# Misclassification at bin level
binned["bin_hifi_b"] = binarize(binned["meth_hifi"])
binned["bin_ont_b"]  = binarize(binned["meth_ont"])

both_bins = binned.dropna(subset=["bin_hifi_b", "bin_ont_b"])
n_both_bins = len(both_bins)

b_conc_meth   = ((both_bins["bin_ont_b"] == 1) & (both_bins["bin_hifi_b"] == 1)).sum()
b_conc_unmeth = ((both_bins["bin_ont_b"] == 0) & (both_bins["bin_hifi_b"] == 0)).sum()
b_ont1_hifi0  = ((both_bins["bin_ont_b"] == 1) & (both_bins["bin_hifi_b"] == 0)).sum()
b_ont0_hifi1  = ((both_bins["bin_ont_b"] == 0) & (both_bins["bin_hifi_b"] == 1)).sum()

print(f"\nBin-level misclassification (70/30 thresholds):")
print(f"  Bins with unambiguous call in both: {n_both_bins:,}")
print(f"  Concordant methylated   (ONT=1, HiFi=1): {b_conc_meth:,}  ({b_conc_meth/n_both_bins*100:.2f}%)")
print(f"  Concordant unmethylated (ONT=0, HiFi=0): {b_conc_unmeth:,}  ({b_conc_unmeth/n_both_bins*100:.2f}%)")
print(f"  ONT overcalls  (ONT=1, HiFi=0): {b_ont1_hifi0:,}  ({b_ont1_hifi0/n_both_bins*100:.2f}%)")
print(f"  ONT undercalls (ONT=0, HiFi=1): {b_ont0_hifi1:,}  ({b_ont0_hifi1/n_both_bins*100:.2f}%)")
print(f"  Total discordant bins: {b_ont1_hifi0+b_ont0_hifi1:,}  ({(b_ont1_hifi0+b_ont0_hifi1)/n_both_bins*100:.2f}%)")

# Save discordant bins list
discordant_bins = both_bins[
    ((both_bins["bin_ont_b"] == 1) & (both_bins["bin_hifi_b"] == 0)) |
    ((both_bins["bin_ont_b"] == 0) & (both_bins["bin_hifi_b"] == 1))
].copy()
discordant_bins["discordance_type"] = np.where(
    (discordant_bins["bin_ont_b"] == 1) & (discordant_bins["bin_hifi_b"] == 0),
    "ONT_overcall", "ONT_undercall"
)
discordant_bins[["chrom", "bin_start", "bin_end",
                  "meth_hifi", "meth_ont", "n_sites",
                  "bin_hifi_b", "bin_ont_b",
                  "discordance_type"]].to_csv(
    OUT_PREFIX + "_discordant_bins.tsv", sep="\t", index=False
)
print(f"  -> Discordant bins saved: {OUT_PREFIX}_discordant_bins.tsv")

# -------------------------------------------------------
# 6. Full misclassification summary to text
# -------------------------------------------------------
with open(OUT_PREFIX + "_misclassification_summary.txt", "w") as f:
    f.write("HiFi vs ONT CG Methylation — Misclassification Summary\n")
    f.write("=" * 60 + "\n")
    f.write(f"Thresholds: methylated > {METH_THRESH}%, unmethylated < {UNMETH_THRESH}%\n\n")
    f.write("--- Per-site ---\n")
    f.write(f"Sites in common: {len(merged):,}\n")
    f.write(f"Sites with unambiguous call in both: {n_both:,}\n")
    f.write(f"Concordant methylated   (ONT=1, HiFi=1): {concordant_meth:,} ({concordant_meth/n_both*100:.2f}%)\n")
    f.write(f"Concordant unmethylated (ONT=0, HiFi=0): {concordant_unmeth:,} ({concordant_unmeth/n_both*100:.2f}%)\n")
    f.write(f"ONT overcalls  (ONT=1, HiFi=0):          {ont1_hifi0:,} ({ont1_hifi0/n_both*100:.2f}%)\n")
    f.write(f"ONT undercalls (ONT=0, HiFi=1):          {ont0_hifi1:,} ({ont0_hifi1/n_both*100:.2f}%)\n")
    f.write(f"Total discordant:                         {ont1_hifi0+ont0_hifi1:,} ({(ont1_hifi0+ont0_hifi1)/n_both*100:.2f}%)\n\n")
    f.write("--- 200bp bins ---\n")
    f.write(f"Total bins: {len(binned):,}\n")
    f.write(f"Bins with unambiguous call in both: {n_both_bins:,}\n")
    f.write(f"Concordant methylated   (ONT=1, HiFi=1): {b_conc_meth:,} ({b_conc_meth/n_both_bins*100:.2f}%)\n")
    f.write(f"Concordant unmethylated (ONT=0, HiFi=0): {b_conc_unmeth:,} ({b_conc_unmeth/n_both_bins*100:.2f}%)\n")
    f.write(f"ONT overcalls  (ONT=1, HiFi=0):          {b_ont1_hifi0:,} ({b_ont1_hifi0/n_both_bins*100:.2f}%)\n")
    f.write(f"ONT undercalls (ONT=0, HiFi=1):          {b_ont0_hifi1:,} ({b_ont0_hifi1/n_both_bins*100:.2f}%)\n")
    f.write(f"Total discordant bins:                    {b_ont1_hifi0+b_ont0_hifi1:,} ({(b_ont1_hifi0+b_ont0_hifi1)/n_both_bins*100:.2f}%)\n")
print(f"  -> Summary saved: {OUT_PREFIX}_misclassification_summary.txt")

# -------------------------------------------------------
# 7. Plot
# -------------------------------------------------------
fig = plt.figure(figsize=(18, 5))
gs_main = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38)

# --- Panel A: per-site hexbin ---
ax1 = fig.add_subplot(gs_main[0])
hb1 = ax1.hexbin(
    merged["meth_hifi"], merged["meth_ont"],
    gridsize=80, cmap="YlOrRd", mincnt=1,
    bins="log", linewidths=0.1
)
ax1.plot([0, 100], [0, 100], color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
cb1 = fig.colorbar(hb1, ax=ax1, pad=0.02)
cb1.set_label("log$_{10}$(CpG sites)", fontsize=9, labelpad=8)
ax1.set_xlabel("HiFi methylation (%)", fontsize=11, fontweight="bold")
ax1.set_ylabel("ONT methylation (%)",  fontsize=11, fontweight="bold")
ax1.set_title("Per-site", fontsize=12, fontweight="bold")
ax1.text(0.05, 0.93,
         f"Pearson r = {r_pearson:.3f}\nSpearman r = {r_spearman:.3f}\nn = {len(merged):,}",
         transform=ax1.transAxes, fontsize=9, verticalalignment="top",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="grey"))
ax1.set_xlim(0, 100); ax1.set_ylim(0, 100)
ax1.set_aspect("equal"); ax1.tick_params(labelsize=9)

# --- Panel B: 200bp binned hexbin ---
ax2 = fig.add_subplot(gs_main[1])
hb2 = ax2.hexbin(
    binned["meth_hifi"], binned["meth_ont"],
    gridsize=80, cmap="YlOrRd", mincnt=1,
    bins="log", linewidths=0.1
)
ax2.plot([0, 100], [0, 100], color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
cb2 = fig.colorbar(hb2, ax=ax2, pad=0.02)
cb2.set_label("log$_{10}$(200bp bins)", fontsize=9, labelpad=8)
ax2.set_xlabel("HiFi methylation (%)", fontsize=11, fontweight="bold")
ax2.set_ylabel("ONT methylation (%)",  fontsize=11, fontweight="bold")
ax2.set_title("200bp binned", fontsize=12, fontweight="bold")
ax2.text(0.05, 0.93,
         f"Pearson r = {r_bin_p:.3f}\nSpearman r = {r_bin_s:.3f}\nn = {len(binned):,} bins",
         transform=ax2.transAxes, fontsize=9, verticalalignment="top",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="grey"))
ax2.set_xlim(0, 100); ax2.set_ylim(0, 100)
ax2.set_aspect("equal"); ax2.tick_params(labelsize=9)

# --- Panel C: distribution comparison ---
ax3 = fig.add_subplot(gs_main[2])
bins_hist = np.arange(0, 101, 2)
ax3.hist(merged["meth_hifi"], bins=bins_hist, density=True,
         alpha=0.6, color="#2c7fb8", label="HiFi", linewidth=0)
ax3.hist(merged["meth_ont"],  bins=bins_hist, density=True,
         alpha=0.6, color="#d7301f", label="ONT",  linewidth=0)
ax3.set_xlabel("CG methylation (%)", fontsize=11, fontweight="bold")
ax3.set_ylabel("Density",            fontsize=11, fontweight="bold")
ax3.set_title("Methylation distribution", fontsize=12, fontweight="bold")
ax3.legend(fontsize=10, framealpha=0.85)
ax3.set_xlim(0, 100)
ax3.tick_params(labelsize=9)

plt.tight_layout()
plt.savefig(OUT_PREFIX + ".pdf", bbox_inches="tight", dpi=300)
plt.savefig(OUT_PREFIX + ".png", bbox_inches="tight", dpi=300)
print(f"\n-> Figure saved: {OUT_PREFIX}.pdf / .png")
print("\nDone.")
