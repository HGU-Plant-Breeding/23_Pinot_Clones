#!/usr/bin/env python3
"""
Step 1: Explore genome-wide methylation distributions per context
Samples 2M sites per clone (post coverage filter >=4) and plots density curves
Usage: python 01_methylation_distributions.py
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import gaussian_kde

# ── CONFIG ────────────────────────────────────────────────────────────────────
CG_DIR  = "/"
CHG_DIR = "/"
CHH_DIR = "/"

MIN_COV    = 4
N_SAMPLE   = 2_000_000   # sites per clone after coverage filter
RANDOM_SEED = 42
OUT_FIG    = "methylation_distributions.pdf"

CONTEXTS = {
    "CG":  CG_DIR,
    "CHG": CHG_DIR,
    "CHH": CHH_DIR,
}

# colour palette — one per clone, cycling
PALETTE = plt.cm.tab20.colors  # 20 colours, cycles for 23 clones

# ── FUNCTIONS ─────────────────────────────────────────────────────────────────
def load_clone(filepath, min_cov=4, n_sample=None, seed=42):
    """Load a single clone file, apply coverage filter, optionally sample."""
    df = pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["chr", "pos", "coverage", "methylation"],
        dtype={"chr": str, "pos": int, "coverage": float, "methylation": float},
    )
    df = df[df["coverage"] >= min_cov]
    if n_sample and len(df) > n_sample:
        df = df.sample(n=n_sample, random_state=seed)
    return df["methylation"].values


def kde_curve(values, bw_method=0.05):
    """Compute KDE curve over 0-100 range."""
    x = np.linspace(0, 100, 500)
    kde = gaussian_kde(values, bw_method=bw_method)
    y = kde(x)
    return x, y


# ── MAIN ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=False)
fig.suptitle(
    "Genome-wide methylation level distributions per sequence context\n"
    f"Coverage ≥ {MIN_COV}x · {N_SAMPLE:,} sites sampled per clone",
    fontsize=11, y=1.01
)

for ax, (context, dirpath) in zip(axes, CONTEXTS.items()):
    files = sorted(glob.glob(os.path.join(dirpath, "*_merged.txt")))
    print(f"\n[{context}] Found {len(files)} files")

    all_values = []

    for i, fp in enumerate(files):
        clone_name = os.path.basename(fp).replace("_merged.txt", "")
        print(f"  Loading {clone_name}...", end=" ", flush=True)

        vals = load_clone(fp, min_cov=MIN_COV, n_sample=N_SAMPLE, seed=RANDOM_SEED)
        print(f"{len(vals):,} sites retained")

        color = PALETTE[i % len(PALETTE)]

        # individual clone KDE
        x, y = kde_curve(vals, bw_method=0.04)
        ax.plot(x, y, color=color, alpha=0.45, linewidth=0.9, label=clone_name)

        all_values.append(vals)

    # aggregate KDE across all clones
    all_vals_concat = np.concatenate(all_values)
    x_agg, y_agg = kde_curve(all_vals_concat, bw_method=0.04)
    ax.plot(x_agg, y_agg, color="black", linewidth=2.0,
            linestyle="--", label="All clones (aggregate)", zorder=5)

    ax.set_title(context, fontsize=13, fontweight="bold", pad=8)
    ax.set_xlabel("Methylation level (%)", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_xlim(0, 100)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(20))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(10))
    ax.tick_params(axis="both", labelsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.3, linewidth=0.5)

    # legend only on last panel to avoid clutter
    if context == "CHH":
        ax.legend(
            fontsize=6.5, ncol=2, loc="upper left",
            framealpha=0.7, edgecolor="none"
        )

plt.tight_layout()
fig.savefig(OUT_FIG, dpi=200, bbox_inches="tight")
print(f"\nFigure saved → {OUT_FIG}")
