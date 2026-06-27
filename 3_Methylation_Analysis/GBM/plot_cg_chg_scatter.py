#!/usr/bin/env python3
"""
Script Name: plot_cg_chg_scatter.py
Description: CG vs CHG methylation scatter plot with marginal histograms,
             colored by gene classification. Reproduces the standard
             gbM/teM/UM/ambiguous diagnostic figure from the output of GBM_classify_CG_CHG.py.

             Each point is one gene. Marginal histograms on top and right
             show the per-class distribution of CG and CHG methylation.

             Low-coverage genes (CG_nSites < min-n-cg OR CHG_nSites < min-n-chg)
             are excluded from the plot and the count is annotated in a corner.
             "Unclassified" genes that pass the coverage filter are relabeled
             as "ambiguous" - these are the genes with intermediate methylation
             that failed both the gbM significance test and the UM effect-size
             guard.

Inputs:
  - <clone>_classification.tsv  produced by GBM_classify_CG_CHG.py
                                 (must contain Classification, CG_frac, CHG_frac,
                                  CG_nSites, CHG_nSites)

Outputs:
  - <o>.png                 the scatter+marginals figure

Usage:
  python plot_cg_chg_scatter.py per_clone/S1-2_classification.tsv \\
      --output S1-2_scatter.png --title "S1-2 all genes" \\
      --min-n-cg 15 --min-n-chg 15
"""

import sys
import os
import argparse

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
except ImportError:
    sys.exit("Error: pip install pandas matplotlib numpy")


# Class colors and order. UM purple and ambiguous salmon match the
# Niederhuth-style reference figure.
DEFAULT_COLORS = {
    "UM":        "#B084CC",   # lavender
    "gbM":       "#7BA428",   # olive green
    "teM":       "#3FB8AF",   # teal cyan
    "ambiguous": "#E8836E",   # salmon
}
# Drawing order: UM behind, then ambiguous, gbM, teM on top
DRAW_ORDER = ["UM", "ambiguous", "gbM", "teM"]


# ---------------- Font setup ----------------

def setup_roboto():
    """
    Try to use Roboto. matplotlib needs the font to be available either through
    fontconfig (system install) or installed in a directory it scans. If Roboto
    is not found, fall back silently to the default sans-serif and warn on stderr.
    Returns the family name actually selected.
    """
    available = {f.name for f in font_manager.fontManager.ttflist}
    if "Roboto" in available:
        plt.rcParams["font.family"]     = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Roboto"] + plt.rcParams["font.sans-serif"]
        return "Roboto"
    # Try refreshing the font cache once in case Roboto was installed after
    # matplotlib was imported on a previous run
    try:
        font_manager._load_fontmanager(try_read_cache=False)
        available = {f.name for f in font_manager.fontManager.ttflist}
        if "Roboto" in available:
            plt.rcParams["font.family"]     = "sans-serif"
            plt.rcParams["font.sans-serif"] = ["Roboto"] + plt.rcParams["font.sans-serif"]
            return "Roboto"
    except Exception:
        pass
    sys.stderr.write(
        "Warning: Roboto font not found in matplotlib's font cache. "
        "Falling back to default sans-serif. To use Roboto, install it with:\n"
        "    pip install --user roboto-font\n"
        "or download the TTFs from https://fonts.google.com/specimen/Roboto, "
        "place them in ~/.fonts/ and run `fc-cache -f`.\n"
    )
    fam = plt.rcParams["font.family"]
    return fam[0] if isinstance(fam, list) else fam


# ---------------- Data ----------------

def load_classification(path, min_n_cg, min_n_chg):
    """
    Load classification TSV and split into:
      kept   : genes to plot
      dropped: low-coverage Unclassified genes (will not appear in plot)

    Coverage filtering is applied ONLY to Unclassified genes. gbM, teM and UM
    genes are always kept: their classification already implies they passed the
    relevant coverage thresholds inside the classifier (or in the case of teM,
    the call is based purely on CHG significance regardless of CG coverage).
    Filtering classified genes by coverage would produce count mismatches between
    the classifier summary and the plot.

    Within `kept`, "Unclassified" genes that survive the coverage filter are
    relabeled as "ambiguous" - these are genes with intermediate methylation
    that failed both the gbM significance test and the UM effect-size guard.
    """
    df = pd.read_csv(path, sep="\t")
    needed = ["Classification", "CG_frac", "CHG_frac", "CG_nSites", "CHG_nSites"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        sys.exit(f"Error: input TSV missing columns: {missing}")

    df["CG_frac"]    = pd.to_numeric(df["CG_frac"],    errors="coerce")
    df["CHG_frac"]   = pd.to_numeric(df["CHG_frac"],   errors="coerce")
    df["CG_nSites"]  = pd.to_numeric(df["CG_nSites"],  errors="coerce").fillna(0).astype(int)
    df["CHG_nSites"] = pd.to_numeric(df["CHG_nSites"], errors="coerce").fillna(0).astype(int)

    # Drop rows with no usable methylation values at all
    df = df.dropna(subset=["CG_frac", "CHG_frac"])

    # Coverage filter applies ONLY to Unclassified genes.
    # gbM, teM, UM are always plotted - their counts must match the classifier.
    is_unclassified = df["Classification"] == "Unclassified"
    cov_ok = (df["CG_nSites"] >= min_n_cg) & (df["CHG_nSites"] >= min_n_chg)

    drop_mask = is_unclassified & ~cov_ok
    kept    = df[~drop_mask].copy()
    dropped = df[drop_mask].copy()

    # Unclassified genes that survived the filter are the genuinely ambiguous
    # ones (intermediate CG, failed UM guard). Relabel for the plot.
    kept["Classification"] = kept["Classification"].replace(
        {"Unclassified": "ambiguous"}
    )

    return kept, dropped


# ---------------- Plot ----------------

def make_plot(df_kept, n_dropped, title, outpath, point_size, alpha,
              colors, max_points_per_class):
    # Per-class data
    data = {cl: df_kept[df_kept["Classification"] == cl] for cl in DRAW_ORDER}
    counts = {cl: len(sub) for cl, sub in data.items()}

    # ---- Layout: main scatter + top hist + right hist ----
    fig = plt.figure(figsize=(8.5, 8.5))
    gs = fig.add_gridspec(2, 2,
                          width_ratios=(5, 1),
                          height_ratios=(1, 5),
                          left=0.10, right=0.97,
                          bottom=0.08, top=0.93,
                          wspace=0.03, hspace=0.03)
    ax_main  = fig.add_subplot(gs[1, 0])
    ax_top   = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    # ---- Main scatter ----
    for cl in DRAW_ORDER:
        sub = data[cl]
        if len(sub) == 0:
            continue
        # Subsample dense classes for plot clarity. Legend always shows
        # the true count.
        if max_points_per_class and len(sub) > max_points_per_class:
            sub = sub.sample(n=max_points_per_class, random_state=0)
        ax_main.scatter(sub["CG_frac"], sub["CHG_frac"],
                        s=point_size, alpha=alpha,
                        c=colors[cl], edgecolors="none",
                        rasterized=True,
                        label=f"{cl} ({counts[cl]:,})")

    ax_main.set_xlim(0, 1)
    ax_main.set_ylim(0, 1)
    ax_main.set_xlabel("CG methylation", fontsize=13)
    ax_main.set_ylabel("CHG methylation", fontsize=13)
    leg = ax_main.legend(title="Gene methylation state", loc="upper left",
                     	 bbox_to_anchor=(0.01, 0.94),
                         fontsize=11, title_fontsize=12,
                         frameon=False, markerscale=3)
    handles = leg.legend_handles if hasattr(leg, "legend_handles") else leg.legendHandles
    for h in handles:
        h.set_alpha(1.0)

    # ---- Annotation: dropped low-coverage genes ----
    if n_dropped > 0:
        ax_main.text(
            0.015, 0.985,
            f"Genes excluded: {n_dropped:,}",
            transform=ax_main.transAxes,
            ha="left", va="top",
            fontsize=9, color="#555555", style="italic",
        )

    # ---- Top marginal: stacked CG histogram per class ----
    bins = np.linspace(0, 1, 60)
    top_data   = [data[cl]["CG_frac"].values for cl in DRAW_ORDER if len(data[cl]) > 0]
    top_colors = [colors[cl]                 for cl in DRAW_ORDER if len(data[cl]) > 0]
    ax_top.hist(top_data, bins=bins, stacked=True,
                color=top_colors, edgecolor="white", linewidth=0.2)
    ax_top.set_yticks([])
    ax_top.tick_params(axis="x", labelbottom=False)
    for spine in ("top", "right", "left"):
        ax_top.spines[spine].set_visible(False)

    # ---- Right marginal: stacked CHG histogram, horizontal ----
    right_data   = [data[cl]["CHG_frac"].values for cl in DRAW_ORDER if len(data[cl]) > 0]
    right_colors = [colors[cl]                  for cl in DRAW_ORDER if len(data[cl]) > 0]
    ax_right.hist(right_data, bins=bins, stacked=True,
                  orientation="horizontal",
                  color=right_colors, edgecolor="white", linewidth=0.2)
    ax_right.set_xticks([])
    ax_right.tick_params(axis="y", labelleft=False)
    for spine in ("top", "right", "bottom"):
        ax_right.spines[spine].set_visible(False)

    # ---- Title ----
    if title:
        ax_top.set_title(title, fontsize=14, fontweight="bold", pad=10)

    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath}")
    print(f"  Plotted: " + "  ".join(f"{cl}={counts[cl]:,}" for cl in DRAW_ORDER))
    print(f"  Excluded (low coverage): {n_dropped:,}")


def main():
    ap = argparse.ArgumentParser(
        description="CG vs CHG scatter + marginal histograms, colored by classification.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("classification_tsv",
                    help="Output TSV from GBM_classify_CG_CHG.py")
    ap.add_argument("--output",      default=None,
                    help="Output PNG path. Default: <input_basename>_scatter.png")
    ap.add_argument("--title",       default=None,
                    help="Figure title (e.g. clone name). Default: input basename.")
    ap.add_argument("--min-n-cg",    type=int, default=15,
                    help="Minimum CG_nSites to keep a gene in the plot. "
                         "Should match the value used in the classifier run.")
    ap.add_argument("--min-n-chg",   type=int, default=15,
                    help="Minimum CHG_nSites to keep a gene in the plot. "
                         "Should match the value used in the classifier run.")
    ap.add_argument("--point-size",  type=float, default=3.0,
                    help="Scatter point size.")
    ap.add_argument("--alpha",       type=float, default=0.55,
                    help="Scatter point alpha (transparency).")
    ap.add_argument("--max-points-per-class", type=int, default=20000,
                    help="Subsample classes larger than this for plot clarity. "
                         "0 = no subsampling. Legend always shows the true count.")
    args = ap.parse_args()

    if not os.path.isfile(args.classification_tsv):
        sys.exit(f"Error: input TSV not found: {args.classification_tsv}")

    fam = setup_roboto()
    print(f"Font family: {fam}")

    kept, dropped = load_classification(args.classification_tsv,
                                         args.min_n_cg, args.min_n_chg)
    print(f"Loaded {len(kept) + len(dropped):,} genes from "
          f"{os.path.basename(args.classification_tsv)}")
    print(f"  Kept (CG_nSites>={args.min_n_cg}, CHG_nSites>={args.min_n_chg}): "
          f"{len(kept):,}")
    print(f"  Dropped (low coverage): {len(dropped):,}")

    out = args.output
    if out is None:
        base = os.path.splitext(os.path.basename(args.classification_tsv))[0]
        out  = f"{base}_scatter.png"

    title = args.title
    if title is None:
        base = os.path.splitext(os.path.basename(args.classification_tsv))[0]
        if base.endswith("_classification"):
            base = base[:-len("_classification")]
        title = f"{base} all genes"

    make_plot(kept, len(dropped), title, out,
              args.point_size, args.alpha,
              DEFAULT_COLORS,
              args.max_points_per_class if args.max_points_per_class > 0 else None)


if __name__ == "__main__":
    main()
