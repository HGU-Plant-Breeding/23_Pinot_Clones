#!/usr/bin/env python3
"""
Script Name: GBM_classify_multiclone.py
Description: Run the per-clone gbM/teM/UM classifier (GBM_classify_CG_CHG.py)
             across many clones from a sample sheet, then merge the per-clone
             results into a single wide table for cross-clone comparison and
             differential-class analysis.

             Uses CG and CHG only (CHH dropped due to noise in ONT calls).

             Each clone is classified independently with its own per-clone
             background (pCG, pCHG), since global methylation drift between
             clones would bias a shared background. The merged output flags
             genes whose classification differs between clones (singleton
             changes, splits, instability) for downstream differential
             analysis.

Inputs:
  - genes BED                    : same as the per-clone classifier expects
  - sample sheet TSV (tab-separated):
        clone_name  cg_file  chg_file  [chrom_rename]
        # comments allowed (lines starting with '#')
        chrom_rename is optional; if present, it's a comma-separated list of
        FROM=TO substring rewrites applied to the chromosome column of each
        methylation BED (longest FROM is matched first). Use this when the
        methylation BEDs come from a RagTag-style assembly with suffixes
        like '*_2_RagTag' that need to be remapped to the genes BED's
        haplotype convention (e.g. '*_HapB').

        Backward compatibility: if a 5-column sheet is supplied (legacy
        format with a CHH column), the CHH file is silently ignored. The
        chrom_rename, if present, is read from column 5.

  - all flags supported by GBM_classify_CG_CHG.py are forwarded to each
    per-clone run (e.g. --format counts, --min-cov, --min-n-cg, etc.)

Outputs (under <output_dir>):
  - per_clone/<clone>_classification.tsv     each clone's full TSV
  - multiclone_classification_wide.tsv       gene_id + one class column per clone
                                              + summary columns
  - multiclone_classification_long.tsv       (gene_id, clone, classification)
  - multiclone_summary.tsv                   per-clone class counts table
  - plot_class_proportions.png               stacked bar of class proportions
  - plot_clones_per_class.png                per-gene: how many clones support each class
  - plot_pattern_top.png                     most common multi-clone patterns

Usage:
  python GBM_classify_multiclone.py samples.tsv genes_cds.bed out_dir \\
      --format counts --min-cov 3 --min-n-cg 15 --min-n-chg 15 \\
      --min-chg-frac 0.10

Notes:
  - Per-clone runs are independent and run sequentially by default. Use
    --threads N to parallelise (each clone gets one thread; set N to the
    number of clones you can afford to run concurrently given RAM).
  - If a clone fails, the script logs the error and continues with the rest.
  - Set --skip-existing to reuse already-completed per-clone TSVs (handy when
    re-running just the merge step or recovering from a partial failure).
"""

import os
import sys
import argparse
import subprocess
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("Error: pip install pandas matplotlib numpy")


# ---------------- Locate the per-clone classifier ----------------

def find_classifier(explicit=None):
    """
    Locate GBM_classify_CG_CHG.py (preferred) or older variants.
    Searches: --classifier flag, same dir as this script, $PATH, cwd.
    """
    candidates = []
    if explicit:
        candidates.append(explicit)
    here = os.path.dirname(os.path.abspath(__file__))
    for name in ("GBM_classify_CG_CHG.py", "GBM_classify_sites.py", "GBM_classify.py"):
        candidates.append(os.path.join(here, name))
        candidates.append(os.path.join(os.getcwd(), name))
        # also try $PATH
        which = shutil.which(name)
        if which:
            candidates.append(which)
    for c in candidates:
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    sys.exit(
        "Error: could not find GBM_classify_CG_CHG.py.\n"
        "Pass it explicitly with --classifier /path/to/GBM_classify_CG_CHG.py"
    )


# ---------------- Sample sheet ----------------

def parse_chrom_rename(spec):
    """
    Parse 'FROM1=TO1,FROM2=TO2' into a list of (FROM, TO) sorted by descending
    FROM length so longer suffixes match before shorter prefixes/substrings.
    """
    if not spec or spec == "" or spec == "-":
        return []
    rules = []
    for piece in spec.split(","):
        piece = piece.strip()
        if not piece:
            continue
        if "=" not in piece:
            sys.exit(f"Error: malformed chrom_rename rule '{piece}' (expected FROM=TO).")
        frm, to = piece.split("=", 1)
        rules.append((frm, to))
    # longest FROM first
    rules.sort(key=lambda x: -len(x[0]))
    return rules


def load_sample_sheet(path):
    """
    Read a TSV with columns: clone_name, cg_file, chg_file, [chrom_rename]
    Comment lines (#...) and blank lines are ignored.

    Backward compatibility: a legacy 5-column sheet
        clone_name, cg_file, chg_file, chh_file, [chrom_rename]
    is also accepted. The CHH file is silently ignored, and the chrom_rename
    is read from column 5 if present.

    To distinguish the two layouts on a 4-column row, we check whether column
    4 looks like an existing file (legacy CHH) or a rename rule
    (e.g. '_RagTag=_HapB' or '-' / empty).

    Returns list of dicts.
    """
    samples = []
    seen_names = set()
    with open(path) as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                sys.exit(f"Error: sample sheet line has <3 columns: {line!r}")

            clone = parts[0].strip()
            cg = parts[1].strip()
            chg = parts[2].strip()

            # Detect layout from columns 4+
            rename_spec = ""
            if len(parts) == 4:
                col4 = parts[3].strip()
                # Legacy 4-col with a CHH file in col 4 is unusual but possible:
                # if col4 is a real file path AND looks like a methylation BED,
                # treat as legacy with no rename. Otherwise treat as rename.
                if col4 and col4 != "-" and os.path.isfile(col4) and "=" not in col4:
                    # legacy: col4 was CHH, no rename
                    print(f"  [{clone}] legacy 5-col layout detected, "
                          f"ignoring CHH file: {col4}")
                    rename_spec = ""
                else:
                    rename_spec = col4
            elif len(parts) >= 5:
                # Legacy 5-col: col4 = CHH (ignored), col5 = rename
                col4 = parts[3].strip()
                if col4 and col4 != "-":
                    print(f"  [{clone}] legacy 5-col layout, ignoring CHH file: {col4}")
                rename_spec = parts[4].strip()

            if clone in seen_names:
                sys.exit(f"Error: duplicate clone_name '{clone}' in sample sheet.")
            seen_names.add(clone)
            for label, p in (("cg", cg), ("chg", chg)):
                if not os.path.isfile(p):
                    sys.exit(f"Error: {label}_file for clone '{clone}' not found: {p}")
            samples.append({
                "clone": clone,
                "cg": cg, "chg": chg,
                "rename_rules": parse_chrom_rename(rename_spec),
            })
    if not samples:
        sys.exit("Error: no samples parsed from sample sheet.")
    print(f"Sample sheet: {len(samples)} clone(s) loaded.")
    return samples


# ---------------- Chrom rename ----------------

def apply_rename_to_chrom(chrom, rules):
    """Apply the first matching rule (rules pre-sorted longest FROM first)."""
    for frm, to in rules:
        if frm in chrom:
            return chrom.replace(frm, to)
    return chrom


def write_renamed_bed(src, dst, rules):
    """
    Stream src -> dst, rewriting column 1 (chrom) by applying the rename rules.
    Pure text passthrough for everything else.
    """
    n = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            if not line.strip() or line.startswith("#"):
                fo.write(line)
                continue
            parts = line.split("\t", 1)
            chrom = apply_rename_to_chrom(parts[0], rules)
            if len(parts) == 2:
                fo.write(chrom + "\t" + parts[1])
            else:
                fo.write(chrom + "\n")
            n += 1
    return n


# ---------------- Per-clone runner ----------------

def run_one_clone(clone_info, genes_bed, classifier_path, out_dir, classifier_args,
                  skip_existing=False, work_root=None):
    """
    Run the per-clone classifier. Returns dict with status & paths.
    Designed to be called from ProcessPoolExecutor, so it must be picklable
    (only stdlib types in args).
    """
    clone = clone_info["clone"]
    per_clone_dir = os.path.join(out_dir, "per_clone")
    os.makedirs(per_clone_dir, exist_ok=True)
    out_prefix = os.path.join(per_clone_dir, clone)
    expected_tsv = f"{out_prefix}_classification.tsv"

    if skip_existing and os.path.isfile(expected_tsv):
        return {"clone": clone, "status": "skipped", "tsv": expected_tsv, "log": ""}

    # If chrom rename is needed, write temp renamed BEDs in the work_root.
    # Otherwise pass the source paths directly to the classifier.
    rules = clone_info["rename_rules"]
    cg_in, chg_in = clone_info["cg"], clone_info["chg"]
    tmpdir = None
    try:
        if rules:
            tmpdir = tempfile.mkdtemp(prefix=f"mc_{clone}_", dir=work_root)
            cg_in_r  = os.path.join(tmpdir, "cg.bed")
            chg_in_r = os.path.join(tmpdir, "chg.bed")
            write_renamed_bed(clone_info["cg"],  cg_in_r,  rules)
            write_renamed_bed(clone_info["chg"], chg_in_r, rules)
            cg_in, chg_in = cg_in_r, chg_in_r

        cmd = [
            sys.executable, classifier_path,
            genes_bed, cg_in, chg_in, out_prefix,
        ] + list(classifier_args)

        log_path = f"{out_prefix}.log"
        with open(log_path, "w") as logf:
            logf.write("CMD: " + " ".join(cmd) + "\n\n")
            logf.flush()
            proc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT)

        if proc.returncode != 0 or not os.path.isfile(expected_tsv):
            return {"clone": clone, "status": "failed",
                    "tsv": None, "log": log_path,
                    "returncode": proc.returncode}
        return {"clone": clone, "status": "ok", "tsv": expected_tsv, "log": log_path}
    finally:
        if tmpdir and os.path.isdir(tmpdir):
            shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------- Merge ----------------

CLASS_ORDER = ["gbM", "teM", "UM", "Unclassified"]
COLOR_MAP   = {"gbM": "#4A90D9", "teM": "#E57373",
               "UM": "#BDBDBD", "Unclassified": "#EEEEEE"}


def load_per_clone(tsv_path, clone):
    """Read per-clone classification TSV, keep gene_id and Classification."""
    df = pd.read_csv(tsv_path, sep="\t", usecols=["Gene_ID", "Classification"])
    df = df.rename(columns={"Gene_ID": "gene_id",
                            "Classification": f"class_{clone}"})
    return df


def make_pattern(class_counts, top_k=2):
    """
    Compact pattern string from a Counter.
    e.g. Counter({'gbM': 22, 'teM': 1}) -> '22xgbM,1xteM'
    Sorted by descending count, then by CLASS_ORDER for ties.
    """
    order_idx = {c: i for i, c in enumerate(CLASS_ORDER)}
    items = sorted(class_counts.items(),
                   key=lambda kv: (-kv[1], order_idx.get(kv[0], 99)))
    return ",".join(f"{n}x{c}" for c, n in items)


def summarise_row(row, clone_cols):
    """Compute summary columns for one gene across clones."""
    classes = [row[c] for c in clone_cols]
    counts = Counter(classes)
    n_total = len(classes)
    n_unc = counts.get("Unclassified", 0)
    n_classified = n_total - n_unc

    # restrict to non-Unclassified for "biological" agreement metrics
    bio_counts = Counter({k: v for k, v in counts.items() if k != "Unclassified"})
    n_distinct = len(bio_counts)

    if bio_counts:
        dominant, n_dominant = bio_counts.most_common(1)[0]
    else:
        dominant, n_dominant = "Unclassified", 0

    is_stable = (n_classified > 0 and n_distinct == 1)
    # singleton: n_classified-1 clones agree on one class, exactly 1 differs (also classified)
    is_singleton = False
    minority_class = ""
    if n_classified >= 2 and n_distinct == 2:
        # one class has n_classified-1, the other has 1
        if n_dominant == n_classified - 1:
            is_singleton = True
            minority_class = next(c for c in bio_counts if c != dominant)

    class_set = "|".join(sorted(bio_counts.keys()))
    pattern = make_pattern(counts)

    return pd.Series({
        "n_clones_total":       n_total,
        "n_clones_classified":  n_classified,
        "n_distinct_classes":   n_distinct,
        "dominant_class":       dominant,
        "n_dominant":           n_dominant,
        "is_stable":            is_stable,
        "is_singleton_change":  is_singleton,
        "minority_class":       minority_class,
        "class_set":            class_set,
        "pattern":              pattern,
    })


def merge_clones(per_clone_results, clone_order):
    """
    Outer-join per-clone classifications on gene_id.
    Genes missing from a clone (shouldn't happen if BED is shared) become NaN
    and are filled with 'Unclassified' for the summary logic.
    """
    if not per_clone_results:
        sys.exit("Error: no successful per-clone classifications to merge.")

    merged = None
    for clone in clone_order:
        if clone not in per_clone_results:
            continue
        df = per_clone_results[clone]
        merged = df if merged is None else merged.merge(df, on="gene_id", how="outer")

    clone_cols = [f"class_{c}" for c in clone_order if c in per_clone_results]
    merged[clone_cols] = merged[clone_cols].fillna("Unclassified")

    summary = merged.apply(lambda row: summarise_row(row, clone_cols), axis=1)
    wide = pd.concat([merged, summary], axis=1)

    # long form for ggplot/seaborn
    long_df = wide[["gene_id"] + clone_cols].melt(
        id_vars="gene_id", var_name="clone", value_name="classification"
    )
    long_df["clone"] = long_df["clone"].str.replace("class_", "", regex=False)

    return wide, long_df, clone_cols


# ---------------- Plots ----------------

def plot_class_proportions(wide, clone_cols, outpath):
    """Stacked bar: per-clone proportion of gbM/teM/UM/Unclassified."""
    counts = pd.DataFrame({
        c.replace("class_", ""): wide[c].value_counts()
        for c in clone_cols
    }).T
    counts = counts.reindex(columns=CLASS_ORDER, fill_value=0)
    props = counts.div(counts.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(max(8, 0.4 * len(props) + 4), 5.5))
    bottom = np.zeros(len(props))
    for cl in CLASS_ORDER:
        vals = props[cl].values
        ax.bar(props.index, vals, bottom=bottom,
               color=COLOR_MAP[cl], label=cl, edgecolor="white", linewidth=0.4)
        bottom += vals
    ax.set_ylabel("Proportion of genes")
    ax.set_title("Per-clone classification proportions", fontweight="bold")
    ax.set_ylim(0, 1)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def plot_clones_per_class(wide, clone_cols, outpath):
    """For each class, distribution of 'how many clones called this gene <class>'."""
    n_clones = len(clone_cols)
    counts_per_gene = {cl: [] for cl in ("gbM", "teM", "UM")}
    for cl in counts_per_gene:
        counts_per_gene[cl] = (wide[clone_cols] == cl).sum(axis=1).values

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for ax, cl in zip(axes, ("gbM", "teM", "UM")):
        vals = counts_per_gene[cl]
        # only show genes that were called this class in at least one clone
        vals_nz = vals[vals > 0]
        ax.hist(vals_nz, bins=np.arange(0.5, n_clones + 1.5, 1),
                color=COLOR_MAP[cl], edgecolor="white", alpha=0.9)
        ax.set_xlabel(f"# clones calling gene as {cl}")
        ax.set_title(f"{cl}: {len(vals_nz):,} genes called in \u22651 clone",
                     fontweight="bold", fontsize=10)
        ax.axvline(n_clones, color="black", linestyle="--", linewidth=1,
                   label=f"all {n_clones} clones")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Number of genes")
    plt.suptitle("Per-class clone-support distribution", fontweight="bold")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


def plot_top_patterns(wide, outpath, top_n=20):
    """Bar of the most common multi-clone patterns."""
    pat = wide["pattern"].value_counts().head(top_n)
    fig, ax = plt.subplots(figsize=(10, max(5, 0.3 * len(pat) + 2)))
    ax.barh(pat.index[::-1], pat.values[::-1],
            color="#4A90D9", edgecolor="white")
    ax.set_xlabel("Number of genes")
    ax.set_title(f"Top {len(pat)} multi-clone classification patterns",
                 fontweight="bold")
    for i, v in enumerate(pat.values[::-1]):
        ax.text(v, i, f" {v:,}", va="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {outpath}")


# ---------------- Main ----------------

def main():
    parser = argparse.ArgumentParser(
        description="Multi-clone gbM/teM/UM classification + differential summary.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("sample_sheet", help="TSV: clone_name, cg, chg, [chrom_rename]")
    parser.add_argument("genes_bed",    help="Shared genes/CDS BED (chrom, start, end, gene_id)")
    parser.add_argument("output_dir",   help="Output directory (created if needed)")
    parser.add_argument("--classifier", default=None,
                        help="Path to GBM_classify_sites.py (auto-detected if omitted).")
    parser.add_argument("--threads",    type=int, default=1,
                        help="Number of clones to run concurrently.")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Reuse per-clone TSVs that already exist on disk.")
    parser.add_argument("--work-root",  default=None,
                        help="Where to put temporary renamed BEDs (default: system tmp).")

    # Forwarded to per-clone classifier
    parser.add_argument("--format",        choices=["percent", "counts"], default="percent")
    parser.add_argument("--min-cov",       type=int,   default=2)
    parser.add_argument("--min-n-cg",      type=int,   default=15)
    parser.add_argument("--min-n-chg",     type=int,   default=15)
    parser.add_argument("--alpha",         type=float, default=0.05)
    parser.add_argument("--min-chg-frac",  type=float, default=0.10)
    parser.add_argument("--max-cg-frac-um", type=float, default=0.10,
                        help="Effect-size guard for UM (forwarded to per-clone classifier).")

    args = parser.parse_args()

    # Locate inputs/outputs
    classifier = find_classifier(args.classifier)
    print(f"Using classifier: {classifier}")
    if not os.path.isfile(args.genes_bed):
        sys.exit(f"Error: genes BED not found: {args.genes_bed}")
    os.makedirs(args.output_dir, exist_ok=True)

    samples = load_sample_sheet(args.sample_sheet)

    # Build forwarded classifier args (only those that differ from defaults
    # are technically necessary, but forwarding them all is harmless and explicit).
    fwd = [
        "--format",         args.format,
        "--min-cov",        str(args.min_cov),
        "--min-n-cg",       str(args.min_n_cg),
        "--min-n-chg",      str(args.min_n_chg),
        "--alpha",          str(args.alpha),
        "--min-chg-frac",   str(args.min_chg_frac),
        "--max-cg-frac-um", str(args.max_cg_frac_um),
    ]

    # ---- Run per-clone ----
    print(f"\nRunning per-clone classification (threads={args.threads})...")
    results = {}
    if args.threads <= 1:
        for s in samples:
            print(f"  [{s['clone']}] starting...")
            r = run_one_clone(s, args.genes_bed, classifier, args.output_dir,
                              fwd, args.skip_existing, args.work_root)
            print(f"  [{s['clone']}] {r['status']}")
            results[s["clone"]] = r
    else:
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            future_to_clone = {
                ex.submit(run_one_clone, s, args.genes_bed, classifier,
                          args.output_dir, fwd, args.skip_existing,
                          args.work_root): s["clone"]
                for s in samples
            }
            for fut in as_completed(future_to_clone):
                clone = future_to_clone[fut]
                try:
                    r = fut.result()
                except Exception as e:
                    r = {"clone": clone, "status": "exception", "tsv": None,
                         "log": "", "error": str(e)}
                print(f"  [{clone}] {r['status']}")
                results[clone] = r

    failed = [c for c, r in results.items() if r["status"] not in ("ok", "skipped")]
    ok = [c for c, r in results.items() if r["status"] in ("ok", "skipped")]
    if failed:
        print(f"\nWARNING: {len(failed)} clone(s) failed: {failed}")
        for c in failed:
            print(f"  See log: {results[c].get('log', '<no log>')}")
    if not ok:
        sys.exit("All clones failed; cannot merge.")

    # ---- Load per-clone TSVs ----
    print(f"\nMerging {len(ok)} clone(s)...")
    per_clone_dfs = {}
    for clone in ok:
        per_clone_dfs[clone] = load_per_clone(results[clone]["tsv"], clone)

    clone_order = [s["clone"] for s in samples if s["clone"] in per_clone_dfs]
    wide, long_df, clone_cols = merge_clones(per_clone_dfs, clone_order)

    # ---- Save merged tables ----
    out_wide = os.path.join(args.output_dir, "multiclone_classification_wide.tsv")
    out_long = os.path.join(args.output_dir, "multiclone_classification_long.tsv")
    wide.to_csv(out_wide, sep="\t", index=False)
    long_df.to_csv(out_long, sep="\t", index=False)
    print(f"  Saved: {out_wide}")
    print(f"  Saved: {out_long}")

    # ---- Per-clone summary table ----
    summ = pd.DataFrame({
        c.replace("class_", ""): wide[c].value_counts()
        for c in clone_cols
    }).T.reindex(columns=CLASS_ORDER, fill_value=0)
    summ.index.name = "clone"
    out_summary = os.path.join(args.output_dir, "multiclone_summary.tsv")
    summ.to_csv(out_summary, sep="\t")
    print(f"  Saved: {out_summary}")

    # ---- Console summary ----
    n_total = len(wide)
    print(f"\n=== Multi-clone summary ({len(clone_cols)} clones, {n_total:,} genes) ===")
    print(f"  Stable across all classified clones: "
          f"{wide['is_stable'].sum():,} ({100*wide['is_stable'].mean():.1f}%)")
    print(f"  Singleton class change (1 clone differs from rest): "
          f"{wide['is_singleton_change'].sum():,} ({100*wide['is_singleton_change'].mean():.1f}%)")
    print(f"  Multi-class (>=2 distinct classes among classified): "
          f"{(wide['n_distinct_classes'] >= 2).sum():,}")
    print(f"  Highly unstable (>=3 distinct classes): "
          f"{(wide['n_distinct_classes'] >= 3).sum():,}")
    print(f"  Genes Unclassified in ALL clones: "
          f"{(wide['n_clones_classified'] == 0).sum():,}")
    print(f"  Genes classified in ALL clones: "
          f"{(wide['n_clones_classified'] == len(clone_cols)).sum():,}")

    print("\nTop 10 multi-clone patterns:")
    for pat, n in wide["pattern"].value_counts().head(10).items():
        print(f"  {n:7,d}  {pat}")

    # ---- Plots ----
    print("\nGenerating plots...")
    plot_class_proportions(wide, clone_cols,
                           os.path.join(args.output_dir, "plot_class_proportions.png"))
    plot_clones_per_class(wide, clone_cols,
                          os.path.join(args.output_dir, "plot_clones_per_class.png"))
    plot_top_patterns(wide,
                      os.path.join(args.output_dir, "plot_pattern_top.png"))

    print("\nDone.")


if __name__ == "__main__":
    main()
