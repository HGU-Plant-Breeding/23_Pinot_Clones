#!/usr/bin/env python3
"""
Script Name: call_switching_genes.py
Description: Apply a conservative "switching gene" definition to the output
             of GBM_classify_multiclone.py.

             Uses CG and CHG only (CHH dropped due to noise in ONT calls).

A gene is a switching gene if and only if:
  1. It is classified (gbM/teM/UM) in ALL clones (no Unclassified anywhere).
     Note: with the --max-cg-frac-um guard in the classifier, Unclassified
     includes both low-coverage genes AND genes with intermediate CG methylation
     that failed the UM effect-size guard. Both are excluded here.
  2. At least 2 distinct biological classes appear across clones.
  3. The minority class is supported by >= MIN_MINORITY clones (default 2).
  4. The mean per-gene methylation fraction in the relevant context differs
     between the majority-class clones and the minority-class clones by at
     least DELTA, where:
        gbM <-> UM   : delta on CG_frac      (default 0.30)
        gbM <-> teM  : delta on CHG_frac     (default 0.10)
        UM  <-> teM  : delta on CG_frac AND CHG_frac
        gbM/teM/UM   : when 3 classes appear, the gene is flagged as
                       'multi_switch' and effect-size guards are applied
                       pairwise to the dominant class vs. each other class.
                       All pairs must pass to be called a switch.

Inputs:
  - <multiclone_dir>/multiclone_classification_wide.tsv  (from the multiclone driver)
  - <multiclone_dir>/per_clone/<clone>_classification.tsv (per-clone TSVs, also
    from the multiclone driver -- needed to pull per-clone fractions)

Outputs:
  - <multiclone_dir>/switching_genes.tsv          conservative switching genes
                                                   with full annotation
  - <multiclone_dir>/switching_genes_summary.txt  count breakdown by transition

Usage:
  python call_switching_genes.py multiclone_out [options]
"""

import os
import sys
import argparse
from collections import Counter

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Error: pip install pandas numpy")


CLASS_ORDER = ["gbM", "teM", "UM"]   # Unclassified excluded by criterion 1
FRAC_COLS   = ["CG_frac", "CHG_frac"]


def load_per_clone_fracs(per_clone_dir, clones):
    """
    For each clone, load gene_id + CG_frac/CHG_frac.
    Returns dict clone -> DataFrame indexed by gene_id.
    """
    out = {}
    for clone in clones:
        path = os.path.join(per_clone_dir, f"{clone}_classification.tsv")
        if not os.path.isfile(path):
            sys.exit(f"Error: per-clone TSV missing for {clone}: {path}")
        df = pd.read_csv(path, sep="\t",
                         usecols=["Gene_ID"] + FRAC_COLS,
                         dtype={"Gene_ID": str})
        for col in FRAC_COLS:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df = df.rename(columns={"Gene_ID": "gene_id"}).set_index("gene_id")
        out[clone] = df
    return out


def gather_fracs(gene_id, clones_in_group, per_clone_fracs, context):
    """Mean of <context>_frac across the listed clones for one gene."""
    col = f"{context}_frac"
    vals = []
    for c in clones_in_group:
        if gene_id in per_clone_fracs[c].index:
            v = per_clone_fracs[c].loc[gene_id, col]
            if pd.notna(v):
                vals.append(float(v))
    return float(np.mean(vals)) if vals else float("nan")


def transition_delta_ok(maj_class, min_class, fracs_maj, fracs_min,
                         d_cg, d_chg):
    """
    Test whether the mean-fraction shift between the majority and minority
    groups is large enough for the given class transition.
    fracs_maj and fracs_min are dicts with keys 'CG', 'CHG'.
    Returns (ok: bool, reason: str) where reason is the tested rule.
    """
    pair = frozenset({maj_class, min_class})

    if pair == frozenset({"gbM", "UM"}):
        delta = abs(fracs_maj["CG"] - fracs_min["CG"])
        return (delta >= d_cg,
                f"|dCG|={delta:.3f} vs >= {d_cg}")

    if pair == frozenset({"gbM", "teM"}):
        d_chg_obs = abs(fracs_maj["CHG"] - fracs_min["CHG"])
        return (d_chg_obs >= d_chg,
                f"|dCHG|={d_chg_obs:.3f} vs >= {d_chg}")

    if pair == frozenset({"UM", "teM"}):
        d_cg_obs  = abs(fracs_maj["CG"]  - fracs_min["CG"])
        d_chg_obs = abs(fracs_maj["CHG"] - fracs_min["CHG"])
        ok_cg  = d_cg_obs  >= d_cg
        ok_chg = d_chg_obs >= d_chg
        return (ok_cg and ok_chg,
                f"|dCG|={d_cg_obs:.3f}>={d_cg} AND "
                f"|dCHG|={d_chg_obs:.3f}>={d_chg}")

    return (False, f"unknown pair {pair}")


def main():
    ap = argparse.ArgumentParser(
        description="Conservative switching-gene caller for multiclone gbM/teM/UM data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("multiclone_dir",
                    help="Directory containing multiclone_classification_wide.tsv "
                         "and per_clone/*.tsv (output of GBM_classify_multiclone.py).")
    ap.add_argument("--min-minority", type=int, default=2,
                    help="Minimum clones supporting the minority class. "
                         "2 = at least one replication; 3 = stricter.")
    ap.add_argument("--delta-cg",  type=float, default=0.30,
                    help="Required |delta CG_frac| between majority and minority "
                         "groups for gbM<->UM and UM<->teM transitions.")
    ap.add_argument("--delta-chg", type=float, default=0.10,
                    help="Required |delta CHG_frac| for gbM<->teM and UM<->teM "
                         "transitions.")
    ap.add_argument("--min-classified", type=int, default=None,
                    help="Minimum number of clones that must give a biological "
                         "classification (gbM/teM/UM, not Unclassified) for the "
                         "gene to be considered. Default: all clones (strict). "
                         "Recommended: total_clones - 3 to allow 1-3 coverage gaps "
                         "or intermediate-CG ambiguities without losing biologically "
                         "clear cases.")
    ap.add_argument("--output", default=None,
                    help="Output TSV path (default: <multiclone_dir>/switching_genes.tsv).")
    args = ap.parse_args()

    wide_path = os.path.join(args.multiclone_dir, "multiclone_classification_wide.tsv")
    if not os.path.isfile(wide_path):
        sys.exit(f"Error: cannot find {wide_path}")

    per_clone_dir = os.path.join(args.multiclone_dir, "per_clone")
    if not os.path.isdir(per_clone_dir):
        sys.exit(f"Error: cannot find {per_clone_dir}")

    print(f"Loading wide table: {wide_path}")
    wide = pd.read_csv(wide_path, sep="\t")
    SUMMARY_COLS_TO_SKIP = {"class_set"}
    clone_cols = [c for c in wide.columns
                  if c.startswith("class_") and c not in SUMMARY_COLS_TO_SKIP]
    clones = [c[len("class_"):] for c in clone_cols]
    n_clones = len(clones)
    print(f"  {len(wide):,} genes, {n_clones} clones.")

    # ---- Criterion 1: classified in at least min_classified clones ----
    min_classified = args.min_classified if args.min_classified is not None else n_clones
    if min_classified > n_clones:
        sys.exit(f"Error: --min-classified ({min_classified}) > total clones ({n_clones}).")

    n_classified_per_gene = (wide[clone_cols] != "Unclassified").sum(axis=1)
    fully_classified = n_classified_per_gene >= min_classified
    if min_classified == n_clones:
        print(f"  Classified in all {n_clones} clones: "
              f"{fully_classified.sum():,} / {len(wide):,} "
              f"({100*fully_classified.mean():.1f}%)")
    else:
        print(f"  Classified in >= {min_classified} of {n_clones} clones: "
              f"{fully_classified.sum():,} / {len(wide):,} "
              f"({100*fully_classified.mean():.1f}%)")

    # ---- Criterion 2: at least 2 distinct biological classes ----
    multi_class = wide["n_distinct_classes"] >= 2
    candidates = wide[fully_classified & multi_class].copy()
    print(f"  Plus >= 2 distinct classes:           {len(candidates):,}")

    if len(candidates) == 0:
        print("\nNo candidate genes after basic filtering. Nothing to do.")
        return

    # ---- Criterion 3: minority supported by >= min_minority clones ----
    def smallest_nonzero(row):
        counts = Counter(row[clone_cols])
        bio = {k: v for k, v in counts.items() if k != "Unclassified" and v > 0}
        return min(bio.values()) if bio else 0
    candidates["n_minority"] = candidates.apply(smallest_nonzero, axis=1)
    crit3 = candidates["n_minority"] >= args.min_minority
    candidates = candidates[crit3].copy()
    print(f"  Plus minority class supported by >= {args.min_minority} clones: "
          f"{len(candidates):,}")

    if len(candidates) == 0:
        print("\nNo candidate genes after replication filter.")
        return

    # ---- Criterion 4: effect-size guard on the relevant fractions ----
    print(f"\nLoading per-clone fractions for effect-size test...")
    per_clone_fracs = load_per_clone_fracs(per_clone_dir, clones)

    print(f"Applying effect-size guards "
          f"(dCG>={args.delta_cg}, dCHG>={args.delta_chg})...")

    keep_rows = []
    extra = []

    for _, row in candidates.iterrows():
        gid = row["gene_id"]
        per_clone_class = {c: row[f"class_{c}"] for c in clones}
        class_to_clones = {}
        for c, cl in per_clone_class.items():
            if cl == "Unclassified":
                continue
            class_to_clones.setdefault(cl, []).append(c)

        if len(class_to_clones) < 2:
            continue

        present_classes = sorted(class_to_clones.keys(),
                                  key=lambda k: -len(class_to_clones[k]))
        majority = present_classes[0]
        majority_clones = class_to_clones[majority]

        all_pairs_pass = True
        rule_details = []
        for other in present_classes[1:]:
            other_clones = class_to_clones[other]
            fracs_maj = {ctx: gather_fracs(gid, majority_clones, per_clone_fracs, ctx)
                         for ctx in ("CG", "CHG")}
            fracs_oth = {ctx: gather_fracs(gid, other_clones, per_clone_fracs, ctx)
                         for ctx in ("CG", "CHG")}
            ok, why = transition_delta_ok(
                majority, other, fracs_maj, fracs_oth,
                args.delta_cg, args.delta_chg
            )
            rule_details.append(f"{majority}vs{other}: {why} -> {'PASS' if ok else 'FAIL'}")
            if not ok:
                all_pairs_pass = False

        if all_pairs_pass:
            keep_rows.append(gid)
            primary_transition = (f"{majority}<->{present_classes[1]}"
                                  if len(present_classes) >= 2 else "")
            extra.append({
                "gene_id":            gid,
                "primary_transition": primary_transition,
                "is_multi_switch":    len(present_classes) >= 3,
                "majority_class":     majority,
                "n_majority":         len(majority_clones),
                "minority_class_set": "|".join(present_classes[1:]),
                "delta_check":        " ; ".join(rule_details),
                "mean_CG_majority":   fracs_maj["CG"],
                "mean_CHG_majority":  fracs_maj["CHG"],
            })

    print(f"  Plus effect-size guards: {len(keep_rows):,}")

    # ---- Compose final output ----
    extra_df = pd.DataFrame(extra)
    if len(extra_df) == 0:
        print("\nNo switching genes survived the conservative filter. "
              "Try relaxing --min-minority or the deltas, but understand the trade-off.")
        return

    final = candidates.merge(extra_df, on="gene_id", how="inner")

    front = ["gene_id", "primary_transition", "is_multi_switch",
             "majority_class", "n_majority",
             "minority_class_set", "n_minority",
             "n_distinct_classes", "pattern", "class_set",
             "mean_CG_majority", "mean_CHG_majority",
             "delta_check"]
    front = [c for c in front if c in final.columns]
    rest  = [c for c in final.columns if c not in front]
    final = final[front + rest]

    out_path = args.output or os.path.join(args.multiclone_dir, "switching_genes.tsv")
    final.to_csv(out_path, sep="\t", index=False, float_format="%.4f")
    print(f"\nSaved: {out_path}  ({len(final):,} switching genes)")

    # ---- Summary by transition ----
    summary_path = os.path.join(args.multiclone_dir, "switching_genes_summary.txt")
    lines = []
    lines.append(f"Conservative switching genes: {len(final):,}")
    lines.append(f"  Filters: min_minority>={args.min_minority}, "
                 f"dCG>={args.delta_cg}, dCHG>={args.delta_chg}")
    lines.append("")
    lines.append("By primary transition (majority <-> largest minority):")
    trans_counts = final["primary_transition"].value_counts()
    for t, n in trans_counts.items():
        lines.append(f"  {t}: {n:,}")
    lines.append("")
    lines.append(f"Multi-switch genes (>=3 classes):   "
                 f"{int(final['is_multi_switch'].sum()):,}")
    lines.append(f"Two-class switches:                 "
                 f"{int((~final['is_multi_switch']).sum()):,}")
    text = "\n".join(lines)
    with open(summary_path, "w") as f:
        f.write(text + "\n")
    print(f"Saved: {summary_path}")
    print("\n" + text)


if __name__ == "__main__":
    main()
