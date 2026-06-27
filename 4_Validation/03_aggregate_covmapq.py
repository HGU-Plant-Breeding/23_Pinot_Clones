#!/usr/bin/env python3
# =============================================================
# 03_aggregate_covmapq.py
# Aggregate coverage and MAPQ results from all 23 clones
# Coverage file columns: clone strategy region mean_depth mean_mapq pct_mapq0
# Callipo et al. 2026
# =============================================================

import pandas as pd
import numpy as np
import os
import glob

WORK_DIR = "/"
COV_DIR  = f"{WORK_DIR}/results/covmapq"
OUT_DIR  = f"{WORK_DIR}/results/aggregated"
os.makedirs(OUT_DIR, exist_ok=True)

CLONES = [
    "1-32", "1-47", "20-13-28-8", "20-16", "20-19-19", "20-20",
    "20-24-3", "20-25-2", "20-27-3", "2-2-13", "2-6-41", "2-9-35",
    "2-9-40", "37-1", "4000-1", "4000-3", "4047-1", "404",
    "Kastenholz", "Probstei", "RC", "S1-1", "S1-2"
]

# -------------------------------------------------------
# 1. Aggregate coverage tables
# Columns: clone strategy region mean_depth mean_mapq pct_mapq0
# -------------------------------------------------------
print("Aggregating coverage data...")

cov_cols = ["clone", "strategy", "region", "mean_depth", "mean_mapq", "pct_mapq0"]
all_cov  = []

for clone in CLONES:
    for strategy in ["masked", "unmasked"]:
        f = f"{COV_DIR}/{clone}_{strategy}_coverage.tsv"
        if not os.path.exists(f):
            print(f"  WARNING: missing {f}")
            continue
        df = pd.read_csv(f, sep="\t", header=None, names=cov_cols)
        all_cov.append(df)

cov_df = pd.concat(all_cov, ignore_index=True)
cov_df.to_csv(f"{OUT_DIR}/coverage_all_clones.tsv", sep="\t", index=False)
print(f"  -> coverage_all_clones.tsv: {len(cov_df)} rows")

# Summary: mean ± sd across clones per strategy x region
cov_summary = cov_df.groupby(["strategy", "region"]).agg(
    mean_depth   = ("mean_depth",  "mean"),
    sd_depth     = ("mean_depth",  "std"),
    mean_mapq    = ("mean_mapq",   "mean"),
    sd_mapq      = ("mean_mapq",   "std"),
    mean_pct_mapq0 = ("pct_mapq0", "mean"),
    sd_pct_mapq0   = ("pct_mapq0", "std"),
    n_clones     = ("clone",       "count")
).reset_index()

cov_summary.to_csv(f"{OUT_DIR}/coverage_summary.tsv", sep="\t", index=False)
print(f"  -> coverage_summary.tsv")
print(cov_summary.to_string(index=False))

# -------------------------------------------------------
# 2. Aggregate MAPQ distributions (for histogram)
# -------------------------------------------------------
print("\nAggregating MAPQ distributions...")

mapq_cols = ["clone", "strategy", "mapq"]

def load_mapq(pattern, region_type):
    files = glob.glob(pattern)
    if not files:
        print(f"  WARNING: no files matching {pattern}")
        return pd.DataFrame()
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep="\t", header=None, names=mapq_cols)
        df["region"] = region_type
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

mapq_masked_roh      = load_mapq(f"{COV_DIR}/*_masked_mapq_RoH.tsv",      "RoH")
mapq_masked_nonroh   = load_mapq(f"{COV_DIR}/*_masked_mapq_nonRoH.tsv",   "nonRoH")
mapq_unmasked_roh    = load_mapq(f"{COV_DIR}/*_unmasked_mapq_RoH.tsv",    "RoH")
mapq_unmasked_nonroh = load_mapq(f"{COV_DIR}/*_unmasked_mapq_nonRoH.tsv", "nonRoH")

mapq_all = pd.concat([
    mapq_masked_roh, mapq_masked_nonroh,
    mapq_unmasked_roh, mapq_unmasked_nonroh
], ignore_index=True)

# MAPQ summary per strategy x region
mapq_summary = mapq_all.groupby(["strategy", "region"]).apply(
    lambda x: pd.Series({
        "n_reads"        : len(x),
        "pct_mapq0"      : (x["mapq"] == 0).mean() * 100,
        "pct_mapq_lt10"  : (x["mapq"] < 10).mean() * 100,
        "median_mapq"    : x["mapq"].median(),
        "mean_mapq"      : x["mapq"].mean()
    })
).reset_index()

mapq_summary.to_csv(f"{OUT_DIR}/mapq_summary.tsv", sep="\t", index=False)
print(f"  -> mapq_summary.tsv")
print(mapq_summary.to_string(index=False))

# MAPQ frequency table for plotting
mapq_freq = mapq_all.groupby(["strategy", "region", "mapq"]).size().reset_index(name="count")
mapq_freq["total"]    = mapq_freq.groupby(["strategy", "region"])["count"].transform("sum")
mapq_freq["fraction"] = mapq_freq["count"] / mapq_freq["total"]
mapq_freq.to_csv(f"{OUT_DIR}/mapq_frequency.tsv", sep="\t", index=False)
print(f"  -> mapq_frequency.tsv: {len(mapq_freq)} rows")

print("\nAggregation complete.")
