#!/usr/bin/env python3
"""
GO enrichment analysis for gbM/teM/UM shifting genes.
Uses Fisher's exact test with BH FDR correction.

Inputs:
  - grapevine_HMM_results.emapper.annotations  (eggNOG-mapper output)
  - 20-13_classification.tsv                   (background: genes with nCG >= 15)
  - switching_genes.tsv                        (shifting genes with transition type)
"""

import sys
import re
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ─── ID normalisation ────────────────────────────────────────────────────────

def normalise_id(gene_id: str) -> str:
    """
    Bring any ID variant to a common key:
      mikado.PN10_HapAG1017   -> mikado.PN10_Hap_1G1017
      mikado.PN10_HapBG1017   -> mikado.PN10_Hap_2G1017
      mikado.PN10_Hap_1G2.1   -> mikado.PN10_Hap_1G2      (strip isoform)
    """
    gid = gene_id.strip()
    # strip isoform suffix (.1, .2, …)
    gid = re.sub(r'\.\d+$', '', gid)
    # HapA -> Hap_1 , HapB -> Hap_2
    gid = re.sub(r'HapA', 'Hap_1', gid)
    gid = re.sub(r'HapB', 'Hap_2', gid)
    return gid


# ─── Load eggNOG GO mappings ──────────────────────────────────────────────────

def load_go_map(emapper_file: str) -> dict[str, set[str]]:
    """
    Returns {normalised_gene_id: {GO:XXXXXXX, …}}
    Skips comment lines (##) and the header line (#query …).
    """
    go_map: dict[str, set[str]] = {}
    with open(emapper_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 10:
                continue
            gene_id = normalise_id(parts[0])
            go_col  = parts[9].strip()
            if go_col and go_col != '-':
                go_map[gene_id] = set(go_col.split(','))
    print(f"  eggNOG: loaded GO terms for {len(go_map):,} genes")
    return go_map


# ─── Load background genes ────────────────────────────────────────────────────

def load_background(classification_file: str, min_sites: int = 15) -> set[str]:
    """
    Background = all genes in 20-13_classification.tsv with CG_nSites >= min_sites
    (mirrors the nCG >= 15 filter used in the revised analysis).
    """
    df = pd.read_csv(classification_file, sep='\t')
    df['gene_norm'] = df['Gene_ID'].apply(normalise_id)
    bg = set(df.loc[df['CG_nSites'] >= min_sites, 'gene_norm'])
    print(f"  Background: {len(bg):,} genes (nCG >= {min_sites})")
    return bg


# ─── Load switching genes ─────────────────────────────────────────────────────

def load_switching_genes(switching_file: str) -> pd.DataFrame:
    df = pd.read_csv(switching_file, sep='\t')
    df['gene_norm'] = df['gene_id'].apply(normalise_id)
    return df


# ─── Fisher's exact test GO enrichment ───────────────────────────────────────

def go_enrichment(
    query_genes:      set[str],
    background_genes: set[str],
    go_map:           dict[str, set[str]],
    min_count:        int = 3,
) -> pd.DataFrame:
    """
    For each GO term, run a 2x2 Fisher exact test:

              | in query | not in query
    annotated |    a     |      b
    not annot.|    c     |      d

    background_genes already includes query_genes (as expected for ORA).
    """
    # restrict to genes that have GO annotations AND are in background
    bg_with_go    = {g for g in background_genes if g in go_map}
    query_with_go = {g for g in query_genes       if g in go_map}

    N  = len(bg_with_go)       # total annotated background
    nq = len(query_with_go)    # annotated query genes

    if nq == 0:
        print("  WARNING: no query genes found in GO map — check ID matching")
        return pd.DataFrame()

    # build term -> gene sets
    term_to_bg    : dict[str, set[str]] = defaultdict(set)
    term_to_query : dict[str, set[str]] = defaultdict(set)

    for g in bg_with_go:
        for term in go_map[g]:
            term_to_bg[term].add(g)

    for g in query_with_go:
        for term in go_map[g]:
            term_to_query[term].add(g)

    rows = []
    for term, q_genes in term_to_query.items():
        a = len(q_genes)
        if a < min_count:
            continue
        b = len(term_to_bg[term]) - a   # bg hits outside query
        c = nq - a                       # query genes without this term
        d = N - nq - b                   # bg genes without term, not in query
        oddsratio, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        rows.append({
            'GO_term':    term,
            'count':      a,
            'query_size': nq,
            'bg_count':   len(term_to_bg[term]),
            'bg_size':    N,
            'fold_enr':   (a / nq) / (len(term_to_bg[term]) / N),
            'odds_ratio': oddsratio,
            'pvalue':     pval,
            'genes':      ','.join(sorted(q_genes)),
        })

    if not rows:
        return pd.DataFrame()

    result = pd.DataFrame(rows)
    # BH FDR correction
    _, qvals, _, _ = multipletests(result['pvalue'], method='fdr_bh')
    result['qvalue'] = qvals
    result = result.sort_values('qvalue')
    return result


# ─── Simple GO namespace lookup (avoids requiring goatools) ──────────────────

# We'll infer namespace from the GO term prefix counts in data (good enough)
# and add a lightweight OBO parser for proper labels if goatools is available.

def try_load_go_labels(obo_path: str | None) -> dict[str, tuple[str, str]]:
    """
    Returns {GO:XXXXXXX: (name, namespace_abbr)} if an OBO file is available.
    namespace_abbr: BP / MF / CC
    """
    if obo_path is None:
        return {}
    labels: dict[str, tuple[str, str]] = {}
    ns_map = {
        'biological_process':  'BP',
        'molecular_function':  'MF',
        'cellular_component':  'CC',
    }
    try:
        current_id = current_name = current_ns = None
        with open(obo_path) as fh:
            for line in fh:
                line = line.rstrip()
                if line == '[Term]':
                    if current_id and current_name and current_ns:
                        labels[current_id] = (current_name, ns_map.get(current_ns, current_ns))
                    current_id = current_name = current_ns = None
                elif line.startswith('id: GO:'):
                    current_id = line[4:]
                elif line.startswith('name: '):
                    current_name = line[6:]
                elif line.startswith('namespace: '):
                    current_ns = line[11:]
        print(f"  OBO: loaded labels for {len(labels):,} terms")
    except Exception as e:
        print(f"  OBO load failed: {e}")
    return labels


# ─── Dot plot ─────────────────────────────────────────────────────────────────

def dot_plot(
    result:     pd.DataFrame,
    title:      str,
    outfile:    str,
    go_labels:  dict,
    top_n:      int = 20,
    fdr_cutoff: float = 0.05,
):
    sig = result[result['qvalue'] < fdr_cutoff].head(top_n).copy()
    if sig.empty:
        print(f"  No significant terms (FDR < {fdr_cutoff}) for: {title}")
        return

    sig = sig.sort_values('fold_enr', ascending=True)
    labels = []
    for t in sig['GO_term']:
        if t in go_labels:
            name, ns = go_labels[t]
            labels.append(f"{name} [{ns}]\n({t})")
        else:
            labels.append(t)

    fig, ax = plt.subplots(figsize=(9, max(4, len(sig) * 0.45 + 1.5)))
    sc = ax.scatter(
        sig['fold_enr'], range(len(sig)),
        c=sig['qvalue'], cmap='RdYlBu_r',
        s=sig['count'] * 12, vmin=0, vmax=fdr_cutoff,
        edgecolors='grey', linewidths=0.4, zorder=3,
    )
    ax.set_yticks(range(len(sig)))
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.axvline(1, color='grey', linestyle='--', linewidth=0.8)
    ax.set_xlabel('Fold enrichment', fontsize=10)
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.grid(axis='x', linestyle=':', linewidth=0.5, alpha=0.6)
    cb = plt.colorbar(sc, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label('FDR (q-value)', fontsize=9)

    # size legend
    for sz in [5, 10, 20]:
        ax.scatter([], [], s=sz * 12, c='grey', alpha=0.6, label=f'n={sz}')
    ax.legend(title='Gene count', fontsize=8, title_fontsize=8,
              loc='lower right', framealpha=0.7)

    plt.tight_layout()
    plt.savefig(outfile, dpi=180, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {outfile}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--emapper',        required=True, help='eggNOG .emapper.annotations file')
    ap.add_argument('--classification', required=True, help='20-13_classification.tsv (background)')
    ap.add_argument('--switching',      required=True, help='switching_genes.tsv')
    ap.add_argument('--obo',            default=None,  help='go-basic.obo for term labels (optional)')
    ap.add_argument('--outdir',         default='.',   help='output directory')
    ap.add_argument('--fdr',            type=float, default=0.05)
    ap.add_argument('--min_count',      type=int,   default=3,
                    help='min genes per GO term to test')
    args = ap.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)

    print("\n=== Loading data ===")
    go_map    = load_go_map(args.emapper)
    bg_genes  = load_background(args.classification, min_sites=15)
    sw        = load_switching_genes(args.switching)
    go_labels = try_load_go_labels(args.obo)

    print(f"  Switching genes total: {len(sw):,}")

    # ── define gene sets ──────────────────────────────────────────────────────
    gene_sets: dict[str, set[str]] = {}

    # 1. All shifters
    gene_sets['All_shifters'] = set(sw['gene_norm'])

    # 2. gbM <-> teM (both directions)
    gbm_tem = sw[sw['primary_transition'].isin(['gbM<->teM', 'teM<->gbM'])]
    gene_sets['gbM_teM_shifters'] = set(gbm_tem['gene_norm'])

    # 3. UM <-> teM
    um_tem = sw[sw['primary_transition'].isin(['UM<->teM', 'teM<->UM'])]
    gene_sets['UM_teM_shifters'] = set(um_tem['gene_norm'])

    # 4. Singleton changes only (is_singleton_change == True)
    singletons = sw[sw['is_singleton_change'] == True]
    gene_sets['Singleton_shifters'] = set(singletons['gene_norm'])

    # ── run enrichment for each set ───────────────────────────────────────────
    print("\n=== Running GO enrichment ===")
    all_results = {}
    for label, genes in gene_sets.items():
        # how many are actually in background?
        in_bg = genes & bg_genes
        print(f"\n  [{label}]  {len(genes):,} genes, {len(in_bg):,} in background")
        if len(in_bg) < 5:
            print("  -> skipping (too few genes in background)")
            continue

        result = go_enrichment(in_bg, bg_genes, go_map, min_count=args.min_count)
        if result.empty:
            print("  -> no terms passed minimum count filter")
            continue

        all_results[label] = result
        n_sig = (result['qvalue'] < args.fdr).sum()
        print(f"  -> {len(result):,} terms tested, {n_sig} significant at FDR < {args.fdr}")

        # save full table
        out_tsv = os.path.join(args.outdir, f'GO_{label}.tsv')
        if go_labels:
            result.insert(1, 'GO_name',      result['GO_term'].map(lambda t: go_labels.get(t, ('',))[0]))
            result.insert(2, 'GO_namespace', result['GO_term'].map(lambda t: go_labels.get(t, ('',''))[1]))
        result.to_csv(out_tsv, sep='\t', index=False)
        print(f"  Saved: {out_tsv}")

        # dot plot
        dot_plot(
            result, title=label.replace('_', ' '),
            outfile=os.path.join(args.outdir, f'GO_{label}_dotplot.png'),
            go_labels=go_labels, fdr_cutoff=args.fdr,
        )

    # ── summary across sets ───────────────────────────────────────────────────
    print("\n=== Summary ===")
    rows = []
    for label, result in all_results.items():
        sig = result[result['qvalue'] < args.fdr]
        rows.append({'gene_set': label,
                     'n_genes': len(gene_sets[label]),
                     'n_terms_tested': len(result),
                     'n_significant': len(sig),
                     'top_term': sig.iloc[0]['GO_term'] if not sig.empty else 'none',
                     'top_term_name': (go_labels.get(sig.iloc[0]['GO_term'], ('',))[0]
                                       if not sig.empty and go_labels else ''),
                     'top_qvalue': sig.iloc[0]['qvalue'] if not sig.empty else np.nan})
    summary = pd.DataFrame(rows)
    print(summary.to_string(index=False))
    summary.to_csv(os.path.join(args.outdir, 'GO_summary.tsv'), sep='\t', index=False)

    print("\nDone.")


if __name__ == '__main__':
    main()
