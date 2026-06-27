#!/usr/bin/env python3
"""
Script: 01_parse_modkit.py
Description:
    Parses raw modkit bedMethyl files. For each context (CG, CHG, CHH):
    - Extracts Nmod (methylated reads) and Nvalid (total valid coverage)
    - Uses parts[9] for Nvalid (true coverage), NOT parts[4] (score, can be capped at 1000)
    - CG:  Merges symmetrical +/- strand pairs. Unpaired sites are KEPT (ONT coverage fluctuates).
    - CHG: Merges symmetrical +/- strand pairs (1bp offset). Unpaired sites are kept.
    - CHH: No merging needed (asymmetric context).
    
    NOTE: Requires modkit to have been run with context flags (--cpg --chg --chh)
    so that context info (CG/CHG/CHH) is embedded in column 4 (mod_code).
    If your mod_code is just 'm' or '5mC', this script will warn and skip those lines.

    Output per context: chr <tab> start <tab> Nmod <tab> Nvalid
    
    Processes line-by-line for memory efficiency.

Usage:
    python 01_parse_modkit.py input.bed output_prefix
    
    Produces: output_prefix.CG.bed, output_prefix.CHG.bed, output_prefix.CHH.bed

Author: Paolo Callipo 
Date: 2026
"""

import sys
import argparse
from collections import deque


def parse_line(line_str):
    """Parse a bedMethyl line, return dict or None if malformed."""
    parts = line_str.strip().split('\t')
    if len(parts) < 17:
        return None
    try:
        # bedMethyl columns (0-indexed):
        #  0: chrom, 1: start, 2: end, 3: mod_code, 4: score (UNRELIABLE - capped at 1000)
        #  5: strand, 6: tStart, 7: tEnd, 8: color, 9: Nvalid_cov (TRUE coverage)
        # 10: percent_modified, 11: Nmod, 12: Ncanonical, 13: Nother_mod
        # 14: Ndelete, 15: Nfail, 16: Ndiff, 17: Nnocall
        return {
            'chrom':  parts[0],
            'start':  int(parts[1]),
            'end':    int(parts[2]),
            'context': extract_context(parts[3]),
            'nvalid': int(parts[9]),   # Nvalid_cov: true coverage, NOT parts[4]
            'strand': parts[5],
            'nmod':   int(parts[11]),
        }
    except (ValueError, IndexError):
        return None


def extract_context(mod_code):
    """
    Extract CG/CHG/CHH from the mod_code field.
    
    Handles two modkit output formats:
    - With --cpg/--chg/--chh flags: "m,CHH,0", "m,CG,0", "h,CHG,0"
    - Default modkit output: just "m" or "5mC" (context not embedded)
    
    If context is not in the mod_code, returns None.
    The caller should check and warn if too many None contexts are seen.
    """
    mod_code_upper = mod_code.upper()
    # Order matters: check CHH before CHG before CG to avoid partial matches
    if 'CHH' in mod_code_upper:
        return 'CHH'
    elif 'CHG' in mod_code_upper:
        return 'CHG'
    elif 'CG' in mod_code_upper:
        return 'CG'
    return None


def write_site(fout, chrom, start, nmod, nvalid):
    """Write a single processed site."""
    fout.write(f"{chrom}\t{start}\t{nmod}\t{nvalid}\n")


def process_file(input_file, output_prefix):
    """
    Single-pass processing of a modkit bedMethyl file.
    
    Strategy:
    - Read all lines, route by context
    - CHH: write directly (no merging)
    - CG/CHG: buffer the previous line per context.
      When a +/- pair is detected at the correct offset, merge and write.
      Otherwise, flush the buffer (skip for CG, write for CHG).
    """
    
    out_cg  = open(f"{output_prefix}.CG.bed", 'w')
    out_chg = open(f"{output_prefix}.CHG.bed", 'w')
    out_chh = open(f"{output_prefix}.CHH.bed", 'w')
    
    context_files = {'CG': out_cg, 'CHG': out_chg, 'CHH': out_chh}
    
    # Buffers for the previous + strand site per context (for merging)
    prev_plus = {'CG': None, 'CHG': None}
    
    counts = {'CG': 0, 'CHG': 0, 'CHH': 0, 'skipped': 0, 'no_context': 0}
    
    with open(input_file, 'r') as fin:
        for line_num, line in enumerate(fin, 1):
            if line.startswith('#'):
                continue
            
            rec = parse_line(line)
            if rec is None:
                counts['skipped'] += 1
                continue
            if rec['context'] is None:
                counts['no_context'] += 1
                if counts['no_context'] == 1:
                    print(f"  WARNING: Line {line_num} has no recognizable context in mod_code: "
                          f"'{line.strip().split(chr(9))[3]}'. If your modkit was run without "
                          f"--cpg/--chg/--chh, the context won't be in column 4 and this script "
                          f"cannot parse it. You may need to re-run modkit with context flags.",
                          file=sys.stderr)
                continue
            
            ctx = rec['context']
            fout = context_files[ctx]
            
            # --- CHH: no merging, write directly ---
            if ctx == 'CHH':
                write_site(fout, rec['chrom'], rec['start'], rec['nmod'], rec['nvalid'])
                counts['CHH'] += 1
                continue
            
            # --- CG / CHG: strand merging ---
            # Both CG and CHG: keep unpaired sites (ONT coverage fluctuates,
            # discarding unpaired CGs loses real data)
            if rec['strand'] == '+':
                # Flush any existing unmerged + site in buffer
                if prev_plus[ctx] is not None:
                    old = prev_plus[ctx]
                    # Keep unpaired + sites for BOTH CG and CHG
                    write_site(fout, old['chrom'], old['start'], old['nmod'], old['nvalid'])
                    counts[ctx] += 1
                
                # Buffer this + site
                prev_plus[ctx] = rec
            
            elif rec['strand'] == '-':
                prev = prev_plus[ctx]
                
                if prev is not None and prev['chrom'] == rec['chrom']:
                    # Check if this - site is the complement of the buffered + site
                    is_pair = False
                    
                    if ctx == 'CG':
                        # CG: + end == - start (adjacent, e.g., pos 100-101 and 101-102)
                        is_pair = (prev['end'] == rec['start'])
                    elif ctx == 'CHG':
                        # CHG: 1bp gap between + end and - start
                        is_pair = (prev['end'] + 1 == rec['start'])
                    
                    if is_pair:
                        # Merge: sum counts from both strands
                        merged_nmod = prev['nmod'] + rec['nmod']
                        merged_nvalid = prev['nvalid'] + rec['nvalid']
                        write_site(fout, prev['chrom'], prev['start'], merged_nmod, merged_nvalid)
                        counts[ctx] += 1
                        prev_plus[ctx] = None  # consumed
                    else:
                        # Not a pair: flush old + site, write this - site standalone
                        write_site(fout, prev['chrom'], prev['start'], prev['nmod'], prev['nvalid'])
                        counts[ctx] += 1
                        write_site(fout, rec['chrom'], rec['start'], rec['nmod'], rec['nvalid'])
                        counts[ctx] += 1
                        prev_plus[ctx] = None
                else:
                    # Different chromosome or no buffered + site
                    if prev is not None:
                        write_site(fout, prev['chrom'], prev['start'], prev['nmod'], prev['nvalid'])
                        counts[ctx] += 1
                    write_site(fout, rec['chrom'], rec['start'], rec['nmod'], rec['nvalid'])
                    counts[ctx] += 1
                    prev_plus[ctx] = None
            
            if line_num % 5_000_000 == 0:
                print(f"  Processed {line_num:,} lines...", file=sys.stderr)
    
    # Flush remaining buffers — keep unpaired sites for both CG and CHG
    for ctx in ['CG', 'CHG']:
        if prev_plus[ctx] is not None:
            old = prev_plus[ctx]
            write_site(context_files[ctx], old['chrom'], old['start'], old['nmod'], old['nvalid'])
            counts[ctx] += 1
    
    out_cg.close()
    out_chg.close()
    out_chh.close()
    
    # Report
    print(f"\nDone processing {input_file}", file=sys.stderr)
    print(f"  CG sites (merged+unpaired):  {counts['CG']:,}", file=sys.stderr)
    print(f"  CHG sites (merged+unpaired): {counts['CHG']:,}", file=sys.stderr)
    print(f"  CHH sites:                   {counts['CHH']:,}", file=sys.stderr)
    print(f"  Skipped (malformed):         {counts['skipped']:,}", file=sys.stderr)
    if counts['no_context'] > 0:
        print(f"  WARNING - No context found:  {counts['no_context']:,} lines", file=sys.stderr)
        print(f"    Your modkit may not have been run with --cpg/--chg/--chh flags.", file=sys.stderr)
    print(f"\nOutput files:", file=sys.stderr)
    print(f"  {output_prefix}.CG.bed", file=sys.stderr)
    print(f"  {output_prefix}.CHG.bed", file=sys.stderr)
    print(f"  {output_prefix}.CHH.bed", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Parse modkit bedMethyl: separate contexts, merge strands, output raw counts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_file", help="Raw modkit bedMethyl file (all contexts, both strands).")
    parser.add_argument("output_prefix", help="Output prefix. Produces <prefix>.CG.bed, .CHG.bed, .CHH.bed")
    
    args = parser.parse_args()
    process_file(args.input_file, args.output_prefix)


if __name__ == '__main__':
    main()
