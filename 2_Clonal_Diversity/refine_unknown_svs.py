#!/usr/bin/env python3
"""
refine_unknown_svs.py

Re-analyses ONLY the SVs currently classified as 'Unknown' by the existing
classify_sv.py pipeline. Existing TE-family and Centromeric/Satellite calls
are NOT touched. The Unknown bucket is decomposed into:

  1. Telomeric         : sequence contains (TTTAGGG){5,} or its reverse complement
  2. rDNA              : significant BLAST hit against V. vinifera rDNA reference
  3. Tandem repeat     : TRF detects a tandem repeat covering >=50% of the
                          sequence whose unit length is NOT in {79,107,135,187}
                          (i.e. not a centromeric satellite already covered)
  4. TE_diverged       : BLAST against EDTA TE library at >=60% identity over
                          >=60% query coverage (more permissive than the
                          original 80/80 threshold)
  5. Complex/Unique    : none of the above

Inputs (CONFIG block):
  - all_svs_with_lengths.tsv          : sv_id \\t SV_Type \\t length
  - insertions.fasta                   : insertion sequences (>chrom:pos)
  - deletions.bed                      : deletion coordinates
  - te_insertion_families.tsv          : insertions already classified as TEs
  - te_deletion_families.tsv           : deletions already classified as TEs
  - reference fasta + .fai             : for extracting deletion sequences
  - EDTA TE library FASTA              : for permissive re-BLAST
  - rDNA reference (auto-fetched)      : V. vinifera 45S+5S rDNA

Outputs:
  - unknown_refined_classification.tsv : per-SV refined category
  - unknown_refined_summary.tsv        : counts per category, by SV type
  - unknown_refined_donut.png          : updated donut chart for Figure 4D
  - work_dir/                          : intermediate FASTA, BLAST tables, TRF output

Tools required in PATH:
  - blastn, makeblastdb (for rDNA + TE re-search)
  - trf (Tandem Repeats Finder)
  - samtools (for deletion sequence extraction)
"""

import gzip
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

# ============================================================
# CONFIG
# ============================================================
ALL_SVS_TSV    = "all_svs_with_lengths.tsv"
INSERTIONS_FA  = "insertions.fasta"
DELETIONS_BED  = "deletions.bed"
TE_INS_TSV     = "te_insertion_families.tsv"
TE_DEL_TSV     = "te_deletion_families.tsv"

REFERENCE_FA   = "/mnt/data/paolo/reference/PN-20-13_DIPLOID_Final_masked.fasta"
EDTA_LIB       = "TElib.fa"

# Centromeric/satellite unit lengths used in the original classification
CENTRO_UNITS   = (79, 107, 135, 187)
CENTRO_TOLER   = 1
CENTRO_MAX_MULT = 20

# Permissive TE BLAST thresholds (less stringent than original 80/80)
TE_MIN_PIDENT  = 60.0
TE_MIN_QCOV    = 60.0

# tblastx pass for residual (catches diverged TEs at protein level).
# tblastx pident is amino-acid identity; we use lower thresholds because
# protein homology can be detected at much higher divergence levels.
TBX_MIN_PIDENT = 35.0
TBX_MIN_QCOV   = 30.0

# Organellar BLAST thresholds (chloroplast and mitochondrion)
# Plant NUPTs/NUMTs typically retain high identity to the source organellar
# genome, so we use moderately stringent thresholds.
ORG_MIN_PIDENT = 85.0
ORG_MIN_QCOV   = 50.0

# Self-match (segmental duplication) thresholds.
# Strict identity ensures we only call true segmental duplications, not
# distant homology. We also require the match to be located at a different
# genomic locus (>10 kb away from the SV's own coordinates) to exclude
# the SV's reference position itself.
SELF_MIN_PIDENT = 90.0
SELF_MIN_QCOV   = 50.0
SELF_MIN_DIST   = 10000

# rDNA pass thresholds (existing)
RDNA_MIN_PIDENT = 75.0
RDNA_MIN_QCOV   = 30.0

# TRF parameters: match mismatch indel P_match P_indel min_score max_period
TRF_PARAMS = ("2", "7", "7", "80", "10", "30", "500")
TRF_MIN_COV = 0.3

# Short-fragment threshold
SHORT_FRAGMENT_BP = 100

# Reference paths for organellar genomes (you must provide these)
CHLOROPLAST_FA = "vitis_chloroplast.fa"   # e.g., NCBI NC_007957.1
MITOCHONDRION_FA = "vitis_mitochondrion.fa" # e.g., NCBI NC_012119.1

# Telomere repeat
TELOMERE_FWD = re.compile(r"(?:TTTAGGG){5,}", re.IGNORECASE)
TELOMERE_REV = re.compile(r"(?:CCCTAAA){5,}", re.IGNORECASE)

# Path to user-provided rDNA reference. If missing, the script will print
# instructions for fetching one and exit cleanly so the user can rerun.
RDNA_REFERENCE = "vitis_rdna.fa"

WORK_DIR = "unknown_refine_work"

OUT_PER_SV    = "unknown_refined_classification.tsv"
OUT_SUMMARY   = "unknown_refined_summary.tsv"
OUT_PLOT      = "unknown_refined_donut.png"
# ============================================================


def run(cmd, **kw):
    """Run shell command, raise on error."""
    print(f"[run] {' '.join(cmd) if isinstance(cmd, list) else cmd}", file=sys.stderr)
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def ensure_tool(name):
    if shutil.which(name) is None:
        sys.exit(f"ERROR: '{name}' not found in PATH. Please load/install it.")


def generate_centro_sizes():
    sizes = set()
    for u in CENTRO_UNITS:
        for i in range(1, CENTRO_MAX_MULT + 1):
            for t in range(-CENTRO_TOLER, CENTRO_TOLER + 1):
                sizes.add(u * i + t)
    return sizes


def load_te_ids(path):
    """Load SV IDs that are already classified as TEs (column 0)."""
    ids = set()
    with open(path) as fh:
        first = fh.readline()
        if "Insertion_ID" not in first and "Deletion_ID" not in first:
            fh.seek(0)
        for line in fh:
            parts = line.strip().split("\t")
            if parts and parts[0]:
                ids.add(parts[0])
    return ids


def identify_unknown_svs():
    """Replicate classify_sv.py's logic to identify SVs classified as Unknown.
    Returns dict: {sv_id -> ("INS"|"DEL", length)}."""
    centro = generate_centro_sizes()
    te_ins = load_te_ids(TE_INS_TSV)
    te_del = load_te_ids(TE_DEL_TSV)
    unknown = {}
    with open(ALL_SVS_TSV) as fh:
        for line in fh:
            sv_id, sv_type, sv_len = line.rstrip("\n").split("\t")
            sv_len_i = abs(int(sv_len))
            if sv_type == "INS":
                if sv_id in te_ins or sv_len_i in centro:
                    continue
                unknown[sv_id] = ("INS", sv_len_i)
            elif sv_type == "DEL":
                if sv_id in te_del or sv_len_i in centro:
                    continue
                unknown[sv_id] = ("DEL", sv_len_i)
    n_ins = sum(1 for v in unknown.values() if v[0] == "INS")
    n_del = sum(1 for v in unknown.values() if v[0] == "DEL")
    print(f"[identify] Unknown insertions: {n_ins:,}", file=sys.stderr)
    print(f"[identify] Unknown deletions : {n_del:,}", file=sys.stderr)
    return unknown


def load_insertion_sequences():
    """Return {sv_id: sequence} from insertions.fasta."""
    seqs = {}
    cur_id = None
    cur_seq = []
    with open(INSERTIONS_FA) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(cur_seq)
    return seqs


def extract_deletion_sequences(unknown_dels, fasta, work_dir):
    """Use samtools faidx to extract reference sequence at each deletion coord.
    
    SV IDs in all_svs_with_lengths.tsv use Sniffles2's 1-based POS:
        e.g. 'PN1_HapA:765787' for a deletion starting at position 765787 (1-based).
    deletions.bed uses 0-based BED coordinates:
        e.g. 'PN1_HapA  765786  769830  4043' for the same event.
    So the SV-ID matching key is f"{chrom}:{bed_start + 1}".
    """
    bed_lookup = {}
    with open(DELETIONS_BED) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            sv_id = f"{chrom}:{start + 1}"   # match the 1-based ID convention
            bed_lookup[sv_id] = (chrom, start, end)

    # Write a region list for samtools (one region per line)
    regions_path = os.path.join(work_dir, "del_regions.txt")
    out_fa_path = os.path.join(work_dir, "deletions.fasta")
    matched = 0
    with open(regions_path, "w") as fh:
        for sv_id in unknown_dels:
            if sv_id in bed_lookup:
                c, s, e = bed_lookup[sv_id]
                # samtools faidx is 1-based inclusive
                fh.write(f"{c}:{s+1}-{e}\n")
                matched += 1
    print(f"[extract] Matched {matched:,} of {len(unknown_dels):,} unknown deletions "
          f"to BED entries", file=sys.stderr)

    # Run samtools
    print("[extract] Running samtools faidx for deletion sequences...", file=sys.stderr)
    with open(out_fa_path, "w") as out:
        subprocess.run(
            ["samtools", "faidx", "-r", regions_path, fasta],
            stdout=out, check=True
        )

    # Parse the output, keying by chrom:pos (1-based start to match our IDs)
    seqs = {}
    cur_header = None
    cur_seq = []
    with open(out_fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_header:
                    # samtools header: ">chrom:start-end" (1-based)
                    m = re.match(r">([^:]+):(\d+)-\d+", cur_header)
                    if m:
                        sv_id = f"{m.group(1)}:{m.group(2)}"
                        seqs[sv_id] = "".join(cur_seq)
                cur_header = line.strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_header:
            m = re.match(r">([^:]+):(\d+)-\d+", cur_header)
            if m:
                sv_id = f"{m.group(1)}:{m.group(2)}"
                seqs[sv_id] = "".join(cur_seq)
    print(f"[extract] Extracted {len(seqs):,} deletion sequences", file=sys.stderr)
    return seqs


def is_telomeric(seq):
    """Return True if sequence contains telomeric repeat ≥5 copies."""
    if TELOMERE_FWD.search(seq) or TELOMERE_REV.search(seq):
        return True
    return False


def fetch_rdna_reference(work_dir):
    """Verify that a user-supplied rDNA reference exists and looks like a real
    plant rRNA sequence. If not, print clear instructions and exit cleanly."""
    rdna_path = os.path.join(work_dir, "vitis_rdna.fa")
    # Allow the user to drop a vitis_rdna.fa anywhere they like
    if os.path.exists(RDNA_REFERENCE):
        if not os.path.exists(rdna_path) or os.path.realpath(RDNA_REFERENCE) != os.path.realpath(rdna_path):
            shutil.copy(RDNA_REFERENCE, rdna_path)
    if not os.path.exists(rdna_path):
        msg = f"""
ERROR: rDNA reference not found.
The script needs a FASTA file of plant rDNA (18S, 5.8S, 26S/28S) to identify
SVs derived from ribosomal DNA. Please obtain one and place it at:

    {os.path.abspath(rdna_path)}
or
    {os.path.abspath(RDNA_REFERENCE)}

Recommended sources:

1. NCBI web interface:
   - Search 'Vitis vinifera 18S ribosomal RNA' at https://www.ncbi.nlm.nih.gov/nuccore/
   - Filter sequence length 1500-2000, click any hit, download as FASTA.
   - Repeat for '26S' (filter 3000-3500) and '5.8S internal transcribed spacer'.
   - Concatenate the three FASTAs into vitis_rdna.fa.

2. Entrez Direct (if installed):
   esearch -db nuccore -query 'Vitis vinifera[ORGN] AND 18S ribosomal RNA[Title] AND 1500:2000[SLEN]' \\
     | efetch -format fasta > vitis_18S.fa
   esearch -db nuccore -query 'Vitis vinifera[ORGN] AND 26S ribosomal RNA[Title] AND 3000:3500[SLEN]' \\
     | efetch -format fasta > vitis_26S.fa
   cat vitis_18S.fa vitis_26S.fa > vitis_rdna.fa

3. Use Arabidopsis rDNA as a fallback (>=95% conserved with Vitis):
   Download from NCBI:
     - 18S: X16077 (~1808 bp)
     - 25S: X52320 (~3376 bp; this entry contains 5.8S + ITS as well)
   esearch -db nuccore -query 'X16077.1' | efetch -format fasta >  arabidopsis_rdna.fa
   esearch -db nuccore -query 'X52320.1' | efetch -format fasta >> arabidopsis_rdna.fa
   cp arabidopsis_rdna.fa vitis_rdna.fa
"""
        sys.exit(msg)

    # Sanity check: the file should be at least 1500 bp and look like DNA
    n_bases = 0
    headers = 0
    with open(rdna_path) as fh:
        for line in fh:
            if line.startswith(">"):
                headers += 1
            else:
                n_bases += len(line.strip())
    if n_bases < 1500 or headers == 0:
        sys.exit(f"ERROR: {rdna_path} looks too small ({n_bases} bp, {headers} headers). "
                 "Please provide a proper plant rDNA FASTA.")
    print(f"[rdna] Using {rdna_path} ({headers} sequences, {n_bases:,} bp total)",
          file=sys.stderr)
    return rdna_path


def make_blast_db(fasta, db_prefix):
    if not (os.path.exists(db_prefix + ".nhr") or os.path.exists(db_prefix + ".ndb")):
        run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", db_prefix])


def run_blastn(query_fa, db_prefix, out_tsv, evalue="1e-5"):
    cmd = [
        "blastn",
        "-query", query_fa,
        "-db", db_prefix,
        "-out", out_tsv,
        "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
        "-evalue", evalue,
        "-max_target_seqs", "1",
        "-num_threads", "4",
    ]
    run(cmd)


def parse_blast_best(tsv_path):
    """Return {qseqid: (pident, length, qlen, qcov, evalue)} keeping the best
    hit per query (max bitscore)."""
    best = {}
    if not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
        return best
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            qid = parts[0]
            try:
                pid = float(parts[2])
                length = int(parts[3])
                qlen = int(parts[4])
                evalue = float(parts[6])
                bitscore = float(parts[7])
            except ValueError:
                continue
            qcov = 100 * length / qlen if qlen else 0
            entry = (pid, length, qlen, qcov, evalue, bitscore)
            if qid not in best or entry[5] > best[qid][5]:
                best[qid] = entry
    # Strip bitscore from the returned tuples
    return {q: v[:5] for q, v in best.items()}


def run_trf(query_fa, work_dir):
    """Run TRF and return {qseqid: max_array_coverage_fraction, max_unit_len}."""
    trf_outdir = os.path.join(work_dir, "trf")
    os.makedirs(trf_outdir, exist_ok=True)
    # Copy the FASTA inside trf_outdir so TRF writes its outputs there
    local_fa = os.path.join(trf_outdir, os.path.basename(query_fa))
    shutil.copy(query_fa, local_fa)
    cmd = ["trf", os.path.basename(query_fa), *TRF_PARAMS, "-d", "-h", "-ngs"]
    print(f"[trf] Running: cd {trf_outdir} && {' '.join(cmd)}", file=sys.stderr)
    # TRF with -ngs writes to stdout. Capture it.
    proc = subprocess.run(cmd, cwd=trf_outdir, check=False,
                          capture_output=True, text=True)
    if proc.returncode not in (0, 1):
        # TRF often returns nonzero exit codes even on success; only fail on
        # genuinely missing output
        print(f"[trf] returncode={proc.returncode}, but continuing.", file=sys.stderr)

    # Parse the -ngs output: a multi-record format starting with "@<seqid>",
    # then one repeat per line (start end period_size copies ... %match %indel
    # score %A %C %G %T entropy consensus repeat_seq).
    arrays = defaultdict(list)
    cur = None
    for line in proc.stdout.splitlines():
        line = line.rstrip()
        if line.startswith("@"):
            cur = line[1:].split()[0]
            continue
        if cur is None or not line:
            continue
        parts = line.split()
        if len(parts) < 14:
            continue
        try:
            start = int(parts[0]); end = int(parts[1])
            period = int(parts[2])
        except ValueError:
            continue
        arrays[cur].append((start, end, period))

    # Compute max array coverage and corresponding unit length per query
    out = {}
    qlens = {}
    # parse the local FASTA to get query lengths
    cur_id = None
    cur_seq = []
    with open(local_fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None:
                    qlens[cur_id] = len("".join(cur_seq))
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id is not None:
            qlens[cur_id] = len("".join(cur_seq))

    for qid, runs in arrays.items():
        # Best array = longest single span
        best_span = 0
        best_period = 0
        for s, e, p in runs:
            span = e - s + 1
            if span > best_span:
                best_span = span
                best_period = p
        qlen = qlens.get(qid, 0)
        cov = best_span / qlen if qlen else 0
        out[qid] = (cov, best_period)
    return out


def write_fasta(seq_dict, path):
    with open(path, "w") as fh:
        for k, v in seq_dict.items():
            fh.write(f">{k}\n{v}\n")


def main():
    # Tool checks
    for tool in ("blastn", "makeblastdb", "trf", "samtools"):
        ensure_tool(tool)

    os.makedirs(WORK_DIR, exist_ok=True)

    # 1. Identify unknown SVs
    unknown = identify_unknown_svs()

    # 2. Load/extract sequences
    print("[seq] Loading insertion sequences...", file=sys.stderr)
    ins_seqs_all = load_insertion_sequences()
    unknown_ins = {sid: ins_seqs_all[sid]
                   for sid, (t, _) in unknown.items()
                   if t == "INS" and sid in ins_seqs_all}
    missing_ins = sum(1 for sid, (t, _) in unknown.items()
                      if t == "INS" and sid not in ins_seqs_all)
    print(f"[seq] Loaded {len(unknown_ins):,} unknown insertion seqs "
          f"({missing_ins} missing in FASTA)", file=sys.stderr)

    unknown_del_ids = [sid for sid, (t, _) in unknown.items() if t == "DEL"]
    unknown_dels = extract_deletion_sequences(unknown_del_ids, REFERENCE_FA, WORK_DIR)

    all_unknown = {**unknown_ins, **unknown_dels}
    print(f"[seq] Total unknown sequences: {len(all_unknown):,}", file=sys.stderr)

    # Write a single combined FASTA
    combined_fa = os.path.join(WORK_DIR, "unknown_combined.fasta")
    write_fasta(all_unknown, combined_fa)

    # 3. Build / load BLAST DBs
    rdna_fa = fetch_rdna_reference(WORK_DIR)
    rdna_db = os.path.join(WORK_DIR, "rdna_db")
    make_blast_db(rdna_fa, rdna_db)

    edta_db = os.path.join(WORK_DIR, "edta_db")
    make_blast_db(EDTA_LIB, edta_db)

    # 4. Run BLAST and TRF
    rdna_blast_out = os.path.join(WORK_DIR, "rdna_blast.tsv")
    run_blastn(combined_fa, rdna_db, rdna_blast_out, evalue="1e-10")
    rdna_hits = parse_blast_best(rdna_blast_out)
    rdna_pass = {q for q, (pid, _, _, qcov, _) in rdna_hits.items()
                 if pid >= RDNA_MIN_PIDENT and qcov >= RDNA_MIN_QCOV}
    print(f"[rdna] {len(rdna_pass):,} sequences with rDNA hit", file=sys.stderr)

    # 4b. Organellar passes (chloroplast and mitochondrion)
    organellar_pass = {}   # sv_id -> "Chloroplast" or "Mitochondrion"
    for fa_path, label, db_suffix in [
        (CHLOROPLAST_FA, "Chloroplast", "cp_db"),
        (MITOCHONDRION_FA, "Mitochondrion", "mt_db"),
    ]:
        if not os.path.exists(fa_path):
            print(f"[org] {label} reference not found at {fa_path} — skipping. "
                  f"To enable this category, place a Vitis vinifera {label.lower()} "
                  f"FASTA at this path.", file=sys.stderr)
            continue
        db_prefix = os.path.join(WORK_DIR, db_suffix)
        make_blast_db(fa_path, db_prefix)
        out_tsv = os.path.join(WORK_DIR, f"{db_suffix}_blast.tsv")
        run_blastn(combined_fa, db_prefix, out_tsv, evalue="1e-10")
        hits = parse_blast_best(out_tsv)
        n_pass = 0
        for q, (pid, _, _, qcov, _) in hits.items():
            if pid >= ORG_MIN_PIDENT and qcov >= ORG_MIN_QCOV:
                # Only assign if not already assigned (chloroplast checked first)
                if q not in organellar_pass:
                    organellar_pass[q] = label
                    n_pass += 1
        print(f"[org] {n_pass:,} sequences with {label} hit "
              f"({ORG_MIN_PIDENT}/{ORG_MIN_QCOV})", file=sys.stderr)

    # 4c. Self-similarity pass: BLAST the unknown sequences against the
    # diploid reference itself to detect segmental duplications. We exclude
    # hits falling at the SV's own coordinates (within +/- SELF_MIN_DIST bp).
    self_pass = set()
    self_db_prefix = os.path.join(WORK_DIR, "self_db")
    if not (os.path.exists(self_db_prefix + ".nhr") or
            os.path.exists(self_db_prefix + ".ndb")):
        run(["makeblastdb", "-in", REFERENCE_FA, "-dbtype", "nucl",
             "-out", self_db_prefix])
    self_blast_out = os.path.join(WORK_DIR, "self_blast.tsv")
    if not os.path.exists(self_blast_out) or os.path.getsize(self_blast_out) == 0:
        cmd = ["blastn",
               "-query", combined_fa, "-db", self_db_prefix, "-out", self_blast_out,
               "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore sstart send",
               "-evalue", "1e-20", "-max_target_seqs", "5", "-num_threads", "4",
               "-perc_identity", "80"]
        run(cmd)
    # Parse hits, keeping any hit that meets thresholds AND is far from the SV's own locus
    with open(self_blast_out) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            qid = parts[0]
            sid = parts[1]
            try:
                pid = float(parts[2])
                length = int(parts[3])
                qlen = int(parts[4])
                sstart = int(parts[8])
                send = int(parts[9])
            except ValueError:
                continue
            qcov = 100 * length / qlen if qlen else 0
            if pid < SELF_MIN_PIDENT or qcov < SELF_MIN_QCOV:
                continue
            # Parse SV's own coordinates from the SV ID (chrom:pos)
            m = re.match(r"([^:]+):(\d+)$", qid)
            if not m:
                continue
            sv_chrom = m.group(1)
            sv_pos = int(m.group(2))
            hit_chrom = sid
            hit_start = min(sstart, send)
            hit_end = max(sstart, send)
            # Exclude self-locus: hit on same chrom within +/- SELF_MIN_DIST
            if hit_chrom == sv_chrom:
                if abs(hit_start - sv_pos) < SELF_MIN_DIST and \
                   abs(hit_end - sv_pos) < SELF_MIN_DIST:
                    continue   # this is the SV's own reference locus
            self_pass.add(qid)
    print(f"[self] {len(self_pass):,} sequences with self-match elsewhere "
          f"in genome (segmental duplications)", file=sys.stderr)

    edta_blast_out = os.path.join(WORK_DIR, "edta_permissive_blast.tsv")
    run_blastn(combined_fa, edta_db, edta_blast_out, evalue="1e-5")
    te_hits = parse_blast_best(edta_blast_out)
    te_pass = {q for q, (pid, _, _, qcov, _) in te_hits.items()
               if pid >= TE_MIN_PIDENT and qcov >= TE_MIN_QCOV}
    print(f"[te_div] {len(te_pass):,} sequences with diverged-TE hit "
          f"({TE_MIN_PIDENT}/{TE_MIN_QCOV})", file=sys.stderr)

    # tblastx pass: catches highly diverged TEs via protein-level homology
    tbx_out = os.path.join(WORK_DIR, "edta_tblastx.tsv")
    if not os.path.exists(tbx_out) or os.path.getsize(tbx_out) == 0:
        cmd = ["tblastx",
               "-query", combined_fa, "-db", edta_db, "-out", tbx_out,
               "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
               "-evalue", "1e-3", "-max_target_seqs", "1", "-num_threads", "4"]
        run(cmd)
    tbx_hits = parse_blast_best(tbx_out)
    tbx_pass = set()
    for q, (pid, length, qlen, qcov_orig, _) in tbx_hits.items():
        if qlen == 0:
            continue
        qcov_nt = 100 * 3 * length / qlen
        if pid >= TBX_MIN_PIDENT and qcov_nt >= TBX_MIN_QCOV:
            tbx_pass.add(q)
    print(f"[tbx_te] {len(tbx_pass):,} sequences with diverged-TE protein hit "
          f"({TBX_MIN_PIDENT}/{TBX_MIN_QCOV})", file=sys.stderr)

    trf_results = run_trf(combined_fa, WORK_DIR)
    centro_units_set = set(CENTRO_UNITS)
    trf_pass = {q for q, (cov, period) in trf_results.items()
                if cov >= TRF_MIN_COV and period not in centro_units_set}
    print(f"[trf] {len(trf_pass):,} sequences with non-centromeric tandem array "
          f"covering ≥{int(TRF_MIN_COV*100)}% of length", file=sys.stderr)

    # 5. Apply priority classification
    # Priority order (most specific first):
    #   Telomeric → rDNA → Chloroplast → Mitochondrion → Tandem repeat
    #   → TE diverged (blastn or tblastx) → Segmental duplication
    #   → Short fragment → Complex/unique
    rows = []
    counts = defaultdict(lambda: defaultdict(int))
    for sv_id, (sv_type, sv_len) in unknown.items():
        seq = all_unknown.get(sv_id, "")
        if not seq:
            cat = "no_sequence_available"
        elif is_telomeric(seq):
            cat = "Telomeric"
        elif sv_id in rdna_pass:
            cat = "rDNA"
        elif sv_id in organellar_pass:
            cat = organellar_pass[sv_id]   # "Chloroplast" or "Mitochondrion"
        elif sv_id in trf_pass:
            cat = "Tandem_repeat_other"
        elif sv_id in te_pass or sv_id in tbx_pass:
            cat = "TE_diverged"
        elif sv_id in self_pass:
            cat = "Segmental_duplication"
        elif sv_len < SHORT_FRAGMENT_BP:
            cat = "Short_fragment"
        else:
            cat = "Complex_unique"
        rows.append((sv_id, sv_type, sv_len, cat, len(seq)))
        counts[sv_type][cat] += 1

    # 6. Write per-SV TSV
    with open(OUT_PER_SV, "w") as fh:
        fh.write("SV_ID\tSV_Type\tLength\tRefined_Category\tSequence_Length\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # 7. Summary table
    cats_in_priority = ["Telomeric", "rDNA", "Chloroplast", "Mitochondrion",
                        "Tandem_repeat_other", "TE_diverged",
                        "Segmental_duplication", "Short_fragment",
                        "Complex_unique", "no_sequence_available"]
    with open(OUT_SUMMARY, "w") as fh:
        fh.write("SV_Type\tCategory\tCount\tPct_of_Unknown\n")
        for t in ("INS", "DEL"):
            tot = sum(counts[t].values())
            for c in cats_in_priority:
                n = counts[t].get(c, 0)
                pct = 100 * n / tot if tot else 0
                fh.write(f"{t}\t{c}\t{n}\t{pct:.2f}\n")

    print("\n=== REFINED UNKNOWN SV BREAKDOWN ===\n")
    for t in ("INS", "DEL"):
        tot = sum(counts[t].values())
        print(f"--- {t} (total Unknown: {tot:,}) ---")
        for c in cats_in_priority:
            n = counts[t].get(c, 0)
            pct = 100 * n / tot if tot else 0
            print(f"  {c:25s}: {n:>5,} ({pct:5.2f}%)")
        print()

    # 8. Donut chart
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        colors = {
            "Telomeric":             "#3498db",
            "rDNA":                  "#9b59b6",
            "Chloroplast":           "#27ae60",
            "Mitochondrion":         "#c0392b",
            "Tandem_repeat_other":   "#e67e22",
            "TE_diverged":           "#16a085",
            "Segmental_duplication": "#f39c12",
            "Short_fragment":        "#bdc3c7",
            "Complex_unique":        "#7f8c8d",
        }
        for ax, t, label in zip(axes, ("INS", "DEL"),
                                 ("Insertions", "Deletions")):
            tot = sum(counts[t].values())
            cats = [c for c in cats_in_priority if c != "no_sequence_available"]
            sizes = [counts[t].get(c, 0) for c in cats]
            cs = [colors[c] for c in cats]
            ax.pie(sizes, labels=cats, colors=cs, autopct="%1.1f%%",
                   pctdistance=0.78, wedgeprops=dict(width=0.4, edgecolor="white"))
            ax.set_title(f"{label}\n(Unknown SVs, n = {tot:,})")
        fig.suptitle("Refined classification of previously 'Unknown' SVs",
                     y=1.02, fontsize=12)
        fig.tight_layout()
        fig.savefig(OUT_PLOT, dpi=200, bbox_inches="tight")
        print(f"[plot] Wrote {OUT_PLOT}", file=sys.stderr)
    except Exception as e:
        print(f"[plot] Could not generate donut chart: {e}", file=sys.stderr)

    print(f"\nWrote: {OUT_PER_SV}", file=sys.stderr)
    print(f"Wrote: {OUT_SUMMARY}", file=sys.stderr)


if __name__ == "__main__":
    main()
