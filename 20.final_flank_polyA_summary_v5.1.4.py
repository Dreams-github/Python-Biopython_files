#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
flank_polyA_summary_v5.1.5.py  (CSV Primary/Secondary + ' no flank' + wasRC)

- Pulls Primary/Secondary from CSV _15692-SYB.csv (no header):
    * Primary  -> Column D (index 3)
    * Secondary-> Column E (index 4)
- Matches FASTQ barcode _S###_ to CSV Column A's 4-digit code (e.g., S101 -> '0101').
- Adds samplePrimary, sampleSecondary, wasRC to the summary.
- Strict A-only polyA; no flank strings for missing flank sequences; zeros for metrics when absent.
"""

import os
import re
import csv
import gzip
from statistics import mean
from typing import Optional, Tuple, List, Dict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align


# ====================== USER SETTINGS =========================

# Placeholder when a flank sequence is absent
NA_STR = " no flank "

RUN_DIR       = r"provide the path"
FASTQ_SUBDIR  = "Fastq files"

# Reference CSV (no header) sitting inside "Fastq files"
REF_CSV_BASENAME = "_14417-SYB.csv"

# CSV column indices (0-based) per your screenshot
BARCODE_KEY_COL       = 0  # Column A, e.g., '0101'
SAMPLE_PRIMARY_COL    = 3  # Column D, e.g., 'Moderna2Large'
SAMPLE_SECONDARY_COL  = 4  # Column E, e.g., 'sq104411_1'

# Vector sequence column (used to extract upstream/downstream reference flanks)
REF_VECTOR_COL        = 8  # Column I

# PolyA detection (fringe-friendly)
WINDOW  = 40
THRESH  = 0.90
MIN_LEN = 20

# Reference flank extraction around polyA
REF_FLANK_LEN = 50

# Aligner thresholds
USE_ALIGNER = True
REF_MIN_IDENTITY      = 0.90
REF_MIN_COVER_FRAC    = 0.70
REF_MIN_SCORE         = 12.0
ALIGN_MATCH_SCORE     = 1.0
ALIGN_MISMATCH_SCORE  = -0.5
ALIGN_OPEN_GAP_SCORE  = -1.0
ALIGN_EXTEND_GAP_SCORE= -0.5

# Output directory
OUTPUT_DIR = os.path.join(RUN_DIR, "flank_output_14417")

# Fallback names if CSV row is missing
UNKNOWN_PRIMARY   = "UNKNOWN_PRIMARY"
UNKNOWN_SECONDARY = "UNKNOWN_SECONDARY"

# ====================== HELPERS ===============================

S_KEY_RE   = re.compile(r"_S(\d{1,4})_")  # ..._S101_...

def reverse_complement(s: str) -> str:
    return str(Seq(s).reverse_complement())

def rc_quals(q: List[int]) -> List[int]:
    return list(reversed(q))

def clamp_span(n: int, start: int, end: int) -> Tuple[int, int]:
    s = max(0, start); e = min(n, end)
    return (s, e) if e >= s else (s, s)

def mean_q(quals: List[int], start: int, end: int) -> Optional[float]:
    s, e = clamp_span(len(quals), start, end)
    if e <= s: return None
    return float(mean(quals[s:e]))

def infer_run_name(file_basename: str) -> str:
    """Return FASTQ base name without .fastq[.gz]."""
    return re.sub(r'\.fastq(\.gz)?$', '', file_basename, flags=re.IGNORECASE)

def longest_A_run_coords(s: str) -> Optional[Tuple[int,int]]:
    best_len=cur=0; best_s=best_e=None; start=None
    for i,ch in enumerate(s):
        if ch=='A':
            if cur==0: start=i
            cur+=1
            if cur>best_len: best_len=cur; best_s=start; best_e=i+1
        else:
            cur=0; start=None
    return (best_s,best_e) if best_s is not None else None

def extract_known_ref_flanks(ref: str, flank_len: int = 50) -> Tuple[str,str,Tuple[int,int]]:
    coords = longest_A_run_coords(ref)
    if coords is None:
        raise ValueError("No A-run found in reference; cannot derive known flanks.")
    pa_s, pa_e = coords
    L = len(ref)
    up_s, up_e = max(0, pa_s - flank_len), pa_s
    dn_s, dn_e = pa_e, min(L, pa_e + flank_len)
    upstream = ref[up_s:up_e]
    downstream = ref[dn_s:dn_e]
    return upstream, downstream, (pa_s, pa_e)

def sliding_poly_runs(seq: str, base: str, window: int, thresh: float) -> List[Tuple[int,int]]:
    n=len(seq)
    if n<window: return []
    is_base=[1 if b==base else 0 for b in seq]
    wsum=sum(is_base[:window]); hits=[]
    if wsum/window>=thresh: hits.append((0,window))
    for i in range(1,n-window+1):
        wsum += is_base[i+window-1] - is_base[i-1]
        if wsum/window>=thresh: hits.append((i,i+window))
    return hits

def merge_intervals(iv):
    if not iv: return []
    iv.sort(); merged=[iv[0]]
    for s,e in iv[1:]:
        ls,le=merged[-1]
        if s<=le: merged[-1]=(ls,max(le,e))
        else: merged.append((s,e))
    return merged

def region_purity(seq,start,end,base):
    if end<=start: return 0.0
    seg=seq[start:end]; return seg.count(base)/(end-start)

def all_poly_runs(seq: str) -> List[dict]:
    a = sliding_poly_runs(seq,'A',WINDOW,THRESH)
    t = sliding_poly_runs(seq,'T',WINDOW,THRESH)
    a_runs=[(s,e) for s,e in merge_intervals(a) if e-s>=MIN_LEN]
    t_runs=[(s,e) for s,e in merge_intervals(t) if e-s>=MIN_LEN]
    res=[]
    for s,e in a_runs:
        res.append({"kind":"polyA","start":s,"end":e,"length":e-s,"purity":region_purity(seq,s,e,"A")})
    for s,e in t_runs:
        res.append({"kind":"polyT","start":s,"end":e,"length":e-s,"purity":region_purity(seq,s,e,"T")})
    res.sort(key=lambda r:(r["start"], -r["length"]))
    return res

def init_aligner():
    al = Align.PairwiseAligner()
    al.mode = "local"
    al.match_score     = ALIGN_MATCH_SCORE
    al.mismatch_score  = ALIGN_MISMATCH_SCORE
    al.open_gap_score  = ALIGN_OPEN_GAP_SCORE
    al.extend_gap_score= ALIGN_EXTEND_GAP_SCORE
    return al

def align_flank_to_ref(flank: str, ref_flank: str) -> Tuple[float,float,float]:
    if not flank: return (0.0,0.0,0.0)
    al = init_aligner()
    alns = al.align(ref_flank, flank)
    if alns is None or len(alns)==0: return (0.0,0.0,0.0)
    aln = alns[0]
    try:
        t_blocks, q_blocks = aln.aligned  # may be numpy arrays
    except Exception:
        return (0.0, float(aln.score), 0.0)
    if hasattr(t_blocks,'tolist'): t_blocks=[(int(a),int(b)) for (a,b) in t_blocks.tolist()]
    else:                           t_blocks=[(int(a),int(b)) for (a,b) in (t_blocks or [])]
    if hasattr(q_blocks,'tolist'): q_blocks=[(int(a),int(b)) for (a,b) in q_blocks.tolist()]
    else:                           q_blocks=[(int(a),int(b)) for (a,b) in (q_blocks or [])]
    if len(t_blocks)==0 or len(q_blocks)==0: return (0.0, float(aln.score), 0.0)
    matches=0; aligned_q_len=0
    for (ts,te),(qs,qe) in zip(t_blocks,q_blocks):
        seg_t = ref_flank[ts:te]; seg_q = flank[qs:qe]
        matches       += sum(1 for a,b in zip(seg_t,seg_q) if a==b)
        aligned_q_len += (qe-qs)
    identity = 0.0 if aligned_q_len==0 else (matches/aligned_q_len)
    cover    = 0.0 if len(flank)==0 else (aligned_q_len/len(flank))
    return (identity, float(aln.score), cover)

def load_ref_sequence(csv_path: str, col_index: int) -> str:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(csv_path)
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as f:
        r = csv.reader(f)
        for row in r:
            if col_index < len(row):
                s = (row[col_index] or "").strip().upper().replace(" ", "")
                s = "".join(ch for ch in s if ch in "ACGTN")
                if s: return s
    raise ValueError(f"No usable sequence found in column index {col_index}")

def build_primary_secondary_map(csv_path: str,
                                key_col: int,
                                p_col: int,
                                s_col: int) -> Dict[str, Tuple[str,str]]:
    """
    Returns mapping: '0101' -> (Primary, Secondary)
    Uses the FIRST occurrence for each key (as Jeff requested).
    """
    mapping = {}
    with open(csv_path, "r", newline="", encoding="utf-8-sig") as f:
        r = csv.reader(f)
        for row in r:
            if key_col >= len(row): continue
            key = (row[key_col] or "").strip()
            if not key: continue
            if key in mapping:  # keep first occurrence only
                continue
            prim = (row[p_col] if p_col < len(row) else "").strip() or UNKNOWN_PRIMARY
            seco = (row[s_col] if s_col < len(row) else "").strip() or UNKNOWN_SECONDARY
            mapping[key] = (prim, seco)
    return mapping

def extract_S_key_from_fastq_name(run_name: str) -> Optional[str]:
    """
    From a FASTQ base name containing '_S###_', return a 4-digit CSV key string: e.g., S101 -> '0101'
    """
    m = S_KEY_RE.search(run_name)
    if not m: return None
    n = int(m.group(1))           # 101
    return f"{n:04d}"             # '0101'

# ========================= MAIN ================================

def main():

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fastq_dir = os.path.join(RUN_DIR, FASTQ_SUBDIR)

    files = [os.path.join(fastq_dir, f)
             for f in os.listdir(fastq_dir)
             if f.endswith(".fastq.gz")]

    if not files:
        print("No FASTQ.GZ files found")
        return

    ref_csv_path = os.path.join(fastq_dir, REF_CSV_BASENAME)

    # Load reference vector & extract canonical 50-bp flanks
    try:
        ref_seq  = load_ref_sequence(ref_csv_path, REF_VECTOR_COL)
        ref_up50, ref_dn50, _ = extract_known_ref_flanks(ref_seq, REF_FLANK_LEN)
        print("Reference loaded; extracted 50-bp flanks.")
    except Exception as e:
        print(f"[WARN] Could not load reference flanks: {e}")
        ref_up50 = ref_dn50 = ""

    # Build key -> (Primary, Secondary) map from CSV
    key_to_names = build_primary_secondary_map(
        ref_csv_path,
        BARCODE_KEY_COL,
        SAMPLE_PRIMARY_COL,
        SAMPLE_SECONDARY_COL
    )
    print(f"Loaded {len(key_to_names)} barcode→(Primary,Secondary) mappings from CSV.")

    # Outputs
    filtered_fastq_path = os.path.join(OUTPUT_DIR, "polyA_14417_filtered.fastq.gz")
    summary_path        = os.path.join(OUTPUT_DIR, "polyA_summary_14417.csv")

    write_header = not os.path.exists(summary_path)
    filtered_fastq = gzip.open(filtered_fastq_path, "wt")

    with open(summary_path, "a", newline="", encoding="utf-8") as sf:
        w = csv.writer(sf)
        if write_header:
            w.writerow([
                "trace_name","run_name","samplePrimary","sampleSecondary","wasRC","trace_length",
                "polyA_length","hasFlanking",
                "isReferenceUp","isReferenceDn",
                "avgQMeanUp","avgQMeanDn",
                "simUp","simDn",
                "alnScoreUp","alnScoreDn",
                "coverUp","coverDn",
                "polyASeq","flankUpSeq","flankDnSeq"
            ])

        total_reads = 0
        kept_reads  = 0

        for fq in sorted(files):
            base     = os.path.basename(fq)
            run_name = infer_run_name(base)

            # Convert '_S###_' → '####' to lookup Column A
            csv_key = extract_S_key_from_fastq_name(run_name)
            if csv_key and (csv_key in key_to_names):
                file_primary, file_secondary = key_to_names[csv_key]
            else:
                # Fallback if key missing: unknowns
                file_primary, file_secondary = UNKNOWN_PRIMARY, UNKNOWN_SECONDARY

            print(f">>> {base}  [key={csv_key}]  [Primary={file_primary} | Secondary={file_secondary}]")

            with gzip.open(fq, "rt") as h:
                for rec in SeqIO.parse(h, "fastq"):
                    total_reads += 1
                    name  = rec.id.split()[0]
                    seq   = str(rec.seq)
                    quals = rec.letter_annotations.get("phred_quality", [])

                    runs = all_poly_runs(seq)
                    if not runs:
                        continue

                    # Normalize minus-only (polyT) to plus frame
                    hasA = any(r["kind"]=="polyA" for r in runs)
                    hasT = any(r["kind"]=="polyT" for r in runs)
                    was_rc = 0
                    if (not hasA) and hasT:
                        seq   = reverse_complement(seq)
                        quals = rc_quals(quals)
                        was_rc = 1
                        runs  = all_poly_runs(seq)
                        if not runs:
                            continue
                        hasA = any(r["kind"]=="polyA" for r in runs)
                        if not hasA:
                            continue

                    a_runs = [r for r in runs if r["kind"]=="polyA"]
                    if not a_runs:
                        continue

                    # Strict A-only polyA
                    longest = max(a_runs, key=lambda r: r["length"])
                    win_s, win_e = longest["start"], longest["end"]
                    # tighten to contiguous A-only segment
                    s, e = win_s, win_e
                    best_len = 0; best_s = s; i = s
                    while i < e:
                        if seq[i] != 'A':
                            i += 1; continue
                        j = i
                        while j < e and seq[j] == 'A':
                            j += 1
                        if (j - i) > best_len:
                            best_len = (j - i); best_s = i
                        i = j
                    tight_s, tight_e = (best_s, best_s + best_len) if best_len > 0 else (win_s, win_s)
                    if tight_e <= tight_s:
                        continue

                    flankUpSeq = seq[:tight_s]
                    flankDnSeq = seq[tight_e:]
                    polyASeq   = seq[tight_s:tight_e]
                    polyA_len  = len(polyASeq)

                    # hasFlanking
                    if flankUpSeq and flankDnSeq: hasFlanking = "both"
                    elif flankUpSeq:              hasFlanking = "up"
                    elif flankDnSeq:              hasFlanking = "down"
                    else:                         hasFlanking = "none"

                    # QMeans (0 if absent)
                    avgQMeanUp = mean_q(quals, 0, tight_s) if flankUpSeq else 0.0
                    avgQMeanDn = mean_q(quals, tight_e, len(seq)) if flankDnSeq else 0.0

                    # Align metrics (zeros if absent)
                    simUp = scoreUp = covUp = 0.0
                    simDn = scoreDn = covDn = 0.0
                    isRefUp = 0; isRefDn = 0

                    if USE_ALIGNER and ref_up50 and flankUpSeq:
                        simUp, scoreUp, covUp = align_flank_to_ref(flankUpSeq, ref_up50)
                        if (simUp >= REF_MIN_IDENTITY) and (covUp >= REF_MIN_COVER_FRAC) and (scoreUp >= REF_MIN_SCORE):
                            isRefUp = 1
                    if USE_ALIGNER and ref_dn50 and flankDnSeq:
                        simDn, scoreDn, covDn = align_flank_to_ref(flankDnSeq, ref_dn50)
                        if (simDn >= REF_MIN_IDENTITY) and (covDn >= REF_MIN_COVER_FRAC) and (scoreDn >= REF_MIN_SCORE):
                            isRefDn = 1

                    # Display sequences (no flank if absent)
                    disp_flankUpSeq = flankUpSeq if flankUpSeq else NA_STR
                    disp_flankDnSeq = flankDnSeq if flankDnSeq else NA_STR

                    # Write summary row
                    w.writerow([
                        name,
                        run_name,
                        file_primary,
                        file_secondary,
                        was_rc,
                        len(seq),
                        polyA_len,
                        hasFlanking,
                        isRefUp, isRefDn,
                        f"{avgQMeanUp:.2f}",
                        f"{avgQMeanDn:.2f}",
                        f"{simUp:.3f}",
                        f"{simDn:.3f}",
                        f"{scoreUp:.2f}",
                        f"{scoreDn:.2f}",
                        f"{covUp:.3f}",
                        f"{covDn:.3f}",
                        polyASeq,
                        disp_flankUpSeq,
                        disp_flankDnSeq
                    ])

                    # Write normalized read to filtered FASTQ
                    rec.seq = Seq(seq)
                    rec.letter_annotations["phred_quality"] = quals
                    SeqIO.write(rec, filtered_fastq, "fastq")
                    kept_reads += 1

    filtered_fastq.close()

    print("\nDONE.")
    print(f"Total reads scanned     : {total_reads}")
    print(f"PolyA reads written     : {kept_reads}  -> {filtered_fastq_path}")
    print(f"Summary CSV             : {summary_path}")
    print(f"Ref thresholds (Up/Dn)  : ID>={REF_MIN_IDENTITY}, COVER>={REF_MIN_COVER_FRAC}, SCORE>={REF_MIN_SCORE}")

if __name__ == "__main__":
    main()
