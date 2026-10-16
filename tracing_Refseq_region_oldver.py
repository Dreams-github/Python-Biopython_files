# om sai ram

"""Next task3:
 identify whether the expected REFSEQ region "TGATAATAGGCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGC"
 occurs immediately upstream of polyA, and if not, what other flanking sequence is present instead, also locate where else this regions present without polyA?

 idea:
 Read the fastq.gz file, identifying the longest polyA and polyT, perform RC and RQ for the polyT and normalize everything to polyA forward direction,
 Extract the upstream flank region, align with that part of ref vecseq region,  perform quality score and written output summary data and normalized fastq file.

 Findings:
 REFSEQ-polyA (length of polyA)
 REFSEQ-NopolyA(some random seq)
 No REFSEQ(some different flank should be the same bp len of refseq)-polyA """

# importing all required libraries

import sys
import gzip, csv

import pandas as pd
from pathlib import Path
from statistics import mean
from typing import Optional, List, Tuple, Dict

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align

# -----------------USER SETTINGS-------------

# defining input and output paths
ROOT_DIR = Path(r"G:\Platform Innovation - Applied Biology\Staff\swapna\Moderna_polyA_NGS_Assessment\MD")
OUTPUT_DIR = ROOT_DIR / "identified_flank_regions"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Providing Refseq region to identify from the seq traces
REFSEQ_REGION = "TGATAATAGGCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGC"  # len(REFSEQ_REGION) is 119 noting for my reference

# Column indices (0-based)
# WELLPOSITION_COL = 0
# PRIMARY_COL = 3
# SECONDARY_COL = 4
# SAMPLE_COL = 5
# REF_COL = 8


# Refseq flank extraction around polyA or without polyA
# REFSEQ_FLANK_MINLEN = 25
# REFSEQ_FLANK_MAXLEN = 119

# flank extraction around polyA
FLANK_LEN = len(REFSEQ_REGION)

# Aligner thresholds - check for other possibilities to see how to implement pairwise2 to find the best results for matching score
USE_ALIGNER = True
#WILDCARD = None  # trying to use WILDCARD for mismatch score
REF_MIN_IDENTITY = 0.90
REF_MIN_COVER_FRAC = 0.70
REF_MIN_SCORE = 12.0
ALIGN_MATCH_SCORE = 1.0
ALIGN_MISMATCH_SCORE = -0.5
ALIGN_OPEN_GAP_SCORE = -1.0
ALIGN_EXTEND_GAP_SCORE = -0.5


# ---------------------------------------------------
# -------------------DEFINING FUNCTIONS--------------
# ---------------------------------------------------

# -------------loading fastq path files---------------
def get_run_id(fastq_path: Path) -> str:
    """
    Extract run folder name, e.g.
    MD/250821_M04095_0651_000000000-GT52R/Fastq files/*.fastq.gz
    """
    for parent in fastq_path.parents:  # climbs up directory tree - .parents all parent directories
        if parent.parent == ROOT_DIR:  # finds run-level folder
            return parent.name  # extract run ID
    return ""

# --------Longest polyA/T rusn coordinates from seq reads-----------
def longest_polyA_run(s: str) -> Optional[Tuple[int, int]]:
    best_len = cur = 0;
    best_s = best_e = None;
    start = None
    for i, ch in enumerate(s):
        if ch == 'A' :
            if cur == 0: start = i
            cur += 1
            if cur > best_len: best_len = cur; best_s = start; best_e = i + 1
        else:
            cur = 0;
            start = None
    return (best_s, best_e) if best_s is not None else None

def longest_polyT_run(s: str) -> Optional[Tuple[int, int]]:
    best_len = cur = 0;
    best_s = best_e = None;
    start = None
    for i, ch in enumerate(s):
        if ch == 'T':
            if cur == 0: start = i
            cur += 1
            if cur > best_len: best_len = cur; best_s = start; best_e = i + 1
        else:
            cur = 0;
            start = None
    return (best_s, best_e) if best_s is not None else None


# ---------if polyT identified performs reverse complement & reverse Q_score------
def reverse_complement(s: str) -> str:
    return str(Seq(s).reverse_complement())


def rc_quals(q: List[int]) -> List[int]:
    return list(reversed(q))


# --------avg quality score------------

# clamp_span ensures that a given (start, end) range stays within valid bounds of a sequence
def clamp_span(n: int, start: int, end: int) -> Tuple[int, int]:
    s = max(0, start);
    e = min(n, end)
    return (s, e) if e >= s else (s, s)


def mean_q(quals: List[int], start: int, end: int) -> Optional[float]:
    s, e = clamp_span(len(quals), start, end)
    if e <= s: return None
    return float(mean(quals[s:e]))

#-------Extract upstream flanking region of polyA after normalized---------

def extract_upstream_flank(seq: str, polyA_start: int,) -> Tuple[int, int, str]:
    start = max(0, polyA_start - FLANK_LEN)
    end = polyA_start
    return start, end, seq[start:end]

# ---------Aligning with Ref_seq___________
def init_aligner():
    al = Align.PairwiseAligner()  # creates PA object and stores the alignment parameters to be used for the PA
    al.mode = "local"
    al.match_score = ALIGN_MATCH_SCORE
    al.mismatch_score = ALIGN_MISMATCH_SCORE
    al.open_gap_score = ALIGN_OPEN_GAP_SCORE
    al.extend_gap_score = ALIGN_EXTEND_GAP_SCORE
    return al


def align_flank_to_ref(flank: str, ref_seq: str) -> Tuple[float, float, float]:
    if not flank: return (0.0, 0.0, 0.0)
    al = init_aligner()
    alns = al.align(ref_seq, flank)
    if alns is None or len(alns) == 0: return (0.0, 0.0, 0.0)
    aln = alns[0]  # takes top best alignment
    try:
        t_blocks, q_blocks = aln.aligned  # may be numpy arrays
    except Exception:
        return (0.0, float(aln.score), 0.0)
    if hasattr(t_blocks, 'tolist'):
        t_blocks = [(int(a), int(b)) for (a, b) in t_blocks.tolist()]
    else:
        t_blocks = [(int(a), int(b)) for (a, b) in (t_blocks or [])]
    if hasattr(q_blocks, 'tolist'):
        q_blocks = [(int(a), int(b)) for (a, b) in q_blocks.tolist()]
    else:
        q_blocks = [(int(a), int(b)) for (a, b) in (q_blocks or [])]
    if len(t_blocks) == 0 or len(q_blocks) == 0: return (0.0, float(aln.score), 0.0)
    matches = 0;
    aligned_q_len = 0
    for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):  # zip loops over two lists together
        seg_t = ref_seq[ts:te];
        seg_q = flank[qs:qe]
        matches += sum(1 for a, b in zip(seg_t, seg_q) if a == b)  # compare base by base
        aligned_q_len += (qe - qs)
    identity = 0.0 if aligned_q_len == 0 else (matches / aligned_q_len)
    cover = 0.0 if len(flank) == 0 else (aligned_q_len / len(flank))
    return (identity, float(aln.score), cover)  # aln.score is used to calculate the alignment score between two seqs.

#---locating ref vecseq regions even in absence of polyA----
def locate_refseq_region(seq: str) -> Optional[Tuple[int, int]]:
    region = seq.find(REFSEQ_REGION)
    if region >= 0:
        return (region, region + len(REFSEQ_REGION))
    return None

#-----classifying the reads-----
def classify_read(polyA_found:bool,ref_match: bool)-> str:
    if ref_match and polyA_found:
        return "REFSEQ_polyA"
    elif ref_match and not polyA_found:
        return "REFSEQ_noPolyA"
    elif polyA_found and not ref_match:
        return "No_REFSEQ_polyA"
    else:
        return "Neither"

# ----------------------------
# ---Defining main function---
# ----------------------------
def main():
    fastq_paths = list(ROOT_DIR.rglob("*.fastq.gz"))  # rglobe()=recursive search-goes through all subfolders
    if not fastq_paths:
        print(f"[WARN] no fastq.gz file found under {ROOT_DIR}")  # print the warning message
        sys.exit(
            0)  # stops the script safely, and '0' means exit normally, no fastq file found is not treat as an error crash

    # output files
    allfiles_filtered_fastq_path = OUTPUT_DIR/"Normalized_fastq.gz"
    allfiles_summary_path = OUTPUT_DIR/"Allfiles_summary.csv"

    #write header for output summary csv file
    with open(allfiles_summary_path, "w", newline="", encoding="utf-8") as outf:
        writer = csv.writer(outf)
        writer.writerow([
            "RunName", "trace name", "trace length", "Performed RC & RQ", "polyA length", "polyA seq",
            "extracted_flankUpseq", "extracted_flankupseq len", "align_Refseq (0/1)",
            "avgQmeanUp", "alnScore", "Refseq regions in absence of polyA"
        ]) #missing some header titles which are not matching with the output data

        total_reads = 0  # counts ALL reads processed
        kept_reads = 0  # counts only valid polyA reads written to output

        # read the fastq file scan_fastq_for_reads
        for fastq_path in fastq_paths:
            run_name = get_run_id(fastq_path)

            with gzip.open(fastq_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    total_reads += 1
                    trace_name = record.id.split()[0]
                    seq = str(record.seq)
                    quals = record.letter_annotations["phred_quality"]

                    #finding longest polyA/T
                    polyA = longest_polyA_run(seq)
                    polyT = longest_polyT_run(seq)

                    was_rc = 0 #Normalizing into forward orientation

                    if (not polyA) and polyT:
                        seq = reverse_complement(seq)
                        quals = rc_quals(quals)

                        was_rc = 1

                        polyA = longest_polyA_run(seq) #recalculating polyA after RC

                    refseq_no_polyA = ""
                    if polyA is None:
                        ref_coords = locate_refseq_region(seq)

                        if ref_coords:
                            ref_start, ref_end = ref_coords
                            refseq_no_polyA = seq[ref_start:ref_end]
                        continue
                            # identifying poly A coordinates
                    polyA_start, polyA_end = polyA
                    polyA_len = polyA_end - polyA_start
                    polyA_seq = seq[polyA_start:polyA_end]

                    #identifying the upstream flanking region and its coordinates
                    flank_start, flank_end, flank_seq = extract_upstream_flank(seq, polyA_start)
                    flank_len = len(flank_seq)

                    # QMeans for flank and polyA
                    avgQMeanup = mean_q(quals,flank_start, flank_end) if flank_seq else 0.0
                    #avgQMeanpolyA =  mean_q(quals, polyA_start, polyA_end) if polyA else 0.0

                    # Align metrics (zeros if absent) labeling
                    simUp = scoreUp = covUp = 0.0
                    isRefUp = 0;

                    if USE_ALIGNER and REFSEQ_REGION and flank_seq:
                        simUp, scoreUp, covUp = align_flank_to_ref(flank_seq, REFSEQ_REGION)
                        if (simUp >= REF_MIN_IDENTITY) and (covUp >= REF_MIN_COVER_FRAC) and (scoreUp >= REF_MIN_SCORE):
                            isRefUp = 1

                    category = classify_read( polyA_found=True,
                        ref_match=(isRefUp == 1)
                    )

                    # Write summary row
                    writer.writerow([
                        run_name,
                        trace_name,
                        len(seq),
                        was_rc,
                        polyA_len,
                        polyA_seq,
                        flank_seq,
                        flank_len,
                        isRefUp,
                        f"{avgQMeanup:.2f}",
                        f"{simUp:.3f}",
                        f"{scoreUp:.2f}",
                        f"{covUp:.3f}",
                        refseq_no_polyA,
                        category
                    ])

                    # Write normalized read to filtered FASTQ
                    record.seq = Seq(seq)
                    record.letter_annotations["phred_quality"] = quals
                    SeqIO.write(record, allfiles_filtered_fastq_path, "fastq")
                    kept_reads += 1

        allfiles_filtered_fastq_path.close()

        print("\nDONE.")
        print(f"Total reads scanned     : {total_reads}")
        print(f"PolyA reads written     : {kept_reads}  -> {allfiles_filtered_fastq_path}")
        print(f"Summary CSV             : {allfiles_summary_path}")
        print(f"Ref thresholds (Up/Dn)  : ID>={REF_MIN_IDENTITY}, COVER>={REF_MIN_COVER_FRAC}, SCORE>={REF_MIN_SCORE}")

if __name__ == "__main__":
    main()
