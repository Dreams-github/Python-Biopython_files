
"""Next task3:
 identify whether the expected REFSEQ region (espl last 70bases) "TGATAATAGGCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGC"
 occurs in the seq reads,  immediately the upstream of polyA, in absent of polyA or anywhere else in the reads?

 idea:
 Read the fastq.gz file, -> define function
 1st align the read against to expected REFSEQ_REGION (119 bases) using local pairwise alignment
 if the read aligning  is not meeting with the alignment criteria, reverse the read seq and its quality and perform local pariwise alignment
 if still not meeting with alignment criteria then return no alignment
 extract all matching read with refseq  -> define function
 classify which part of REFSEQ_REGION is aligning with read( full or last 70bases or 49Frist bases or partial bases(anywher from full))
 perform mean quality score (phread quality) and alignment metrics (identity, score, coverage) -> defined function
 then written output summary data and normalized fastq file.
 CSV file output: "RunName","trace_name","orientation",
                  "read_start","read_end","read_len","read_seq",
                  "ref_start","ref_end","ref_len","ref_seq",
                  "aligned_bases", "classify_align",
                  "similarity","alnScore","read_coverage","ref_coverage",
                  "avgQmean_read", "read_Qstring", "avgQmean_ref",
                  "full_read_seq", "seq_len", "full_read_Qstring"

 Findings/ classifying:
 Full_REFSEQ_Align
 Last70_REFSEQ_Align
 Last49_REFSEQ_Align
 Partial_REFSEQ_Align
------------------------------------------------------------------------------------------------------------------------------------------------
 #8/20/2025: After discussion with Jeff
 1. Identify and extract the bases or string of read seqs that matches from 70 to 120 bases of the reference sequence
    - min 40bp
    - give +5 bases of tolerance for the end seq
    - & perform its quality score
    note: remove 0 based index and perform everything with 1 based index for ref and read seqs
    - modify the script accordingly to implement logics
 2. not now but later: identify and extract the regions of  1 to 50bp, 30 to 40bp, 40 to 80bp and 70 to 120bp of the reference sequence from the reads
    - give +5 bases of tolerance for the end seq
    - & perform its quality scores and also the quality score of the overlap regions -> (keeping this task 2 for later)
----------------------------------------------------------------------------------------------------------------------------------------------------
#8/28/2025:
For what you have underlined, we care that the alignment is a sizable alignment >=40 or 50 for a respective ref_len value.
Thus, the check should shift from ref_start to refLength and check accordingly.
To aligning near the 3' terminal end. For the 2nd half ref_end >= ref_len - coord_tolerance should either be just '>' or include a '+1' to the check.
note will try this on monday after confirming with Jeff on 0-based or 1-based indexing
  """

# importing all required libraries

import sys
import gzip, csv

#import pandas as pd
from pathlib import Path
from statistics import mean
from typing import Optional, List, Tuple, Dict

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align
from matplotlib import use
from numpy.ma.core import identity

# -----------------USER SETTINGS-------------

# defining input and output paths
ROOT_DIR = Path(r"G:\Platform Innovation - Applied Biology\Staff\swapna\Moderna_polyA_NGS_Assessment\MD")
OUTPUT_DIR = ROOT_DIR / "identified_refseq_regions"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# Providing Refseq region to identify from the read seq traces
REFSEQ_REGION = ("TGATAATAGGCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGC")  # len(REFSEQ_REGION) is 119

# locating last 70 adn first 49 bases of the Refseq region extraction from seq reads
#REFSEQ_LAST70BASES = REFSEQ_REGION[-70:]
#REFSEQ_FIRST49BASES = REFSEQ_REGION[:49]
REFSEQ_LAST_BASES = REFSEQ_REGION[69:] # using this variable for al.align as a refseq to extract the aligned ref bases (70 to all the 3'end of refseq) from reads - 8/27/26
                                    # 'TTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGC'
# Aligner thresholds - check for other possibilities to see how to implement pairwise2 to find the best results for matching score
USE_ALIGNER = True
#WILDCARD = None  # try how to use WILDCARD for mismatch score
REF_MIN_IDENTITY = 0.80
REF_MIN_COVER_FRAC = 0.70
REF_MIN_SCORE = 12.0
REF_MIN_ALIGNED_BASES = 40 #changed from 20 to 30 to 40 8/26th - 50 8/28th changing to 40 9/08
ALIGN_MATCH_SCORE = 1.0
ALIGN_MISMATCH_SCORE = -0.5
ALIGN_OPEN_GAP_SCORE = -1.0
ALIGN_EXTEND_GAP_SCORE = -0.5

RUN_SANITY_CHECK = True #keep true for sanity_check, turn false for full run
TEST_READ_LIMIT = 10 #keep 10 for testing, set to None for full run

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

# --------- reverse complement & reverse Q_score------
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

# ---------Aligning with Ref_seq___________
def init_aligner():
    al = Align.PairwiseAligner()  # creates PA object and stores the alignment parameters to be used for the PA
    al.mode = "local"
    al.match_score = ALIGN_MATCH_SCORE
    al.mismatch_score = ALIGN_MISMATCH_SCORE
    al.open_gap_score = ALIGN_OPEN_GAP_SCORE
    al.extend_gap_score = ALIGN_EXTEND_GAP_SCORE
    return al

#---Modifying the existing function ---

#---normalizing the blocks before defining alignment---
def normalize_blocks(blocks):
    """
       Convert PairwiseAligner aligned blocks into normal Python list of int tuples.
       """
    if hasattr(blocks, "tolist"):
        return [(int(a), int(b)) for (a, b) in blocks.tolist()]
    return [(int(a), int(b)) for (a, b) in (blocks or [])]

#---aligning full read seq to the expected REFSEQ_REGION to locate the candidate_refseq
def locate_read_by_alignment(seq: str) -> Optional[Tuple[int, int, str, int, int, str,int, float, float, float, float, float, int]]:
    """
       Align the full read sequence to the expected REFSEQ_REGION.

       Returns:
           read_start,
           read_end,
           read,
           ref_start,
           ref_end,
           ref,
           matches,
           identity,
           identity_maxlen,
           score,
           read_coverage,
           ref_coverage,
           aligned_bases

       read_start/read_end:
           extracted location inside sequencing read

       ref_start/ref_end:
           aligned location inside expected REFSEQ_REGION coordinate window 70 to 120(+5 bases tolerance)
       """

    if not seq:
        return None

    al = init_aligner()
    # target = expected REFSEQ_REGION
    # query = full sequencing read
    alns = al.align(REFSEQ_LAST_BASES, seq) #replaced all REFSEQ_REGION (119) to  REFSEQ_LAST_BASES = REFSEQ_REGION[69:] as a target aligner
    if alns is None or len(alns) == 0:
        return None
    aln = alns[0]  # takes top best alignment

    try:
        t_blocks, q_blocks = aln.aligned  # may be numpy arrays
    except Exception:
        return None

    t_blocks = normalize_blocks(t_blocks)
    q_blocks = normalize_blocks(q_blocks)

    if len(t_blocks) == 0 or len(q_blocks) == 0:
        return None

    ref_start = t_blocks[0][0]
    ref_end = t_blocks[-1][1]
    ref = REFSEQ_LAST_BASES[ref_start:ref_end]

    read_start = q_blocks[0][0]
    read_end = q_blocks[-1][1]
    read = seq[read_start:read_end]

    matches = 0;
    aligned_bases = 0
    for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):  # zip loops over two lists together
        ref_part = REFSEQ_LAST_BASES[ts:te]
        read_part = seq[qs:qe]

        matches += sum(1 for a, b in zip(ref_part, read_part) if a == b)  # compare base by base
        aligned_bases += (qe - qs)
        #OR
        #if len(ref_part) > len(read_part):
            #aligned_bases += (te - ts)
        #else:
            #aligned_bases += (qe - qs) both ways are giving same results

    read_len = len(read)
    ref_len = len(ref)

    #if ref_len > read_len:
        #identity_maxlen = ref_len
    #else:
        #identity_maxlen = read_len

    identity = 0.0 if aligned_bases == 0 else (matches / aligned_bases)
    #identity_maxlen = 0.0 if identity_maxlen == 0 else (matches / identity_maxlen) #disabling 9/18
    read_coverage = 0.0 if len(seq) == 0 else (aligned_bases / len(seq)) #9/8 - modified as len(seq) to calculate read_coverage based on full trace length
    ref_coverage = 0.0 if len(REFSEQ_REGION) == 0 else (aligned_bases / len(REFSEQ_REGION)) #9/8 - modified as len(REFSEQ_REGION) to calculate ref_coverage based on full ref length
    score = float(aln.score)

    #Filters to avoid weak alignments
    if aligned_bases < REF_MIN_ALIGNED_BASES:
        return None
    if identity < REF_MIN_IDENTITY:
        return None
    #if read_coverage < REF_MIN_COVER_FRAC: 9/8 commented this out
        #return None
    if score < REF_MIN_SCORE:
        return None
    return (
        read_start,
        read_end,
        read,
        ref_start,
        ref_end,
        ref,
        matches,
        identity,
        identity_maxlen,
        score,
        read_coverage,
        ref_coverage,
        aligned_bases
    )

#classify the alignment into two ways: last region of the refseq (70 to 119bp) or remaining part of the refseq (1 to 69bp) --
#def classify_align(ref_start: int, ref_end: int) -> str:
    #if ref_start == -1 or ref_end == -1:
        #return "No_REFSEQ_Align"
    #refseq_len = len(REFSEQ_REGION)
    #coord_tolerance = 5
    #if ref_start >= 69 - coord_tolerance and ref_end >= refseq_len - coord_tolerance:
        #return "REFSEQ_Align_from_70-119"
    #else:
        #return "Align_Anywhere_REFSEQ_Region"

#---Sanity check----
def sanity_check():

    print("\n===== SANITY CHECK =====")

    #Test_classify_align
    test_seq_full = REFSEQ_LAST_BASES + ("A" * 30)


    test_cases = [
        ("anywhere_refseqalign_plus_polyA", test_seq_full)
   ]

    for name, test_seq in test_cases:
        result = locate_read_by_alignment(test_seq)

        print("\nTest:", name)

        if result:
            read_start, read_end, read, ref_start, ref_end, ref, matches, identity, identity_maxlen, score, read_coverage, ref_coverage, aligned_bases = result

            print("read_start:", read_start)
            print("read_end:", read_end)
            print("read:", read)
            print("ref_start:", ref_start)
            print("ref_end:", ref_end)
            print("ref:", ref)
            print("matches:", matches)
            print("identity:", identity)
            print("identity_maxlen:", identity_maxlen)
            print("score:", score)
            print("read_coverage:", read_coverage)
            print("ref_coverage:", ref_coverage)
            print("aligned_bases:", aligned_bases)

        else:
            print("No extracted REFSEQ detected")

    print("\n===== END SANITY CHECK =====\n")

# ----------------------------
# ---Defining main function---
# ----------------------------
def main():
    if RUN_SANITY_CHECK:
        sanity_check()

    fastq_paths = list(ROOT_DIR.rglob("*_001.fastq.gz"))  # rglobe()=recursive search-goes through all subfolders
    if not fastq_paths:
        print(f"[WARN] no fastq.gz file found under: {ROOT_DIR}")  # print the warning message
        sys.exit(0)  # stops the script safely, and '0' means exit normally, no fastq file found is not treat as an error crash

    # output files
    allfiles_filtered_fastq_path = OUTPUT_DIR/"Allfiles_Normalized_WIP.fastq"
    allfiles_summary_path = OUTPUT_DIR/"Allfiles_summary_WIP.csv"

    total_reads = 0  # counts ALL reads processed
    total_refseq_reads = 0  # counts only valid REFSEQ + polyA and REFSEQ + no_polyA reads written to output


    #write output normalized fastq file and header for output summary csv file
    with open(allfiles_summary_path, "w", newline="")as csv_handle, open(allfiles_filtered_fastq_path, "w") as norm_fastq_handle:
        writer = csv.writer(csv_handle)
        writer.writerow([
            "RunName","file_name","trace_name","orientation",
            "read_start","read_end","read_len","read_seq",
            "ref_start","ref_end","ref_len","ref_seq",
            "matches", "aligned_bases",
            "identity", "identity_maxlen","alnScore","read_coverage","ref_coverage",
            "avgQmean_read", "read_Qstring",
            "full_read_seq", "seq_len", "full_read_Qstring"
        ])

        # read the fastq file scan_fastq_for_reads
        for fastq_path in fastq_paths:
            run_name = get_run_id(fastq_path)
            fastq_file = fastq_path.name
            print("\nProcessing file:", fastq_path)

            with gzip.open(fastq_path, "rt") as handle:
                for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                    if TEST_READ_LIMIT is not None and i >= TEST_READ_LIMIT:
                        break
                    total_reads += 1
                    trace_name = record.id.split()[0]
                    seq = str(record.seq)
                    quals = record.letter_annotations["phred_quality"]
                    quals_string = record.format("fastq").splitlines()[3] #***check this
                    seq_len = len(seq)

                    #Read Alignment against REFSEQ extraction
                    read_result = locate_read_by_alignment(seq) # this will align with original read i.e, fwd orientation
                    if read_result is None: # this means alignment fail then below steps perform rc and rq
                        rc_seq = reverse_complement(seq)
                        rc_q = rc_quals(quals)
                        rc_q_string = quals_string[::-1]

                        read_result = locate_read_by_alignment(rc_seq) #nw aligne the normalized read against ref

                        if read_result is not None:
                            seq = rc_seq
                            quals = rc_q
                            quals_string = rc_q_string
                            orientation = "RC"
                        else:
                            orientation = "NO_ALIGNMENT"

                    else:
                        orientation = "FWD"

                    if read_result:
                        (
                            read_start,
                            read_end,
                            read,
                            ref_start,
                            ref_end,
                            ref,
                            matches,
                            identity,
                            identity_maxlen,
                            score,
                            read_coverage,
                            ref_coverage,
                            aligned_bases
                        ) = read_result #this is unpacking tuple

                        read_start_1BI = read_start + 1 # changing to 1 based index
                        read_end_1BI = read_end
                        read_len = len(read)
                        ref_start_1BI = ref_start + 70 # changing to 1 based index
                        ref_end_1BI = ref_end + 69
                        ref_len = len(ref)
                        isRef = 1 #aligned with refseq
                        total_refseq_reads += 1 #count that read
                    else:
                        read_start_1BI = read_end_1BI = 0
                        read = ""
                        read_len = 0

                        ref_start_1BI = ref_end_1BI = -1 #no alignment match
                        ref = "NONE"
                        ref_len = 0
                        identity = score = read_coverage = ref_coverage = 0.0
                        identity_maxlen = matches = 0.0
                        aligned_bases = 0
                        isRef = 0


                    # Quality metrics for read only not ref
                    avgQMean_read= mean_q(quals,read_start, read_end) if read else 0.0
                    avgQMean_read_string = quals_string[read_start:read_end] if read else ""
                    #avgQMean_ref =  mean_q(quals, ref_start, ref_end) if ref != "NONE" else 0.0
                    #avgQMean_ref_string = quals_sting(ref_start, ref_end) if ref != "NONE" else ""

                    #classifying alignment and read
                    #category_align = classify_align(ref_start, ref_end)


                    # Temporary debug prints for first few reads
                    if TEST_READ_LIMIT is not None:
                        print("\n----- REAL READ TEST -----")
                        print("Read:", record.id)
                        print("trace_name:", trace_name)
                        print("orientation:", orientation)
                        print("read_start:", read_start)
                        print("read_end:", read_end)
                        print("read_len:", read_len)
                        print("read_seq:", read)
                        print("ref_start:", ref_start)
                        print("ref_end:", ref_end)
                        print("ref_len:", ref_len)
                        print("ref:", ref)
                        print("aligned_bases:", aligned_bases)
                        #print("category_align:", category_align)
                        print("sim:", identity)
                        print("score:", score)
                        print("read_coverage:", read_coverage)
                        print("ref_coverage:", ref_coverage)
                        print("seq:", seq)
                        print("seq_len:", len(seq))


                    # Write summary row
                    if read_len > 0:
                        writer.writerow([
                            run_name,
                            fastq_file,
                            trace_name,
                            orientation,
                            read_start_1BI,
                            read_end_1BI,
                            read_len,
                            read,
                            ref_start_1BI,
                            ref_end_1BI,
                            ref_len,
                            ref,
                            matches,
                            aligned_bases,
                            #category_align,
                            f"{identity:.3f}",
                            f"{identity_maxlen:.3f}",
                            f"{score:.2f}",
                            f"{read_coverage:.2f}",
                            f"{ref_coverage:.2f}",
                            f"{avgQMean_read:.2f}",
                            avgQMean_read_string,
                            #f"{avgQMean_ref:.2f}",
                            seq,
                            seq_len,
                            quals_string
                        ])

                    # Write normalized read to filtered FASTQ
                    record.seq = Seq(seq)
                    record.letter_annotations = {}
                    record.letter_annotations["phred_quality"] = quals
                    SeqIO.write(record, norm_fastq_handle, "fastq")

        print("\nDONE.")
        print("Total reads processed:", total_reads)
        print("Total REFSEQ-like reads:", total_refseq_reads)
        print("Allfiles_summary_WIP.csv:", allfiles_summary_path)
        print("Allfiles_Normalized_WIP.fastq:", allfiles_filtered_fastq_path)
if __name__ == "__main__":
    main()
