## om sai ram

"""
Task 3:
Identify whether the expected REFSEQ region occurs immediately upstream of polyA, and
is present in the seq read without polyA, or is absent.

New logic:
1. Read fastq.gz files.
2. Detect longest polyA and polyT.
3. If polyT is stronger than polyA, reverse-complement read and reverse Q scores.
4. Align the full read against expected REFSEQ_REGION.
5. Extract the read segment that aligns to REFSEQ_REGION coordinates 0 to 119.
6. Classify alignment as Full, Last70, First49, Partial, or No_REFSEQ.
7. Extract downstream sequence after the extracted REFSEQ-like region.
8. Write output CSV and normalized FASTQ.
"""

import sys
import gzip
import csv
from pathlib import Path
from statistics import mean
from typing import Optional, List, Tuple

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Align


# ----------------- USER SETTINGS -----------------

ROOT_DIR = Path(r"G:\Platform Innovation - Applied Biology\Staff\swapna\Moderna_polyA_NGS_Assessment\MD")
OUTPUT_DIR = ROOT_DIR / "identified_refseq_regions"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

REFSEQ_REGION = (
    "TGATAATAGGCTGGAGCCTCGGTGGCCATGCTTCTTGCCCCTTGGGCCTCCCCCCAGCCCCTCCTCCCCTTCCTGCACCCGTACCCCCGTGGTCTTTGAATAAAGTCTGAGTGGGCGGC"
)

REFSEQ_LAST70BASES = REFSEQ_REGION[-70:]
REFSEQ_FIRST49BASES = REFSEQ_REGION[:49]

MIN_PolyT_FOR_RC = 10
ext_seq_len = 100

USE_ALIGNER = True

REF_MIN_IDENTITY = 0.80
REF_MIN_COVER_FRAC = 0.70
REF_MIN_SCORE = 12.0
REF_MIN_ALIGNED_BASES = 20

ALIGN_MATCH_SCORE = 1.0
ALIGN_MISMATCH_SCORE = -0.5
ALIGN_OPEN_GAP_SCORE = -1.0
ALIGN_EXTEND_GAP_SCORE = -0.5

RUN_SANITY_CHECK = False
TEST_READ_LIMIT = None          # keep 10 for testing, set to None for full run


# ----------------- FUNCTIONS -----------------

def get_run_id(fastq_path: Path) -> str:
    """
    Extract run folder name.
    Example:
    MD/250821_M04095_0651_000000000-GT52R/Fastq files/*.fastq.gz
    """
    for parent in fastq_path.parents:
        if parent.parent == ROOT_DIR:
            return parent.name
    return ""


def longest_polyA_run(s: str) -> Optional[Tuple[int, int]]:
    best_len = cur = 0
    best_s = best_e = None
    start = None

    for i, ch in enumerate(s):
        if ch == "A":
            if cur == 0:
                start = i
            cur += 1

            if cur > best_len:
                best_len = cur
                best_s = start
                best_e = i + 1
        else:
            cur = 0
            start = None

    return (best_s, best_e) if best_s is not None else None


def longest_polyT_run(s: str) -> Optional[Tuple[int, int]]:
    best_len = cur = 0
    best_s = best_e = None
    start = None

    for i, ch in enumerate(s):
        if ch == "T":
            if cur == 0:
                start = i
            cur += 1

            if cur > best_len:
                best_len = cur
                best_s = start
                best_e = i + 1
        else:
            cur = 0
            start = None

    return (best_s, best_e) if best_s is not None else None


def run_len(coord: Optional[Tuple[int, int]]) -> int:
    if coord is None:
        return 0
    return coord[1] - coord[0]


def reverse_complement(s: str) -> str:
    return str(Seq(s).reverse_complement())


def rc_quals(q: List[int]) -> Listreturn: list(reversed(q))


def clamp_span(n: int, start: int, end: int) -> Tuple[int, int]:
    s = max(0, start)
    e = min(n, end)
    return (s, e) if e >= s else (s, s)


def mean_q(quals: List[int], start: int, end: int) -> float:
    s, e = clamp_span(len(quals), start, end)

    if e <= s:
        return 0.0

    return float(mean(quals[s:e]))


def init_aligner():
    al = Align.PairwiseAligner()
    al.mode = "local"
    al.match_score = ALIGN_MATCH_SCORE
    al.mismatch_score = ALIGN_MISMATCH_SCORE
    al.open_gap_score = ALIGN_OPEN_GAP_SCORE
    al.extend_gap_score = ALIGN_EXTEND_GAP_SCORE
    return al


def normalize_blocks(blocks):
    """
    Convert PairwiseAligner aligned blocks into normal Python list of int tuples.
    """
    if hasattr(blocks, "tolist"):
        return [(int(a), int(b)) for (a, b) in blocks.tolist()]

    return [(int(a), int(b)) for (a, b) in (blocks or [])]


def locate_extracted_refseq_by_alignment(
    seq: str
) -> Optional[Tuple[int, int, str, int, int, float, float, float, float, int]]:
    """
    Align the full read sequence to the expected REFSEQ_REGION.

    Returns:
        read_start,
        read_end,
        extracted_refseq,
        ref_start,
        ref_end,
        identity,
        score,
        extracted_refseq_coverage,
        ref_coverage,
        aligned_bases

    read_start/read_end:
        extracted location inside sequencing read

    ref_start/ref_end:
        aligned location inside expected REFSEQ_REGION coordinate window 0 to 119
    """

    if not seq:
        return None

    al = init_aligner()

    # target = expected REFSEQ_REGION
    # query = full sequencing read
    alns = al.align(REFSEQ_REGION, seq)

    if alns is None or len(alns) == 0:
        return None

    aln = alns[0]

    try:
        t_blocks, q_blocks = aln.aligned
    except Exception:
        return None

    t_blocks = normalize_blocks(t_blocks)
    q_blocks = normalize_blocks(q_blocks)

    if len(t_blocks) == 0 or len(q_blocks) == 0:
        return None

    ref_start = t_blocks[0][0]
    ref_end = t_blocks[-1][1]

    read_start = q_blocks[0][0]
    read_end = q_blocks[-1][1]

    extracted_refseq = seq[read_start:read_end]

    matches = 0
    aligned_bases = 0

    for (ts, te), (qs, qe) in zip(t_blocks, q_blocks):
        ref_part = REFSEQ_REGION[ts:te]
        read_part = seq[qs:qe]

        matches += sum(1 for a, b in zip(ref_part, read_part) if a == b)
        aligned_bases += qe - qs

    identity = 0.0 if aligned_bases == 0 else matches / aligned_bases

    extracted_refseq_coverage = (
        0.0 if len(extracted_refseq) == 0 else aligned_bases / len(extracted_refseq)
    )

    ref_coverage = (
        0.0 if len(REFSEQ_REGION) == 0 else aligned_bases / len(REFSEQ_REGION)
    )

    score = float(aln.score)

    # Filters to avoid extracting random weak alignments
    if aligned_bases < REF_MIN_ALIGNED_BASES:
        return None

    if identity < REF_MIN_IDENTITY:
        return None

    if extracted_refseq_coverage < REF_MIN_COVER_FRAC:
        return None

    if score < REF_MIN_SCORE:
        return None

    return (
        read_start,
        read_end,
        extracted_refseq,
        ref_start,
        ref_end,
        identity,
        score,
        extracted_refseq_coverage,
        ref_coverage,
        aligned_bases
    )


def classify_align(ref_start: int, ref_end: int) -> str:
    """
    Classify where the extracted extracted aligns within REFSEQ_REGION coordinates 0 to 119.
    """

    if ref_start == -1 or ref_end == -1:
        return "No_REFSEQ_Align"

    ref_len = len(REFSEQ_REGION)
    coord_tolerance = 2

    if ref_start <= coord_tolerance and ref_end >= ref_len - coord_tolerance:
        return "Full_REFSEQ_Align"

    elif ref_start >= 49 - coord_tolerance and ref_end >= ref_len - coord_tolerance:
        return "Last70_REFSEQ_Align"

    elif ref_start <= coord_tolerance and ref_end <= 49 + coord_tolerance:
        return "First49_REFSEQ_Align"

    else:
        return "Partial_REFSEQ_Align"


def classify_read(polyA_found: bool, ref_match: bool) -> str:
    if ref_match and polyA_found:
        return "REFSEQ_polyA"

    elif ref_match and not polyA_found:
        return "REFSEQ_noPolyA"

    else:
        return "No_REFSEQ"


def sanity_check():
    print("\n===== SANITY CHECK =====")

    test_seq_full = REFSEQ_REGION + ("A" * 30)
    test_seq_last70 = REFSEQ_LAST70BASES + ("A" * 30)
    test_seq_first49 = REFSEQ_FIRST49BASES + ("C" * 30)

    test_cases = [
        ("full_refseq_plus_polyA", test_seq_full),
        ("last70_plus_polyA", test_seq_last70),
        ("first49_plus_C", test_seq_first49)
    ]

    for name, test_seq in test_cases:
        result = locate_extracted_refseq_by_alignment(test_seq)

        print("\nTest:", name)

        if result:
            read_s, read_e, extracted_refseq, ref_start, ref_end, sim, score, cov, ref_cov, aligned_bases = result
            print("extracted_refseq_len:", len(extracted_refseq))
            print("read_start:", read_s)
            print("read_end:", read_e)
            print("ref_start:", ref_start)
            print("ref_end:", ref_end)
            print("identity:", sim)
            print("score:", score)
            print("extracted_refseq_coverage:", cov)
            print("ref_coverage:", ref_cov)
            print("aligned_bases:", aligned_bases)
            print("category:", classify_align(ref_start, ref_end))
        else:
            print("No extracted REFSEQ detected")

    print("\n===== END SANITY CHECK =====\n")


# ----------------- MAIN -----------------

def main():

    if RUN_SANITY_CHECK:
        sanity_check()

    fastq_paths = list(ROOT_DIR.rglob("*.fastq.gz"))

    if not fastq_paths:
        print("[WARN] No fastq.gz files found under:", ROOT_DIR)
        sys.exit(0)

    summary_csv_path = OUTPUT_DIR / "refseq_alignment_summary.csv"
    normalized_fastq_path = OUTPUT_DIR / "allfiles_normalized.fastq"

    total_reads = 0
    total_refseq_reads = 0

    with open(summary_csv_path, "w", newline="") as csv_handle, open(normalized_fastq_path, "w") as norm_fastq_handle:

        writer = csv.writer(csv_handle)

        writer.writerow([
            "RunName",
            "trace_name",
            "trace_len",
            "was_rc",
            "extracted_read_start",
            "extracted_read_end",
            "extracted_refseq_len",
            "extracted_refseq",
            "dnstream_start",
            "dnstream_end",
            "dnstream_seq_len",
            "dnstream_seq",
            "avgQMean_extracted_refseq",
            "avgQMean_dnstream_seq",
            "similarity",
            "alnScore",
            "extracted_refseq_coverage",
            "ref_coverage",
            "aligned_bases",
            "refseq_align_start",
            "refseq_align_end",
            "classify_align",
            "classify_read"
        ])

        for fastq_path in fastq_paths:

            run_name = get_run_id(fastq_path)

            print("\nProcessing file:", fastq_path)

            with gzip.open(fastq_path, "rt") as handle:

                for i, record in enumerate(SeqIO.parse(handle, "fastq")):

                    if TEST_READ_LIMIT is not None and i >= TEST_READ_LIMIT:
                        break

                    total_reads += 1

                    trace_name = record.id.split()[0]
                    seq = str(record.seq)
                    quals = record.letter_annotations["phred_quality"]

                    trace_len = len(seq)

                    # Find longest polyA and polyT in original read
                    polyA = longest_polyA_run(seq)
                    polyT = longest_polyT_run(seq)

                    polyA_len_before = run_len(polyA)
                    polyT_len_before = run_len(polyT)

                    was_rc = 0

                    # Normalize into forward polyA orientation.
                    # If polyT run is stronger than polyA run, reverse-complement.
                    if (not polyA) and polyT:
                        seq = reverse_complement(seq)
                        quals = rc_quals(quals)
                        was_rc = 1

                        # Recalculate after RC
                        polyA = longest_polyA_run(seq)

                    polyA_len = run_len(polyA)

                    if polyA:
                        polyA_start, polyA_end = polyA
                    else:
                        polyA_start = polyA_end = 0

                    polyA_found = polyA is not None

                    # Alignment-based extracted REFSEQ extraction
                    extracted_refseq_result = locate_extracted_refseq_by_alignment(seq)

                    if extracted_refseq_result:
                        (
                            ref_s,
                            ref_e,
                            extracted_refseq,
                            ref_start,
                            ref_end,
                            sim,
                            score,
                            cov,
                            ref_cov,
                            aligned_bases
                        ) = extracted_refseq_result

                        extracted_refseq_len = len(extracted_refseq)
                        isRef = 1
                        total_refseq_reads += 1

                    else:
                        ref_s = ref_e = 0
                        extracted_refseq = ""
                        extracted_refseq_len = 0

                        ref_start = ref_end = -1
                        sim = score = cov = ref_cov = 0.0
                        aligned_bases = 0
                        isRef = 0

                    # Extract downstream sequence directly from extracted end in read
                    if extracted_refseq:
                        seq_s = ref_e
                        seq_e = min(len(seq), ref_e + ext_seq_len)
                        dnstream_seq = seq[seq_s:seq_e]
                        dnstream_seq_len = len(dnstream_seq)
                    else:
                        seq_s = seq_e = 0
                        dnstream_seq = ""
                        dnstream_seq_len = 0

                    # Quality metrics
                    avgQMean_extracted_refseq = (
                        mean_q(quals, ref_s, ref_e) if extracted_refseq else 0.0
                    )

                    avgQMean_dnstream_seq = (
                        mean_q(quals, seq_s, seq_e) if dnstream_seq else 0.0
                    )

                    # Classify alignment and read
                    category_align = classify_align(ref_start, ref_end)

                    category_read = classify_read(
                        polyA_found=polyA_found,
                        ref_match=(isRef == 1)
                    )

                    # Temporary debug prints for first few reads
                    if TEST_READ_LIMIT is not None:
                        print("\n----- REAL READ TEST -----")
                        print("Read:", record.id)
                        print("Seq len:", len(seq))
                        print("was_rc:", was_rc)
                        print("extracted_refseq_len:", extracted_refseq_len)
                        print("dnstream_seq_len:", dnstream_seq_len)
                        print("sim:", sim)
                        print("score:", score)
                        print("extracted_refseq_coverage:", cov)
                        print("ref_coverage:", ref_cov)
                        print("aligned_bases:", aligned_bases)
                        print("ref_start:", ref_start)
                        print("ref_end:", ref_end)
                        print("category_align:", category_align)
                        print("category_read:", category_read)

                    if extracted_refseq_len > 0:
                        writer.writerow([
                            run_name,
                            trace_name,
                            trace_len,
                            was_rc,
                            ref_s,
                            ref_e,
                            extracted_refseq_len,
                            extracted_refseq,
                            seq_s,
                            seq_e,
                            dnstream_seq_len,
                            dnstream_seq,
                            avgQMean_extracted_refseq,
                            avgQMean_dnstream_seq,
                            sim,
                            score,
                            cov,
                            ref_cov,
                            aligned_bases,
                            ref_start,
                            ref_end,
                            category_align,
                            category_read
                        ])

                    # Write normalized sequence to FASTQ
                    record.seq = Seq(seq)
                    record.letter_annotations = {}
                    record.letter_annotations["phred_quality"] = quals
                    SeqIO.write(record, norm_fastq_handle, "fastq")

    print("\nDone.")
    print("Total reads processed:", total_reads)
    print("Total REFSEQ-like reads:", total_refseq_reads)
    print("Summary CSV:", summary_csv_path)
    print("Normalized FASTQ:", normalized_fastq_path)


if __name__ == "__main__":
    main()