
#This code is fining only the Moderna NGS runs that containing the polyA stretch =75
import os
import csv
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

# ------------------- CONFIG -------------------

# Network root where runs live
ROOT = r"O:\MiSeq Data"

# Limit to top-level run folders starting with any of these prefixes (old MiSeq vs new i100)
RUN_PREFIXES = ("25", "26", "2025", "2026")

# Where manifests live under each run folder
FASTQ_DIRNAME = "FASTQ Files"

# Output report file
OUTPUT_CSV = "Moderna_hits3.csv"

# Column indices for headerless manifest.csv (0-based):
#   A=0, B=1, C=2, D=3, E=4, F=5, ...
# D=PrimaryName, E=SecondaryName, F=Reference.
# WellPosition is LOCKED to column A (index 0) and written as "WellPositionName".
WELLPOSITION_COL   = 0  # A (LOCKED)
PRIMARY_COL   = 3  # D
SECONDARY_COL = 4  # E
REF_COL       = 5  # F

# Detection threshold for homopolymers
POLY_MIN = 75

# What to report:
#   Default: Only poly-A runs (>= POLY_MIN). Set True to also include poly-T runs.
INCLUDE_T_RUNS = False

# Compile regex(es) for detection ANYWHERE in the sequence
POLY_A_ANY_RE = re.compile(rf"(A{{{POLY_MIN},}})", re.IGNORECASE)
POLY_T_ANY_RE = re.compile(rf"(T{{{POLY_MIN},}})", re.IGNORECASE) if INCLUDE_T_RUNS else None

# Full IUPAC DNA set (uppercase) used for trailing trim
IUPAC_DNA = "ACGTRYSWKMBDHVN"

# Cleaning options for the reference string:
# - Remove surrounding square brackets like [SEQUENCE]
# - Uppercase and remove spaces
# - Trim trailing **non-IUPAC** characters (punctuation/annotations at the end)
TRIM_TRAILING_NON_IUPAC_RE = re.compile(rf"[^{IUPAC_DNA}]+$", re.IGNORECASE)


# --- Moderna filtering ---

# Patterns we consider "Moderna" (case-insensitive).
# Includes explicit 'moderna', product aliases, and common program tags starting with "mod".
MODERNA_PATTERNS = [
    r"moderna",                # substring anywhere (e.g., ModernaPilot, SRT_ModernaA1)
    r"mrna[-_ ]?1273",         # vaccine code
    r"spikevax",               # brand
    r"(^|[^a-z])mod[a-z]+",    # mod* program tags: modvec, modrework, modsdm, modpilot, etc.
]
MODERNA_REGEXES = [re.compile(p, re.IGNORECASE) for p in MODERNA_PATTERNS]

def is_moderna(primary: str, secondary: str) -> tuple[bool, str]:
    """
    Return (True, matched_pattern) if either PrimaryName or SecondaryName
    indicates 'Moderna' by any of MODERNA_REGEXES; otherwise (False, "").
    """
    text = f"{str(primary)} {str(secondary)}"
    for rx in MODERNA_REGEXES:
        if rx.search(text):
            return True, rx.pattern
    return False, ""


# ------------------- HELPERS -------------------

def file_times(path: Path) -> Tuple[float, float]:
    """
    Return (birth_time, mod_time) in seconds since epoch.
    On Windows NTFS, st_ctime is creation time; on some SMB shares it may be server-defined.
    We attempt st_birthtime if present; else fallback to st_ctime; always have mtime.
    """
    try:
        st = path.stat()
        birth = getattr(st, "st_birthtime", None)
        if birth is None:
            birth = st.st_ctime  # Windows creation time; on POSIX this is ctime (metadata change)
        return float(birth), float(st.st_mtime)
    except Exception:
        # If anything fails, return "huge" times so such files are deprioritized
        return float("inf"), float("inf")

def choose_oldest_csv(csv_files: List[Path]) -> Optional[Path]:
    """
    Choose the original CSV as the oldest by creation time (fallback: mtime).
    If multiple tie within 1 second, prefer a name containing 'manifest'.
    """
    if not csv_files:
        return None

    scored = []
    for f in csv_files:
        b, m = file_times(f)
        scored.append((f, b, m))

    # Sort by (birth_time, mod_time) ascending (oldest first)
    scored.sort(key=lambda x: (x[1], x[2]))

    oldest_birth = scored[0][1]
    oldest_mod = scored[0][2]
    TOL = 1.0  # seconds tolerance for ties
    tied = [t for t in scored if abs(t[1] - oldest_birth) <= TOL and abs(t[2] - oldest_mod) <= TOL]

    if len(tied) == 1:
        return tied[0][0]

    # Tie-breaker: prefer filenames containing 'manifest'
    with_manifest = [t for t in tied if "manifest" in t[0].name.lower()]
    if with_manifest:
        return with_manifest[0][0]

    # Otherwise, pick the first among the tied (still the oldest)
    return tied[0][0]

def find_csvs_in_fastq(run_dir: Path) -> List[Path]:
    fq = run_dir / FASTQ_DIRNAME
    if not fq.is_dir():
        return []
    return [p for p in fq.iterdir() if p.is_file() and p.suffix.lower() == ".csv"]

def clean_reference(raw: str) -> str:
    """
    Normalize reference sequence:
      - Uppercase
      - Remove spaces
      - Strip surrounding [ ] if present (e.g., "[GGAGG...]")
      - Trim trailing non-IUPAC characters (e.g., "*", ".", adapter markers, etc.)
    """
    if raw is None:
        return ""
    s = raw.strip().upper().replace(" ", "")
    # Strip simple surrounding brackets like [SEQUENCE]
    if s.startswith("[") and s.endswith("]") and len(s) >= 2:
        s = s[1:-1]
    # Trim any trailing non-IUPAC characters at the 3' end
    s = TRIM_TRAILING_NON_IUPAC_RE.sub("", s)
    return s

def find_all_poly_runs(seq: str) -> List[Tuple[str, int, int, int]]:
    """
    Find ALL homopolymer runs of A (and optionally T) with length >= POLY_MIN ANYWHERE.
    Returns a list of tuples:
        (match_type, start_1based, end_1based, run_len)
      where match_type in {"polyA", "polyT"}.
    Coordinates are 1-based, inclusive, relative to the CLEANED sequence.
    """
    hits: List[Tuple[str, int, int, int]] = []
    if not seq:
        return hits

    # Poly-A runs
    for m in POLY_A_ANY_RE.finditer(seq):
        run_len = len(m.group(0))
        start_1 = m.start() + 1
        end_1 = m.end()
        hits.append(("polyA", start_1, end_1, run_len))

    # Optional: Poly-T runs
    if INCLUDE_T_RUNS and POLY_T_ANY_RE is not None:
        for m in POLY_T_ANY_RE.finditer(seq):
            run_len = len(m.group(0))
            start_1 = m.start() + 1
            end_1 = m.end()
            hits.append(("polyT", start_1, end_1, run_len))

    # Sort runs by start coordinate for readability
    hits.sort(key=lambda t: t[1])
    return hits

def primary_matches_moderna(primary: str) -> Tuple[bool, str]:
    """
    Returns (is_moderna, matched_pattern) where matched_pattern is the regex pattern string that matched.
    """
    if not primary:
        return False, ""
    for rx in MODERNA_REGEXES:
        if rx.search(primary):
            return True, rx.pattern
    return False, ""

def csv_row_iter(path: Path):
    """
    Yields rows from CSV, trying UTF-8 variants then Latin-1 if needed.
    """
    for enc in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            with open(path, "r", encoding=enc, newline="") as f:
                for row in csv.reader(f):
                    yield row
            return
        except Exception:
            continue
    # If all encodings failed, re-raise last exception by opening to fail
    with open(path, "r", encoding="utf-8", newline="") as f:
        for row in csv.reader(f):
            yield row

# ------------------- MAIN -------------------

def main():
    print("Column mapping (0-based): WellPosition=A(0) [LOCKED], Primary=D(3), Secondary=E(4), Reference=F(5)")
    print(f"Homopolymer detection: poly-A ANYWHERE with length ≥ {POLY_MIN}"
          f"{' (and poly-T too)' if INCLUDE_T_RUNS else ''}")
    print("Moderna aliases (PrimaryName):", ", ".join(MODERNA_PATTERNS))
    print(f"IUPAC allowed set for trailing trim: {IUPAC_DNA}")

    root = Path(ROOT)
    if not root.is_dir():
        print(f"ERROR: Root path not found: {ROOT}")
        sys.exit(1)

    # Non-recursive: only top-level directories matching prefixes
    top_level_dirs = [d for d in root.iterdir() if d.is_dir()]
    candidate_runs = [d for d in top_level_dirs if d.name.startswith(RUN_PREFIXES)]

    print(f"Top-level entries: {len(top_level_dirs)}")
    print(f"Candidate run folders (prefix-filtered): {len(candidate_runs)}")

    chosen_count = 0
    total_rows = 0
    total_hits = 0

    out_path = Path(OUTPUT_CSV)
    with open(out_path, "w", newline="", encoding="utf-8") as outf:
        writer = csv.writer(outf)
        writer.writerow([
            "RunName", "ManifestFile", "WellPosition", "PrimaryName", "SecondaryName",
            "PrimaryMatch", "Type", "PolyALength", "PolyAStart", "PolyAEnd", "RefLength"
        ])

        for i, run in enumerate(sorted(candidate_runs, key=lambda p: p.name), start=1):
            csv_files = find_csvs_in_fastq(run)
            if not csv_files:
                print(f"[{i}/{len(candidate_runs)}] {run.name}: no CSV files in '{FASTQ_DIRNAME}'")
                continue

            chosen = choose_oldest_csv(csv_files)
            if not chosen:
                print(f"[{i}/{len(candidate_runs)}] {run.name}: could not choose original CSV")
                continue

            chosen_count += 1
            b, m = file_times(chosen)
            print(f"[{i}/{len(candidate_runs)}] {run.name}: ORIGINAL CSV -> {chosen.name} (created={b:.0f}, modified={m:.0f})")

            rows_this = 0
            hits_this = 0

            try:
                for row in csv_row_iter(chosen):
                    if not row:
                        continue
                    rows_this += 1

                    # Extract fields (headerless)
                    primary = row[PRIMARY_COL].strip()   if len(row) > PRIMARY_COL   else ""
                    is_mod, matched_pat = primary_matches_moderna(primary)
                    if not is_mod:
                        continue  # <-- Only process Moderna rows

                    ref_raw = row[REF_COL].strip() if len(row) > REF_COL else ""
                    ref = clean_reference(ref_raw)
                    if not ref:
                        continue

                    # Find ALL qualifying runs (one output row per run)
                    runs = find_all_poly_runs(ref)
                    if not runs:
                        continue

                    wellposition_name = row[WELLPOSITION_COL].strip()   if len(row) > WELLPOSITION_COL   else ""
                    secondary    = row[SECONDARY_COL].strip() if len(row) > SECONDARY_COL else ""

                    for match_type, start_1, end_1, run_len in runs:
                        writer.writerow([
                            run.name, chosen.name, wellposition_name, primary, secondary,
                            matched_pat, match_type, run_len, start_1, end_1, len(ref)
                        ])
                        hits_this += 1

            except Exception as e:
                print(f"[{i}/{len(candidate_runs)}] {run.name}: ERROR reading {chosen.name}: {e}")
                continue

            total_rows += rows_this
            total_hits += hits_this
            print(f"[{i}/{len(candidate_runs)}] {run.name}: rows={rows_this}, hits={hits_this}")

    print("----- SUMMARY -----")
    print(f"Candidate run folders: {len(candidate_runs)}")
    print(f"Original CSVs chosen:  {chosen_count}")
    print(f"Total rows read:       {total_rows}")
    print(f"Total hits written:    {total_hits}")
    print(f"Output CSV:            {out_path.resolve()}")

if __name__ == "__main__":
    main()