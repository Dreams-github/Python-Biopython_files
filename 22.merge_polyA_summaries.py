import sys
from pathlib import Path
import pandas as pd

# ---------------- USER SETTINGS ----------------
ROOT_DIR = Path(r"G:\Platform Innovation - Applied Biology\Staff\swapna\MD")
OUT_FILE = ROOT_DIR / "polyA_summary_ALL_runs.csv"

FINAL_COLS = [
    "trace_name","run_id","run_name","samplePrimary","sampleSecondary","wasRC","trace_length",
    "polyA_length","hasFlanking",
    "isReferenceUp","isReferenceDn",
    "avgQMeanUp","avgQMeanDn",
    "simUp","simDn",
    "alnScoreUp","alnScoreDn",
    "coverUp","coverDn",
    "polyASeq","flankUpSeq","flankDnSeq"
]
# ------------------------------------------------

def get_run_id(csv_path: Path) -> str:
    """
    Extract run folder name, e.g.
    MD/250821_M04095_0651_000000000-GT52R/flank_output_xxxx/polyA_summary.csv
    """
    for parent in csv_path.parents:
        if parent.parent == ROOT_DIR:
            return parent.name
    return ""

def main():
    csv_files = list(ROOT_DIR.rglob("polyA_summary_*.csv"))
    if not csv_files:
        print(f"[WARN] No polyA_summary_*.csv found under {ROOT_DIR}")
        sys.exit(0)

    frames = []

    for csv_file in csv_files:
        try:
            df = pd.read_csv(csv_file, dtype=str)
        except Exception as e:
            print(f"[WARN] Skipping {csv_file}: {e}")
            continue

        # Add run_id
        
        df["run_id"] = get_run_id(csv_file)

        # Ensure all required columns exist
        for col in FINAL_COLS:
            if col not in df.columns:
                df[col] = ""

        # Reorder columns exactly
        df = df[FINAL_COLS]
        frames.append(df)

    if not frames:
        print("[WARN] No valid files to merge.")
        sys.exit(0)

    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(OUT_FILE, index=False)

    print(f"[DONE] Merged {len(frames)} runs")
    print(f"[INFO] Output: {OUT_FILE}")
    print(f"[INFO] Rows: {len(combined)}")

if __name__ == "__main__":
    main()