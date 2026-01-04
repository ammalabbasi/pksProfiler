#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Optional
import pandas as pd

def build_matrix(files, suffix_to_strip: Optional[str] = None) -> pd.DataFrame:
    files = sorted([Path(f) for f in files])
    if not files:
        raise SystemExit("No input files provided.")

    columns = []
    for f in files:
        name = f.name
        sample = name
        if suffix_to_strip and name.endswith(suffix_to_strip):
            sample = name[:-len(suffix_to_strip)]
        else:
            sample = f.stem

        df = pd.read_csv(f, sep="\t", usecols=["Gene", "Count"], dtype={"Gene": str, "Count": float})
        df = df.dropna(subset=["Gene"])
        s = df.groupby("Gene", as_index=True)["Count"].sum()
        s.name = sample
        columns.append(s)

    matrix = pd.concat(columns, axis=1).fillna(0)

    if (matrix % 1 == 0).all().all():
        matrix = matrix.astype(int)

    matrix = matrix.reindex(sorted(matrix.columns), axis=1)
    return matrix

def main():
    ap = argparse.ArgumentParser(description="Build a Gene x Sample count matrix from per-sample HMM count TSVs.")
    ap.add_argument("--inputs", nargs="+", required=True, help="Input TSV files with columns: Gene, Count")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--strip-suffix", default=None, help="Suffix to strip from filename to form sample name")
    args = ap.parse_args()

    m = build_matrix(args.inputs, suffix_to_strip=args.strip_suffix)
    m.to_csv(args.out, sep="\t", index=True)

if __name__ == "__main__":
    main()
