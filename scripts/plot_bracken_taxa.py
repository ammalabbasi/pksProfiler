#!/usr/bin/env python3
import argparse
from pathlib import Path
import re

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def read_bracken_report(path: Path) -> pd.DataFrame:
    if (not path.exists()) or path.stat().st_size == 0:
        return pd.DataFrame(columns=["taxon", "reads"])

    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.DataFrame(columns=["taxon", "reads"])

    if df.empty:
        return pd.DataFrame(columns=["taxon", "reads"])

    # Bracken typical columns include 'name' and 'new_est_reads'
    col_lower = {c.lower(): c for c in df.columns}
    taxon_col = col_lower.get("name", df.columns[0])

    if "new_est_reads" in col_lower:
        reads_col = col_lower["new_est_reads"]
    else:
        # fallback: find a column that looks numeric
        reads_col = None
        for c in reversed(df.columns):
            s = pd.to_numeric(df[c], errors="coerce")
            if s.notna().any():
                reads_col = c
                break
        if reads_col is None:
            return pd.DataFrame(columns=["taxon", "reads"])

    out = pd.DataFrame({
        "taxon": df[taxon_col].astype(str),
        "reads": pd.to_numeric(df[reads_col], errors="coerce"),
    })
    out = out.dropna(subset=["taxon", "reads"])
    out = out[out["reads"] > 0]
    return out


def barplot(ax, df: pd.DataFrame, title: str, top_n: int):
    ax.set_title(title)
    if df.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "No taxa (empty or 0 reads)", ha="center", va="center")
        return

    df = df.sort_values("reads", ascending=False)
    if top_n > 0 and len(df) > top_n:
        df = df.iloc[:top_n].copy()

    ax.bar(df["taxon"], df["reads"])
    ax.set_ylabel("Estimated PKS-island reads")
    ax.set_xlabel("Taxon")
    ax.tick_params(axis="x", labelrotation=60)
    ax.margins(x=0.01)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, help="Sample ID for titles")
    ap.add_argument("--genus", required=True, help="Bracken genus report (*.bracken.G.report.txt)")
    ap.add_argument("--species", required=True, help="Bracken species report (*.bracken.S.report.txt)")
    ap.add_argument("--out", required=True, help="Output PDF")
    ap.add_argument("--top", type=int, default=20, help="Top N taxa to plot (0 = all)")
    args = ap.parse_args()

    sample = args.sample
    genus_path = Path(args.genus)
    species_path = Path(args.species)

    gdf = read_bracken_report(genus_path)
    sdf = read_bracken_report(species_path)

    out_pdf = Path(args.out)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(out_pdf) as pdf:
        # Genus page
        fig, ax = plt.subplots(figsize=(12, 6))
        barplot(ax, gdf, f"{sample} — Bracken Genus (top {args.top if args.top>0 else 'all'})", args.top)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        # Species page
        fig, ax = plt.subplots(figsize=(12, 6))
        barplot(ax, sdf, f"{sample} — Bracken Species (top {args.top if args.top>0 else 'all'})", args.top)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


if __name__ == "__main__":
    main()
