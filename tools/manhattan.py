#!/usr/bin/env python3
import argparse, os, numpy as np, pandas as pd
import matplotlib.pyplot as plt

def parse_args():
    p = argparse.ArgumentParser(description="Manhattan Plot Generator")
    p.add_argument("--assoc", required=True, help="TSV: association results with snp_id, p_value (+ chrom,pos if available)")
    p.add_argument("--snp-annot", help="TSV: optional, to add chrom/pos if assoc lacks them")
    p.add_argument("--p-threshold", type=float, default=5e-8, help="Genome-wide threshold (default 5e-8)")
    p.add_argument("--top-n", type=int, default=20, help="Top N hits to output (default 20)")
    p.add_argument("--out-prefix", required=True, help="Prefix for PNG and top hits TSV")
    return p.parse_args()

def chr_key(x):
    s = str(x).lower().replace("chr", "")
    try:
        return int(s)
    except:
        return 1_000_000_000

def main():
    a = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)) or ".", exist_ok=True)

    df = pd.read_csv(a.assoc, sep="\t")
    if "snp_id" not in df.columns or "p_value" not in df.columns:
        raise SystemExit("Association file must have 'snp_id' and 'p_value'")

    if (("chrom" not in df.columns) or ("pos" not in df.columns)) and a.snp_annot:
        annot = pd.read_csv(a.snp_annot, sep="\t", dtype=str)
        if "snp_id" not in annot.columns:
            raise SystemExit("snp_annotation must contain 'snp_id'")
        df = annot.merge(df, on="snp_id", how="right")

    # Drop rows without p_value or position info
    df = df.dropna(subset=["p_value"]).copy()
    if "chrom" not in df.columns or "pos" not in df.columns:
        raise SystemExit("Need 'chrom' and 'pos' columns (provide --snp-annot if assoc lacks them)")

    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"])
    df["chr_order"] = df["chrom"].map(chr_key)
    df = df.sort_values(["chr_order", "pos"])
    df["neglog10p"] = -np.log10(pd.to_numeric(df["p_value"], errors="coerce"))

    # Build cumulative x position by chromosome
    cum = 0
    xs = []
    ticks = []
    ticklabels = []
    for chrom, sub in df.groupby("chrom", sort=False):
        start = cum
        xvals = start + sub["pos"].values
        xs.append(pd.Series(xvals, index=sub.index))
        mid = start + (sub["pos"].min() + sub["pos"].max()) / 2.0
        ticks.append(mid)
        ticklabels.append(str(chrom))
        cum = xvals.max() + 1
    df["x"] = pd.concat(xs).sort_index()

    # Plot
    plt.figure(figsize=(11, 4))
    plt.scatter(df["x"], df["neglog10p"], s=6)
    plt.axhline(-np.log10(a.p_threshold), linestyle="--")
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(p)")
    plt.xticks(ticks, ticklabels, rotation=0)
    plt.tight_layout()
    png = f"{a.out_prefix}.png"
    plt.savefig(png, dpi=160)
    plt.close()

    # Top hits
    keep_cols = [c for c in ["snp_id", "chrom", "pos", "p_value", "beta", "se", "z", "maf"] if c in df.columns]
    top = df.nsmallest(a.top_n, "p_value")[keep_cols]
    top_path = f"{a.out_prefix}.top_hits.tsv"
    top.to_csv(top_path, sep="\t", index=False)

    print("Done.")
    print(f"  {png}")
    print(f"  {top_path}")

if __name__ == "__main__":
    main()
