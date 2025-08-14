#!/usr/bin/env python3
import argparse, os
import numpy as np
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="Allele Frequency Calculator")
    p.add_argument("--genotypes", required=True, help="TSV: samples x SNPs (0/1/2 or NA)")
    p.add_argument("--snp-annot", help="TSV: SNP annotation (optional; must contain snp_id)")
    p.add_argument("--sample-subset", help="TSV (optional) with a column 'sample_id' to keep")
    p.add_argument("--missing-code", default="NA", help="Missing value code (default 'NA')")
    p.add_argument("--out", required=True, help="Output TSV (allele frequencies)")
    return p.parse_args()

def main():
    a = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)

    geno = pd.read_csv(a.genotypes, sep="\t", index_col=0, dtype=str)
    geno = geno.replace(a.missing_code, np.nan).apply(pd.to_numeric, errors="coerce")

    if a.sample_subset:
        keep = pd.read_csv(a.sample_subset, sep="\t")
        if "sample_id" not in keep.columns:
            raise SystemExit("sample subset must have a 'sample_id' column")
        keep_ids = set(keep["sample_id"].astype(str))
        geno = geno.loc[geno.index.astype(str).isin(keep_ids), :]

    rows = []
    for snp, col in geno.items():
        col = col.astype("float32")
        n0 = int((col == 0).sum())
        n1 = int((col == 1).sum())
        n2 = int((col == 2).sum())
        n_non_missing = int(col.notna().sum())
        maf = np.nan
        if n_non_missing >= 10:
            p = col.dropna().mean() / 2.0
            maf = float(min(p, 1 - p))
        rows.append({"snp_id": snp, "n0": n0, "n1": n1, "n2": n2,
                     "n_non_missing": n_non_missing, "maf": maf})

    out = pd.DataFrame(rows)
    if a.snp_annot:
        annot = pd.read_csv(a.snp_annot, sep="\t", dtype=str)
        if "snp_id" not in annot.columns:
            raise SystemExit("snp_annotation must contain 'snp_id'")
        out = annot.merge(out, on="snp_id", how="left")

    out.to_csv(a.out, sep="\t", index=False)
    print(f"Done. Wrote {a.out}")

if __name__ == "__main__":
    main()
