#!/usr/bin/env python3
import argparse, sys, os
import numpy as np
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(description="SNP & Sample QC")
    p.add_argument("--genotypes", required=True, help="TSV: samples x SNPs, genotypes 0/1/2 or missing")
    p.add_argument("--snp-annot", required=True, help="TSV: SNP annotation (snp_id, chrom, pos, ...)")
    p.add_argument("--samples", help="TSV: phenotypes_covariates with column 'sample_id' (optional)")
    p.add_argument("--maf", type=float, default=0.01, help="MAF threshold (default 0.01)")
    p.add_argument("--snp-miss", type=float, default=0.05, help="SNP missingness threshold (default 0.05)")
    p.add_argument("--sample-miss", type=float, default=0.10, help="Sample missingness threshold (default 0.10)")
    p.add_argument("--missing-code", default="NA", help="Missing value code in input (default 'NA')")
    p.add_argument("--out-prefix", required=True, help="Prefix for outputs, e.g., 'qc'")
    return p.parse_args()

def coerce_numeric(df, missing_code):
    if missing_code is not None:
        df = df.replace(missing_code, np.nan)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.where(df.isin([0,1,2]))
    return df.astype("float32")

def main():
    a = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.out_prefix)) or ".", exist_ok=True)

    geno = pd.read_csv(a.genotypes, sep="\t", header=0, index_col=0, dtype=str)
    geno = coerce_numeric(geno, a.missing_code)
    snp_annot = pd.read_csv(a.snp_annot, sep="\t", header=0, dtype=str)
    if "snp_id" not in snp_annot.columns:
        sys.exit("snp_annotation.tsv must have a 'snp_id' column")

    n_samples_before, n_snps_before = geno.shape

    if geno.index.duplicated().any():
        geno = geno[~geno.index.duplicated(keep="first")]
    if geno.columns.duplicated().any():
        geno = geno.loc[:, ~geno.columns.duplicated(keep="first")]

    # 1) Sample missingness
    miss_sample = geno.isna().mean(axis=1)
    keep_samples = miss_sample <= a.sample_miss
    dropped_samples = int((~keep_samples).sum())
    geno = geno.loc[keep_samples]
    kept_sample_ids = geno.index.to_series()

    # 2) SNP missingness
    miss_snp = geno.isna().mean(axis=0)
    keep_snps_miss = miss_snp <= a.snp_miss
    dropped_snps_miss = int((~keep_snps_miss).sum())
    geno = geno.loc[:, keep_snps_miss]

    # 3) MAF
    nonmissing_n = geno.count(axis=0)
    valid = nonmissing_n > 0
    alt_sum = geno.sum(axis=0, skipna=True)
    p = pd.Series(np.nan, index=geno.columns)
    p.loc[valid] = (alt_sum[valid] / (2.0 * nonmissing_n[valid]))
    maf = pd.Series(np.minimum(p, 1 - p), index=geno.columns)
    keep_snps_maf = (maf >= a.maf).fillna(False)
    dropped_snps_maf = int((~keep_snps_maf).sum())
    geno = geno.loc[:, keep_snps_maf]

    # Annot kept + recomputed MAF
    snp_annot = snp_annot.set_index("snp_id")
    annot_kept = snp_annot.reindex(geno.columns).dropna(how="all")
    annot_kept["maf_recomputed"] = maf.reindex(geno.columns).values

    # Outputs
    g_path = f"{a.out_prefix}.genotypes.filtered.tsv"
    s_path = f"{a.out_prefix}.snp_annotation.filtered.tsv"
    k_path = f"{a.out_prefix}.kept_samples.tsv"
    r_path = f"{a.out_prefix}.qc_report.txt"

    geno.to_csv(g_path, sep="\t", header=True, index=True)
    annot_kept.reset_index().to_csv(s_path, sep="\t", index=False)
    kept_sample_ids.to_frame(name="sample_id").to_csv(k_path, sep="\t", index=False)

    with open(r_path, "w") as f:
        f.write("SNP & Sample QC Report\n")
        f.write("======================\n\n")
        f.write("Inputs:\n")
        f.write(f"  genotypes: {a.genotypes}\n")
        f.write(f"  snp_annotation: {a.snp_annot}\n")
        if a.samples: f.write(f"  samples: {a.samples}\n")
        f.write("\nThresholds:\n")
        f.write(f"  MAF >= {a.maf}\n")
        f.write(f"  SNP missingness <= {a.snp_miss}\n")
        f.write(f"  Sample missingness <= {a.sample_miss}\n")
        f.write(f"  Missing value code: {a.missing_code}\n\n")
        f.write("Counts:\n")
        f.write(f"  Samples: before={n_samples_before}, kept={geno.shape[0]}, dropped={dropped_samples}\n")
        f.write(f"  SNPs:    before={n_snps_before}, "
                f"dropped_by_snp_missingness={dropped_snps_miss}, "
                f"dropped_by_maf={dropped_snps_maf}, kept={geno.shape[1]}\n\n")
        f.write("Distributions (before filtering):\n")
        f.write(f"  median sample missingness: {miss_sample.median():.4f}\n")
        f.write(f"  median SNP missingness: {miss_snp.median():.4f}\n")

    print("Done.")
    print(f"  {g_path}")
    print(f"  {s_path}")
    print(f"  {k_path}")
    print(f"  {r_path}")

if __name__ == "__main__":
    main()
