#!/usr/bin/env python3
import argparse, os, numpy as np, pandas as pd
import statsmodels.api as sm

def parse_args():
    p = argparse.ArgumentParser(description="Per-SNP logistic regression: phenotype ~ genotype + covariates")
    p.add_argument("--genotypes", required=True, help="TSV: samples x SNPs (0/1/2 or NA)")
    p.add_argument("--phenos", required=True, help="TSV: must contain sample_id, phenotype, and covariates")
    p.add_argument("--snp-annot", help="TSV (optional): contains snp_id, chrom, pos, ...")
    p.add_argument("--covariates", default="", help="Comma-separated covariate names, e.g., age,sex")
    p.add_argument("--min-maf", type=float, default=0.01, help="Min MAF to include (default 0.01)")
    p.add_argument("--missing-code", default="NA", help="Missing value code (default 'NA')")
    p.add_argument("--snp-limit", type=int, help="Optional: limit number of SNPs for quick runs")
    p.add_argument("--out", required=True, help="Output TSV for association results")
    return p.parse_args()

def fit_one(y, X):
    try:
        model = sm.Logit(y, X, missing='drop')
        res = model.fit(disp=False, maxiter=100)
        beta = float(res.params.get('geno', np.nan))
        se   = float(res.bse.get('geno', np.nan))
        z    = float(beta / se) if se and not np.isnan(se) else np.nan
        pval = float(res.pvalues.get('geno', np.nan))
        return beta, se, z, pval
    except Exception:
        return np.nan, np.nan, np.nan, np.nan

def main():
    a = parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)

    geno = pd.read_csv(a.genotypes, sep="\t", index_col=0, dtype=str)
    geno = geno.replace(a.missing_code, np.nan).apply(pd.to_numeric, errors="coerce")

    ph = pd.read_csv(a.phenos, sep="\t")
    if "sample_id" not in ph.columns or "phenotype" not in ph.columns:
        raise SystemExit("phenotypes_covariates.tsv must have 'sample_id' and 'phenotype'")
    ph.index = ph["sample_id"].astype(str)

    # design matrix base (covariates)
    covs = [c.strip() for c in a.covariates.split(",") if c.strip()]
    X_base = pd.DataFrame(index=ph.index)
    for c in covs:
        if c not in ph.columns:
            raise SystemExit(f"Missing covariate in phenotypes file: {c}")
        X_base[c] = pd.to_numeric(ph[c], errors="coerce")
    X_base = sm.add_constant(X_base, has_constant='add')

    snps = list(geno.columns)
    if a.snp_limit is not None:
        snps = snps[:a.snp_limit]

    rows = []
    for snp in snps:
        g = pd.to_numeric(geno[snp], errors="coerce")
        n_non_missing = g.notna().sum()
        if n_non_missing < 20:
            continue
        p = g.dropna().mean() / 2.0
        maf = float(min(p, 1 - p))
        if maf < a.min_maf:
            continue

        X = X_base.copy()
        X["geno"] = g
        y = pd.to_numeric(ph.loc[X.index, "phenotype"], errors="coerce")

        beta, se, z, pval = fit_one(y, X)
        rows.append({"snp_id": snp, "beta": beta, "se": se, "z": z, "p_value": pval, "maf": maf})

    out = pd.DataFrame(rows)

    if a.snp_annot:
        annot = pd.read_csv(a.snp_annot, sep="\t", dtype=str)
        if "snp_id" not in annot.columns:
            raise SystemExit("snp_annotation must contain 'snp_id'")
        out = annot.merge(out, on="snp_id", how="right")

    out.to_csv(a.out, sep="\t", index=False)
    print(f"Done. Wrote {a.out}")

if __name__ == "__main__":
    main()
