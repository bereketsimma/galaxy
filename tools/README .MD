# GWAS Galaxy Tools

This package contains 5 Galaxy tools for basic Genome-Wide Association Study (GWAS) analysis.

## Tools Included

### Tool 1 – SNP & Sample QC
**Purpose:** Filter SNPs and samples based on missingness and minor allele frequency (MAF).  
**Inputs:**
- `genotypes.tsv` (samples x SNPs)
- `snp_annotation.tsv`
- `phenotypes_covariates.tsv` (optional)
- Parameters: MAF threshold, SNP missingness threshold, sample missingness threshold, missing value code  
**Outputs:**
- Filtered genotypes (`.genotypes.filtered.tsv`)
- Filtered SNP annotation (`.snp_annotation.filtered.tsv`)
- Kept sample list (`.kept_samples.tsv`)
- QC report (`.qc_report.txt`)

---

### Tool 2 – Allele Frequency Calculator
**Purpose:** Compute genotype counts and per-SNP MAF.  
**Inputs:**
- `genotypes.tsv`
- `snp_annotation.tsv` (optional)
- Sample subset (optional)  
**Outputs:**
- TSV with SNP ID, chrom, pos, counts for 0/1/2, number of non-missing samples, MAF.

---

### Tool 3 – Association Test (Logistic Regression)
**Purpose:** Fit logistic regression for each SNP: `phenotype ~ genotype + covariates`.  
**Inputs:**
- `genotypes.tsv`
- `phenotypes_covariates.tsv`
- `snp_annotation.tsv` (optional)
- Covariates (comma-separated)  
**Outputs:**
- Association results table with SNP info, beta, SE, z, p-value, MAF.

---

### Tool 4 – Manhattan Plot Generator
**Purpose:** Create a Manhattan plot and list top hits.  
**Inputs:**
- Association results (Tool 3 output)
- `snp_annotation.tsv` (optional)
- Parameters: p-value threshold, number of top hits  
**Outputs:**
- PNG plot
- TSV with top N hits.

---

### Tool 5 – Windowed LD Calculator
**Purpose:** Compute r² between SNPs within a genomic window and plot heatmap.  
**Inputs:**
- `genotypes.tsv`
- `snp_annotation.tsv`
- Focal SNP ID or chromosome + position  
**Outputs:**
- LD matrix TSV
- LD heatmap PNG.

---

## Dependencies
- Python 3.8+
- pandas
- numpy
- scipy
- matplotlib
- statsmodels (optional, for logistic regression)

---

## How to Run Locally
Example for Tool 1:
```bash
python tools/tool1_qc/qc.py \
  --genotypes gwas_data/genotypes.tsv \
  --snp-annot gwas_data/snp_annotation.tsv \
  --phenos gwas_data/phenotypes_covariates.tsv \
  --maf 0.01 \
  --snp-miss 0.05 \
  --sample-miss 0.10 \
  --out-prefix sample_outputs/qc
