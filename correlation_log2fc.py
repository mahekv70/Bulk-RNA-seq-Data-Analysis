#!/usr/bin/env python3
"""
Correlation analysis of Cuffdiff log2 fold-changes between conditions.

This script compares log2FC values derived from two Cuffdiff gene
expression tables, applies biologically motivated filters, computes
Pearson and Spearman correlations, and generates publication-quality
scatter plots.

Author: Mahekdeep Kaur
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
# === input files ===
file1 = "/Users/mahek/Desktop/LCRB/projects/rrna_depletion/analysis/cuffdiff/e_coli/5%/gene_exp.annotated.xlsx"   # e.g., NEB
file2 = "/Users/mahek/Desktop/LCRB/projects/rrna_depletion/analysis/cuffdiff/e_coli/10%/gene_exp.annotated.xlsx"   # e.g., 5%

# === where to save outputs ===
OUTPUT_DIR = "/Users/mahek/Desktop/LCRB/projects/rrna_depletion/analysis/cuffdiff/e_coli/correlation/log2fc"  # <-- change me
FILEROOT   = "e_coli_5%_vs_10%"                            # base name
os.makedirs(OUTPUT_DIR, exist_ok=True)

# === knobs ===
PSEUDO = 1e-3         # pseudocount
USE_STATUS_OK = False  # filter to status == "OK" if column exists
MIN_EXPR = 1.0        # require value_1 + value_2 >= MIN_EXPR
P_LO, P_HI = 1, 99    # percentile clip for log2FC (winsorize)

def load_with_filters(fp):
    # choose read function depending on file extension
    if fp.endswith(".xlsx") or fp.endswith(".xls"):
        df = pd.read_excel(fp)
    else:
        df = pd.read_csv(fp, sep="\t")

    gene_col = "gene_id" if "gene_id" in df.columns else "gene"
    cols = [gene_col, "value_1", "value_2"] + \
           (["status"] if "status" in df.columns else []) + \
           (["gene_type"] if "gene_type" in df.columns else [])
    df = df[cols].copy().rename(columns={gene_col: "gene"})

    # --- filters ---
    if USE_STATUS_OK and "status" in df.columns:
        df = df[df["status"] == "OK"]

    if "gene_type" in df.columns:
        df = df[df["gene_type"] == "protein_coding"]

    df["value_1"] = pd.to_numeric(df["value_1"], errors="coerce")
    df["value_2"] = pd.to_numeric(df["value_2"], errors="coerce")
    df = df.dropna(subset=["value_1", "value_2"])
    df = df[(df["value_1"] + df["value_2"]) >= MIN_EXPR]

    df["log2fc"] = np.log2((df["value_2"] + PSEUDO) / (df["value_1"] + PSEUDO))
    return df[["gene", "log2fc"]]

# load both, inner join on common genes
a = load_with_filters(file1).rename(columns={"log2fc":"log2fc_1"})
b = load_with_filters(file2).rename(columns={"log2fc":"log2fc_2"})
merged = a.merge(b, on="gene", how="inner")

# winsorize extremes
low1, high1 = np.percentile(merged["log2fc_1"], [P_LO, P_HI])
low2, high2 = np.percentile(merged["log2fc_2"], [P_LO, P_HI])
merged["log2fc_1"] = merged["log2fc_1"].clip(low1, high1)
merged["log2fc_2"] = merged["log2fc_2"].clip(low2, high2)

# correlations
pearson_r, pearson_p = pearsonr(merged["log2fc_1"], merged["log2fc_2"])
spearman_r, spearman_p = spearmanr(merged["log2fc_1"], merged["log2fc_2"])
print(f"Common genes used after filters: {len(merged):,}")
print(f"Pearson r = {pearson_r:.3f} (p={pearson_p:.3e})")
print(f"Spearman ρ = {spearman_r:.3f} (p={spearman_p:.3e})")

# save the table to your chosen folder
csv_path = os.path.join(OUTPUT_DIR, f"{FILEROOT}.csv")
merged.to_csv(csv_path, index=False)

# plot and SAVE AS PDF to your chosen folder
plt.figure(figsize=(7,6))
plt.scatter(merged["log2fc_1"], merged["log2fc_2"], alpha=0.5, edgecolor="none")
plt.xlabel("5%")
plt.ylabel("10%")
plt.title(
    "Correlation of Cuffdiff log2FC\n"
    f"Pearson r={pearson_r:.2f} | Spearman ρ={spearman_r:.2f}"
)
plt.grid(True, linestyle="--", alpha=0.4)
plt.tight_layout()

pdf_path = os.path.join(OUTPUT_DIR, f"{FILEROOT}.pdf")  # <-- your PDF
plt.savefig(pdf_path, bbox_inches="tight")              # vector PDF
plt.show()

print(f"Saved CSV to: {csv_path}")
print(f"Saved PDF to: {pdf_path}")
