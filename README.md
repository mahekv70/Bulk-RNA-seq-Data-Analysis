# Bulk-RNA-seq-Data-Analysis

This repository contains a Python-based workflow for comparative analysis
of bulk RNA-seq differential expression data, with a focus on assessing
the reproducibility and concordance of log2 fold-change estimates across
experimental conditions.

The analysis is designed for publication-quality visualization and
statistical comparison of RNA-seq datasets generated using Cuffdiff or
similar differential expression pipelines.

---

## Repository Contents

| File | Description |
|------|-------------|
| `correlation_log2fc.py` | Computes log2 fold-change correlations between two RNA-seq conditions, applies biologically motivated filters, calculates Pearson and Spearman correlations, and generates publication-quality scatter plots. |

---

## Analysis Overview

The workflow performs the following steps:

- Loads gene-level differential expression tables (Excel or TSV)
- Filters genes based on expression thresholds and annotation
- Computes log2 fold-changes with pseudocount stabilization
- Applies percentile-based winsorization to reduce outlier effects
- Calculates Pearson and Spearman correlation coefficients
- Outputs a filtered results table and a vector-quality PDF scatter plot

---

## Dependencies

### Python packages

```bash
pip install pandas numpy matplotlib scipy
