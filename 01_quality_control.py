"""
01_quality_control.py

Quality-control assessment of the monocyte/macrophage AnnData object.

This step inspects the quality-control metrics already present in the object
(counts, detected genes, mitochondrial and haemoglobin content) to verify the
filtering performed upstream by the host laboratory, rather than to define new
thresholds. No cells are removed at this stage.

Input : mono_macro.h5ad  (processed AnnData object from the host laboratory)
Output: QC summary tables and violin/scatter plots
"""

import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# File paths — update before running
# ---------------------------------------------------------------------------
DATA_PATH = "data/mono_macro.h5ad"     # path to the input .h5ad object
OUTPUT_DIR = "outputs/qc"              # folder for QC outputs
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Load the data
# ---------------------------------------------------------------------------
adata = sc.read_h5ad(DATA_PATH)

# QC metrics recorded in the object, and the key grouping variables
qc_cols = ["nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb"]
key_cols = qc_cols + ["Tissue", "Patient_ID", "StageGroup", "cell_type_with_cluster"]

# ---------------------------------------------------------------------------
# 1. Summary tables
# ---------------------------------------------------------------------------
# Overall QC summary (median, quartiles, etc.)
adata.obs[qc_cols].describe().T.to_csv(f"{OUTPUT_DIR}/qc_overall_summary.csv")

# Missing-value check (21 liver cells lack mito/hb values)
adata.obs[key_cols].isna().sum().to_csv(f"{OUTPUT_DIR}/missing_values.csv")

# QC metrics broken down by tissue and by broad cell type
adata.obs.groupby("Tissue")[qc_cols].describe().to_csv(f"{OUTPUT_DIR}/qc_by_tissue.csv")
adata.obs.groupby("cell_type_with_cluster")[qc_cols].describe().to_csv(
    f"{OUTPUT_DIR}/qc_by_cell_type.csv"
)

# ---------------------------------------------------------------------------
# 2. Diagnostic plots
# ---------------------------------------------------------------------------
# Violin plots of QC metrics, grouped by tissue
sc.pl.violin(adata, qc_cols, groupby="Tissue", rotation=45, show=False)
plt.savefig(f"{OUTPUT_DIR}/qc_violin_by_tissue.png", dpi=300, bbox_inches="tight")
plt.close()

# Scatter of counts vs detected genes (expected positive relationship)
sc.pl.scatter(adata, x="nCount_RNA", y="nFeature_RNA", show=False)
plt.savefig(f"{OUTPUT_DIR}/scatter_counts_features.png", dpi=300, bbox_inches="tight")
plt.close()

# Scatter of counts vs mitochondrial percentage
sc.pl.scatter(adata, x="nCount_RNA", y="percent.mt", show=False)
plt.savefig(f"{OUTPUT_DIR}/scatter_counts_mt.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"QC outputs written to: {OUTPUT_DIR}")
