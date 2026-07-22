"""
06_downsampling.py

Sequencing-depth sensitivity analysis.

Subcutaneous adipose tissue was sequenced more deeply than the other
compartments, which could in principle inflate the compartmental differences.
To check that the results are not driven by depth, all cells are downsampled to
a common depth (the lowest tissue median), the signatures are recalculated, and
the Friedman test is re-run. If the compartmental differences and their ranking
are unchanged, the findings are robust to sequencing depth.

Input : mono_macro_with_scores.h5ad
Output: Friedman results on the depth-matched data
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare
from statsmodels.stats.multitest import multipletests

# ---------------------------------------------------------------------------
# Paths and settings
# ---------------------------------------------------------------------------
DATA_PATH = "outputs/signatures/mono_macro_with_scores.h5ad"
OUTPUT_DIR = "outputs/downsampling"
os.makedirs(OUTPUT_DIR, exist_ok=True)

adata = sc.read_h5ad(DATA_PATH)

sig_cols = ["CD14_monocyte_score", "CD16_monocyte_score",
            "Macrophage_complement_score", "Disease_associated_macrophage_score"]
tissues = ["PBMC", "LIVER", "VAT", "SAT"]

signatures = {
    "CD14_monocyte_score": ["CD14", "FCN1", "S100A8", "S100A9", "VCAN", "CCR2"],
    "CD16_monocyte_score": ["FCGR3A", "MS4A7", "LILRB1", "CX3CR1"],
    "Macrophage_complement_score": ["C1QA", "C1QB", "C1QC", "APOE", "APOC1"],
    "Disease_associated_macrophage_score": ["TREM2", "CD9", "SPP1", "GPNMB", "LGALS3"],
}

# ---------------------------------------------------------------------------
# 1. Quantify the depth difference between compartments
# ---------------------------------------------------------------------------
total_counts = np.asarray(adata.layers["counts_RNA"].sum(axis=1)).ravel()
depth = pd.DataFrame({"Tissue": adata.obs["Tissue"].values, "counts": total_counts})
print("Median counts per cell by tissue:")
print(depth.groupby("Tissue")["counts"].median().round(0).to_string())

# Downsample target = lowest tissue median
target = int(depth.groupby("Tissue")["counts"].median().min())
print(f"\nDownsampling target: {target} counts per cell")

# ---------------------------------------------------------------------------
# 2. Downsample to the common depth and recompute the signatures
# ---------------------------------------------------------------------------
ad = adata.copy()
ad.X = ad.layers["counts_RNA"].copy()
sc.pp.downsample_counts(ad, counts_per_cell=target, random_state=0)
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)

for name, genes in signatures.items():
    sc.tl.score_genes(ad, gene_list=[g for g in genes if g in ad.var_names],
                      score_name=name + "_ds", use_raw=False)

# ---------------------------------------------------------------------------
# 3. Re-run the Friedman test on the depth-matched scores
# ---------------------------------------------------------------------------
ds_cols = [c + "_ds" for c in sig_cols]
pt = (ad.obs.groupby(["Patient_ID", "Tissue"], observed=True)[ds_cols]
      .mean().reset_index())

rows = []
for c, s in zip(sig_cols, ds_cols):
    wide = pt.pivot(index="Patient_ID", columns="Tissue", values=s)[tissues].dropna()
    chi2, p = friedmanchisquare(*[wide[t].values for t in tissues])
    rows.append({"signature": c, "chi2": round(chi2, 2), "p_raw": p,
                 **wide.mean().round(2).to_dict()})
res = pd.DataFrame(rows)
res["FDR"] = multipletests(res["p_raw"], method="fdr_bh")[1]

print("\nFriedman results after downsampling:")
print(res.to_string(index=False))
res.to_csv(f"{OUTPUT_DIR}/downsampled_friedman.csv", index=False)

print(f"\nDownsampling outputs written to: {OUTPUT_DIR}")
