"""
06_downsampling.py
 
Sequencing-depth downsampling sensitivity analysis. Sequencing depth differed
between compartments, which could in principle inflate gene detection in the
more deeply sequenced samples. The per-cell depth of each compartment is first
quantified; cells above the median depth of the shallowest compartment are then
downsampled to that depth with sc.pp.downsample_counts (random_state=0),
re-normalised and re-scored; and the cross-compartment Friedman tests are
repeated on the downsampled scores.
 
The procedure is specified in Section 2.7 of the dissertation and its results
are reported in Section 4.5.
 
Dissertation sections    : 2.7, 4.5
Supplementary Appendix   : A9
 
Prerequisites
-------------
Depends on 03 (the signature definitions are re-used) and, through it, on 02.
Continues from a previous script in the same Python session; do not re-load the
object. Only the imports and the metadata column names are repeated below.
"""
 
# --- Metadata column names (repeated from 01 for readability) ----------------
PATIENT_COL  = "Patient_ID"
SAMPLE_COL   = "Sample"
TISSUE_COL   = "Tissue"
STAGE_COL    = "StageGroup"
STAGESEP_COL = "StageSep"
CELLTYPE_COL = "cell_type_with_cluster"
CLUSTER_COL  = "leiden_0.5"
 
 
# ============================================================================
# A9. Sensitivity analysis: sequencing-depth downsampling
# ============================================================================
 
import scanpy as sc, numpy as np, pandas as pd
from scipy.stats import friedmanchisquare
from statsmodels.stats.multitest import multipletests
 
sig_cols = ["CD14_monocyte_score", "CD16_monocyte_score",
            "Macrophage_complement_score", "Disease_associated_macrophage_score"]
tissues = ["PBMC", "LIVER", "VAT", "SAT"]
 
sigs = {
    "CD14_monocyte_score":                 ["CD14", "FCN1", "S100A8",
                                            "S100A9", "VCAN", "CCR2"],
    "CD16_monocyte_score":                 ["FCGR3A", "MS4A7", "LILRB1", "CX3CR1"],
    "Macrophage_complement_score":         ["C1QA", "C1QB", "C1QC", "APOE", "APOC1"],
    "Disease_associated_macrophage_score": ["TREM2", "CD9", "SPP1",
                                            "GPNMB", "LGALS3"],
}
 
# --- 1. Quantify the depth difference between compartments --------------------
raw = adata.layers["counts_RNA"]
tot = np.asarray(raw.sum(axis=1)).ravel()
depth = pd.DataFrame({"Tissue": adata.obs[TISSUE_COL].values, "counts": tot})
print("=== Total counts per cell, by compartment ===")
print(depth.groupby("Tissue")["counts"].describe()[["25%", "50%", "75%"]]
      .round(0).to_string())
 
target = int(depth.groupby("Tissue")["counts"].median().min())
print(f"\nDownsampling target (median depth of shallowest compartment): {target}")
 
# --- 2. Downsample to a common depth, re-normalise and re-score ---------------
ad_ds = adata.copy()
ad_ds.X = ad_ds.layers["counts_RNA"].copy()
sc.pp.downsample_counts(ad_ds, counts_per_cell=target, random_state=0)
sc.pp.normalize_total(ad_ds, target_sum=1e4)
sc.pp.log1p(ad_ds)
 
for name, genes in sigs.items():
    sc.tl.score_genes(ad_ds, gene_list=[g for g in genes if g in ad_ds.var_names],
                      score_name=name + "_ds", use_raw=False)
 
# --- 3. Repeat the Friedman tests on the downsampled scores -------------------
ds_cols = [c + "_ds" for c in sig_cols]
pt_ds = (ad_ds.obs.groupby([PATIENT_COL, TISSUE_COL], observed=True)[ds_cols]
         .mean().reset_index())
 
rows = []
for c, s in zip(sig_cols, ds_cols):
    wide = pt_ds.pivot(index=PATIENT_COL, columns=TISSUE_COL,
                       values=s)[tissues].dropna()
    chi2, p = friedmanchisquare(*[wide[t].values for t in tissues])
    means = wide.mean().round(2).to_dict()
    rows.append({"signature": c, "chi2": round(chi2, 2), "p_raw": p, **means})
 
res_ds = pd.DataFrame(rows)
res_ds["FDR"] = multipletests(res_ds["p_raw"], method="fdr_bh")[1]
print("\n=== Friedman results after downsampling ===")
print(res_ds.to_string(index=False))
res_ds.to_csv(OUTPUT_DIR / "downsampled_friedman.csv", index=False)
 
