"""
04_doublet_detection.py
 
Applies Scrublet to the raw counts separately within each sequencing run, since
doublets form within a run, and summarises the resulting scores per subcluster.
The number of principal components is adapted to the size of each run, and runs
with fewer than 30 cells are skipped.
 
Dissertation sections    : 2.4
Tables produced          : Table 4
Supplementary Appendix   : A1.6
 
Prerequisites
-------------
Requires the subcluster labels from 02; independent of 03. Run it at any point
after 02. Continues from a previous script in the same Python session; do not
re-load the object. Only the imports and the metadata column names are repeated
below.
"""
 
import scanpy as sc
import pandas as pd
import numpy as np
 
# --- Metadata column names (repeated from 01 for readability) ----------------
PATIENT_COL  = "Patient_ID"
SAMPLE_COL   = "Sample"
TISSUE_COL   = "Tissue"
STAGE_COL    = "StageGroup"
STAGESEP_COL = "StageSep"
CELLTYPE_COL = "cell_type_with_cluster"
CLUSTER_COL  = "leiden_0.5"
 
 
# ============================================================================
# A1.6 Doublet scoring with Scrublet (Table 4)
# ============================================================================
 
adb = adata.copy()
adb.X = adb.layers["counts_RNA"].copy()
 
scores  = pd.Series(np.nan, index=adata.obs_names, dtype=float)
preds   = pd.Series(np.nan, index=adata.obs_names, dtype=float)
skipped = []
 
for s in adb.obs[SAMPLE_COL].unique():
    sub = adb[adb.obs[SAMPLE_COL] == s].copy()
    n = sub.n_obs
    if n < 30:                                   # too few cells to score
        skipped.append((s, n))
        continue
    npc = min(30, n // 3, sub.n_vars - 1)        # adapt n_prin_comps to run size
    npc = max(npc, 2)
    try:
        sc.pp.scrublet(sub, n_prin_comps=npc, random_state=0, verbose=False)
        scores[sub.obs_names] = sub.obs["doublet_score"].values
        preds[sub.obs_names]  = sub.obs["predicted_doublet"].astype(float).values
    except Exception as e:
        skipped.append((s, n))
        print(f"skip {s} (n={n}): {type(e).__name__}")
 
adata.obs["doublet_score"]     = scores.values
adata.obs["predicted_doublet"] = preds.values
 
print("\nSkipped runs:", skipped)
print("Cells scored: %d / %d" % (scores.notna().sum(), adata.n_obs))
 
summ = (adata.obs.dropna(subset=["doublet_score"])
        .groupby(CLUSTER_COL, observed=True)
        .agg(n_cells=("doublet_score", "size"),
             median_score=("doublet_score", "median"),
             mean_score=("doublet_score", "mean"),
             pct_predicted_doublet=("predicted_doublet", lambda x: 100 * np.mean(x)))
        .round(3))
print(summ.to_string())
summ.to_csv(OUTPUT_DIR / "doublet_by_subcluster.csv")
