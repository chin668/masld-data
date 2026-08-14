"""
01_quality_control.py
 
Loads the processed AnnData object, defines the metadata column names used
throughout the analysis, summarises dataset composition and the patient cohort,
reports quality-control metrics, verifies the normalisation performed upstream
by the host laboratory, and validates the existing cell-type annotation with
canonical markers and UMAP projections.
 
Dissertation sections    : 2.2, 2.4, 2.5, 3.1, 3.2
Tables produced          : Table 1, Table 2 (lower part), Table 3, Table 6
Figures produced         : Figure 1A, Figure 1B, Figure 2, Figure 6
Supplementary Appendix   : A0, A1.1-A1.5, A1.7, A2
 
Prerequisites
-------------
This is the first script in the sequence; run it before any other. The scripts
in this repository operate on a single AnnData object held in memory and are
intended to be executed sequentially within one Python session. Later scripts
depend on columns added by earlier ones and must not re-load the object.
 
Before running, set INPUT_OBJECT below to the path of the processed AnnData
object and OUTPUT_DIR to a folder for the tables and figures. Both are left as
placeholders here: the dataset is not public and no path to it is recorded in
this repository.
 
The object is expected to have log1p-transformed, library-size-normalised
values in .X (target sum 10,000), raw counts in the counts_RNA layer, a Harmony
embedding in .obsm, and an existing cell-type annotation. Normalisation,
integration and that annotation were performed upstream by the host laboratory
and are not re-run here; Section A1.5 below only verifies them.
 
Note on Table 2: the code here produces the lower part of that table (patients
per stage group, type 2 diabetes counts and percentages, and the F0 to F3
counts). The age, sex, BMI and serum ALT rows were compiled from patient-level
clinical data supplied separately by the host laboratory and have no
corresponding code.
"""
 
 
# ============================================================================
# A0. Software environment
# ============================================================================
 
import scanpy, anndata, pandas, numpy, scipy, statsmodels
import skimage, matplotlib, seaborn
 
for m in [scanpy, anndata, pandas, numpy, scipy, statsmodels,
          skimage, matplotlib, seaborn]:
    print(f"{m.__name__}=={m.__version__}")
 
# scanpy==1.11.5      anndata==0.12.16    pandas==2.3.3
# numpy==2.4.6        scipy==1.17.1       statsmodels==0.14.6
# scikit-image==0.26.0                    matplotlib==3.10.9
# seaborn==0.13.2    (used for the faceted boxplots in A8.2)
 
 
# ============================================================================
# A1.1 Loading and setup
# ============================================================================
 
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
 
# Local paths removed; supply the processed object and an output folder here.
INPUT_OBJECT = Path("<processed AnnData object provided by the host laboratory>")
OUTPUT_DIR   = Path("<output directory>")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
 
adata = sc.read_h5ad(INPUT_OBJECT)
print(adata)
print("cells x genes:", adata.shape)
 
# Metadata columns used throughout
PATIENT_COL  = "Patient_ID"
SAMPLE_COL   = "Sample"
TISSUE_COL   = "Tissue"
STAGE_COL    = "StageGroup"
STAGESEP_COL = "StageSep"
CELLTYPE_COL = "cell_type_with_cluster"
CLUSTER_COL  = "leiden_0.5"
 
# Drop categories with no cells, so that empty groups do not appear in tables
for col in [CELLTYPE_COL, TISSUE_COL, STAGE_COL]:
    if hasattr(adata.obs[col], "cat"):
        adata.obs[col] = adata.obs[col].cat.remove_unused_categories()
 
 
# ============================================================================
# A1.2 Dataset composition (Table 1)
# ============================================================================
 
obs = adata.obs.copy()
obs[TISSUE_COL] = obs[TISSUE_COL].astype(str)
obs[STAGE_COL]  = obs[STAGE_COL].astype(str)
 
def clean_celltype(x):
    """Standardise the three annotated population labels."""
    x_low = str(x).lower()
    if "cd14" in x_low:
        return "CD14 monocytes"
    elif "cd16" in x_low or "fcgr3a" in x_low:
        return "CD16 monocytes"
    elif "macro" in x_low:
        return "Macrophages"
    else:
        return str(x)
 
obs["celltype_clean"] = obs[CELLTYPE_COL].apply(clean_celltype)
 
def count_percent(series, category, denominator):
    counts = series.value_counts(dropna=False)
    df = counts.rename_axis("Group").reset_index(name="n")
    df.insert(0, "Category", category)
    df["% of total cells"] = (df["n"] / denominator * 100).round(1)
    return df
 
rows = [pd.DataFrame([
    {"Category": "Dataset", "Group": "Total cells",
     "n": adata.n_obs, "% of total cells": 100.0},
    {"Category": "Dataset", "Group": "Total genes",
     "n": adata.n_vars, "% of total cells": np.nan},
    {"Category": "Dataset", "Group": "Unique Patient_ID",
     "n": obs[PATIENT_COL].nunique(), "% of total cells": np.nan},
])]
 
rows.append(count_percent(obs[TISSUE_COL],        "Tissue",          adata.n_obs))
rows.append(count_percent(obs[STAGE_COL],         "Stage group",     adata.n_obs))
rows.append(count_percent(obs["celltype_clean"],  "Broad cell type", adata.n_obs))
 
table1 = pd.concat(rows, ignore_index=True)
category_order = ["Dataset", "Tissue", "Stage group", "Broad cell type"]
table1["Category"] = pd.Categorical(table1["Category"],
                                    categories=category_order, ordered=True)
table1 = table1.sort_values(["Category", "Group"]).reset_index(drop=True)
table1.to_csv(OUTPUT_DIR / "Table1_dataset_composition.csv", index=False)
print(table1.to_string(index=False))
 
 
# ============================================================================
# A1.3 Patient-level clinical summary (Table 2, lower part)
# ============================================================================
 
pat_clinical = (adata.obs
                .groupby(PATIENT_COL, observed=True)
                .agg(StageGroup=(STAGE_COL, "first"),
                     Stage=("Stage", "first"),
                     StageSep=(STAGESEP_COL, "first"),
                     Diabetic=("Diabetic", "first"),
                     n_cells=(TISSUE_COL, "size"))
                .reset_index())
 
print("Patients per stage group:")
print(pat_clinical["StageGroup"].value_counts().to_string())
print("\nDiabetes status x stage group:")
print(pd.crosstab(pat_clinical["StageGroup"], pat_clinical["Diabetic"]).to_string())
print("\nDetailed fibrosis stage x stage group:")
print(pd.crosstab(pat_clinical["StageGroup"], pat_clinical["Stage"]).to_string())
 
 
# ============================================================================
# A1.4 Quality-control metrics (Table 3)
# ============================================================================
 
qc_cols  = ["nCount_RNA", "nFeature_RNA", "percent.mt", "percent.hb"]
key_cols = qc_cols + [TISSUE_COL, PATIENT_COL, STAGE_COL, CELLTYPE_COL]
 
# Overall QC summary -> Table 3
adata.obs[qc_cols].describe().T.to_csv(OUTPUT_DIR / "qc_overall_summary.csv")
 
# Missing-value check and QC stratified by tissue and by cell type
adata.obs[key_cols].isna().sum().to_csv(OUTPUT_DIR / "missing_values.csv")
adata.obs.groupby(TISSUE_COL)[qc_cols].describe().to_csv(
    OUTPUT_DIR / "qc_by_tissue.csv")
adata.obs.groupby(CELLTYPE_COL)[qc_cols].describe().to_csv(
    OUTPUT_DIR / "qc_by_cell_type.csv")
 
# Distribution plots
sc.pl.violin(adata, qc_cols, groupby=TISSUE_COL, rotation=45, show=False)
plt.savefig(OUTPUT_DIR / "qc_violin_by_tissue.png", dpi=300, bbox_inches="tight")
plt.close()
 
sc.pl.violin(adata, qc_cols, groupby=CELLTYPE_COL, rotation=45, show=False)
plt.savefig(OUTPUT_DIR / "qc_violin_by_cell_type.png", dpi=300, bbox_inches="tight")
plt.close()
 
sc.pl.scatter(adata, x="nCount_RNA", y="nFeature_RNA", show=False)
plt.savefig(OUTPUT_DIR / "scatter_counts_features.png", dpi=300, bbox_inches="tight")
plt.close()
 
sc.pl.scatter(adata, x="nCount_RNA", y="percent.mt", show=False)
plt.savefig(OUTPUT_DIR / "scatter_counts_mt.png", dpi=300, bbox_inches="tight")
plt.close()
 
# Haemoglobin-high cells, checked across compartments (Section 2.4)
hb = pd.to_numeric(adata.obs["percent.hb"], errors="coerce")
for t in [1, 5, 10, 20]:
    print(f"percent.hb > {t}%: {(hb > t).sum()} cells",
          "| by tissue:", adata.obs.loc[hb > t, TISSUE_COL].value_counts().to_dict())
print("total cells:", adata.n_obs)
 
 
# ============================================================================
# A1.5 Normalisation and provenance checks
# ============================================================================
 
import scipy.sparse as sp
 
def densify(M, n=None):
    if n is not None:
        M = M[:n]
    return M.toarray() if sp.issparse(M) else np.asarray(M)
 
# --- Check the identity of .X (verification only; not a normalisation step) ---
Xs = densify(adata.X, 1000)
print("min:", Xs.min(), "max:", Xs.max())
print("has negatives (=> scaled/z-scored):", bool((Xs < 0).any()))
print("all ~integer (=> raw counts):", bool(np.allclose(Xs, np.round(Xs))))
print("expm1(.X) row sums (first 8):", np.round(np.expm1(Xs).sum(1)[:8], 1))
#   -> if approximately constant (e.g. ~1e4), .X is already log1p-normalised
for ln in adata.layers:
    L = densify(adata.layers[ln], 1000)
    print(f"layer '{ln}': min {L.min():.3f} max {L.max():.3f} "
          f"integer={np.allclose(L, np.round(L))}")
 
# --- Confirm the candidate marker genes are present in the object -------------
# NOTE: these are the candidate genes checked for presence at this stage, not
# the final signature definitions. The complement list here additionally
# includes the tissue-macrophage markers MRC1, MARCO and TIMD4, which were
# used for annotation checking only. The four signatures as scored consist of
# 6, 4, 5 and 5 genes respectively; see Section A5 and Table 5.
sigs = {"CD14":           ["CD14", "FCN1", "S100A8", "S100A9", "VCAN", "CCR2"],
        "CD16":           ["FCGR3A", "MS4A7", "LILRB1", "CX3CR1"],
        "Mac_complement": ["C1QA", "C1QB", "C1QC", "APOE", "APOC1",
                           "MRC1", "MARCO", "TIMD4"],
        "Mac_DAM":        ["TREM2", "CD9", "SPP1", "GPNMB", "LGALS3"]}
vn = set(adata.var_names)
for s, g in sigs.items():
    print(s, "missing:", [x for x in g if x not in vn] or "none")
 
# --- Patient versus sample, and patient x tissue pairing ----------------------
print("n patient:", adata.obs[PATIENT_COL].nunique(),
      "| n sample:", adata.obs[SAMPLE_COL].nunique())
print(adata.obs.drop_duplicates(SAMPLE_COL)
      .groupby(PATIENT_COL)[SAMPLE_COL].nunique())
 
ct = pd.crosstab(adata.obs[PATIENT_COL], adata.obs[TISSUE_COL])
print(ct)
print("tissues per patient:\n", ct.gt(0).sum(1).value_counts())
 
# --- Upstream QC, hashing and doublet traces ----------------------------------
print("HTO_classification.global:\n",
      adata.obs["HTO_classification.global"].value_counts(dropna=False))
print("\nLIVER HTO global:\n",
      adata.obs.loc[adata.obs[TISSUE_COL] == "LIVER",
                    "HTO_classification.global"].value_counts(dropna=False))
mt = pd.to_numeric(adata.obs["percent.mt"], errors="coerce")
hb = pd.to_numeric(adata.obs["percent.hb"], errors="coerce")
print("percent.mt max:", mt.max(), "| percent.mt median:", mt.median())
print("percent.hb max:", hb.max())
print("nCount_RNA median:", adata.obs["nCount_RNA"].median(),
      "| nFeature_RNA median:", adata.obs["nFeature_RNA"].median())
 
# --- Cells with missing mitochondrial / haemoglobin values --------------------
for col in ["percent.mt", "percent.hb"]:
    s = pd.to_numeric(adata.obs[col], errors="coerce")
    na = s.isna()
    print(f"{col}: {int(na.sum())} missing | by tissue:",
          adata.obs.loc[na, TISSUE_COL].value_counts().to_dict())
 
# --- Embeddings and stored parameters computed upstream -----------------------
for k in adata.obsm:
    print(f"obsm['{k}'] shape {adata.obsm[k].shape} (n_obs={adata.n_obs})")
print("uns['neighbors_subcluster']:", adata.uns.get("neighbors_subcluster"))
print("rank_genes_groups params:",
      adata.uns.get("rank_genes_groups", {}).get("params"))
 
 
# ============================================================================
# A1.7 Cell-type distribution across compartments (Table 6)
# ============================================================================
 
tab6 = pd.crosstab(adata.obs[TISSUE_COL], adata.obs[CELLTYPE_COL])
print(tab6)
 
# Percentage within each compartment
print((pd.crosstab(adata.obs[TISSUE_COL], adata.obs[CELLTYPE_COL],
                   normalize="index") * 100).round(1))
 
 
# ============================================================================
# A2. Annotation validation and dimensionality reduction (Figures 1A, 1B, 2, 6)
# ============================================================================
 
# --- Cell-type counts and cross-tabulations -----------------------------------
adata.obs[CELLTYPE_COL].value_counts().to_csv(OUTPUT_DIR / "cell_type_counts.csv")
pd.crosstab(adata.obs[TISSUE_COL], adata.obs[CELLTYPE_COL]).to_csv(
    OUTPUT_DIR / "cell_type_by_tissue.csv")
pd.crosstab(adata.obs[STAGE_COL], adata.obs[CELLTYPE_COL]).to_csv(
    OUTPUT_DIR / "cell_type_by_stagegroup.csv")
 
# --- UMAP by annotated cell type (Figure 1A) ----------------------------------
sc.pl.umap(adata, color=CELLTYPE_COL, show=False)
plt.savefig(OUTPUT_DIR / "umap_by_cell_type.png", dpi=300, bbox_inches="tight")
plt.close()
 
# --- UMAP by tissue compartment (Figure 2) ------------------------------------
sc.pl.umap(adata, color=TISSUE_COL, show=False)
plt.savefig(OUTPUT_DIR / "umap_by_tissue.png", dpi=300, bbox_inches="tight")
plt.close()
 
# --- UMAP by fibrosis-stage group (Figure 6) ----------------------------------
sc.pl.umap(adata, color=STAGE_COL, show=False)
plt.savefig(OUTPUT_DIR / "umap_by_stagegroup.png", dpi=300, bbox_inches="tight")
plt.close()
 
# --- Canonical marker panel (27 genes; Section 2.5) ---------------------------
marker_genes = [
    # general myeloid / monocyte markers
    "LYZ", "LST1", "TYROBP", "AIF1",
 
    # CD14 monocyte markers
    "CD14", "FCN1", "S100A8", "S100A9", "VCAN", "CCR2",
 
    # CD16 monocyte markers
    "FCGR3A", "MS4A7", "LILRB1", "CX3CR1",
 
    # macrophage / Kupffer-like / tissue macrophage markers
    "C1QA", "C1QB", "C1QC", "APOE", "APOC1", "MRC1", "MARCO", "TIMD4",
 
    # disease-associated / scar-associated macrophage markers
    "TREM2", "CD9", "SPP1", "GPNMB", "LGALS3"
]
 
marker_genes = list(dict.fromkeys(marker_genes))   # de-duplicate, keep order
present_markers = [g for g in marker_genes if g in adata.var_names]
missing_markers = [g for g in marker_genes if g not in adata.var_names]
 
pd.DataFrame({"present_markers": pd.Series(present_markers)}).to_csv(
    OUTPUT_DIR / "present_marker_genes.csv", index=False)
pd.DataFrame({"missing_markers": pd.Series(missing_markers)}).to_csv(
    OUTPUT_DIR / "missing_marker_genes.csv", index=False)
 
# --- Marker UMAP feature plots, in batches to keep panels legible -------------
batch_size = 6
for i in range(0, len(present_markers), batch_size):
    batch = present_markers[i:i + batch_size]
    sc.pl.umap(adata, color=batch, use_raw=False, show=False)
    plt.savefig(OUTPUT_DIR / f"umap_marker_genes_batch_{i // batch_size + 1}.png",
                dpi=300, bbox_inches="tight")
    plt.close()
 
# --- Dot plot by annotated cell type (Figure 1B) ------------------------------
dp = sc.pl.dotplot(adata, var_names=present_markers, groupby=CELLTYPE_COL,
                   standard_scale="var", use_raw=False, return_fig=True)
dp.savefig(OUTPUT_DIR / "dotplot_markers_by_cell_type.png")
 
# --- Matrix plot by annotated cell type ---------------------------------------
mp = sc.pl.matrixplot(adata, var_names=present_markers, groupby=CELLTYPE_COL,
                      standard_scale="var", use_raw=False, return_fig=True)
mp.savefig(OUTPUT_DIR / "matrixplot_markers_by_cell_type.png")
