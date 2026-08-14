"""
02_subclustering.py
 
Builds a 15-nearest-neighbour graph on the pre-computed Harmony embedding, runs
Leiden clustering at six resolutions (resolution 0.5 is used for all downstream
summaries), assembles the subcluster composition table, and ranks marker genes
at both the broad cell-type level and the subcluster level.
 
Dissertation sections    : 2.6, 3.1, 3.3, 3.5
Tables produced          : Table 7
Figures produced         : Figure 3A, Figure 3B
Supplementary Appendix   : A3.1-A3.3, A4.1-A4.3
 
Prerequisites
-------------
Depends on 01. Continues from the previous script in the same Python session;
do not re-load the object. Only the imports and the metadata column names are
repeated below.
 
Note on rank_genes_groups: sc.tl.rank_genes_groups overwrites
adata.uns["rank_genes_groups"] on each call. A4.1 (broad cell types) is run
before A4.2 (subclusters) so that the subcluster result is the one left in .uns
for the identity checks in A4.3. Keep this order if re-running.
"""
 
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
 
# --- Metadata column names (repeated from 01 for readability) ----------------
PATIENT_COL  = "Patient_ID"
SAMPLE_COL   = "Sample"
TISSUE_COL   = "Tissue"
STAGE_COL    = "StageGroup"
STAGESEP_COL = "StageSep"
CELLTYPE_COL = "cell_type_with_cluster"
CLUSTER_COL  = "leiden_0.5"
 
 
# ============================================================================
# A3.1 Neighbour graph and Leiden clustering
# ============================================================================
 
# Select the low-dimensional representation for subclustering
if "X_harmony" in adata.obsm.keys():
    use_rep = "X_harmony"
elif "X_pca" in adata.obsm.keys():
    use_rep = "X_pca"
else:
    raise ValueError("Neither X_harmony nor X_pca found; cannot build the "
                     "neighbour graph for Leiden subclustering.")
 
# 1. Neighbour graph
sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15,
                key_added="neighbors_subcluster")
 
# 2. Leiden clustering across several resolutions, to compare granularity
resolutions = [0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}",
                 neighbors_key="neighbors_subcluster")
 
cluster_summary = {f"leiden_{res}": adata.obs[f"leiden_{res}"]
                   .value_counts().sort_index() for res in resolutions}
pd.DataFrame(cluster_summary).fillna(0).astype(int).to_csv(
    OUTPUT_DIR / "leiden_resolution_cluster_counts.csv")
 
# 3. Resolution 0.5 is the key used for all downstream summaries
final_cluster_key = CLUSTER_COL          # "leiden_0.5"
adata.obs[final_cluster_key].value_counts().sort_index().to_csv(
    OUTPUT_DIR / "final_subcluster_counts.csv")
 
# 4. UMAP coloured by each resolution, for comparison
for res in resolutions:
    sc.pl.umap(adata, color=f"leiden_{res}", legend_loc="right margin", show=False)
    plt.savefig(OUTPUT_DIR / f"umap_leiden_{res}.png", dpi=300, bbox_inches="tight")
    plt.close()
 
 
# ============================================================================
# A3.2 Subcluster cross-tabulation and marker plots (Figure 3A)
# ============================================================================
 
# Subcluster x broad cell type
pd.crosstab(adata.obs[final_cluster_key], adata.obs[CELLTYPE_COL]).to_csv(
    OUTPUT_DIR / "subcluster_by_cell_type_counts.csv")
(pd.crosstab(adata.obs[final_cluster_key], adata.obs[CELLTYPE_COL],
             normalize="index") * 100).to_csv(
    OUTPUT_DIR / "subcluster_by_cell_type_percent.csv")
 
# Subcluster x tissue
pd.crosstab(adata.obs[final_cluster_key], adata.obs[TISSUE_COL]).to_csv(
    OUTPUT_DIR / "subcluster_by_tissue_counts.csv")
(pd.crosstab(adata.obs[final_cluster_key], adata.obs[TISSUE_COL],
             normalize="index") * 100).to_csv(
    OUTPUT_DIR / "subcluster_by_tissue_percent.csv")
 
# Subcluster x fibrosis-stage group
pd.crosstab(adata.obs[final_cluster_key], adata.obs[STAGE_COL]).to_csv(
    OUTPUT_DIR / "subcluster_by_stagegroup_counts.csv")
(pd.crosstab(adata.obs[final_cluster_key], adata.obs[STAGE_COL],
             normalize="index") * 100).to_csv(
    OUTPUT_DIR / "subcluster_by_stagegroup_percent.csv")
 
# UMAP by final subcluster (Figure 3A)
sc.pl.umap(adata, color=final_cluster_key, legend_loc="right margin", show=False)
plt.savefig(OUTPUT_DIR / "umap_leiden_0.5.png", dpi=300, bbox_inches="tight")
plt.close()
 
# Canonical marker panel by subcluster
dp = sc.pl.dotplot(adata, var_names=present_markers, groupby=final_cluster_key,
                   standard_scale="var", use_raw=False, return_fig=True)
dp.savefig(OUTPUT_DIR / "dotplot_markers_by_final_subcluster.png")
 
mp = sc.pl.matrixplot(adata, var_names=present_markers, groupby=final_cluster_key,
                      standard_scale="var", use_raw=False, return_fig=True)
mp.savefig(OUTPUT_DIR / "matrixplot_markers_by_final_subcluster.png")
 
 
# ============================================================================
# A3.3 Subcluster composition table (Table 7)
# ============================================================================
 
obs = adata.obs.copy()
 
def clean_cluster_label(x):
    """Format cluster labels as C0, C1, ... rather than 0, 1, ..."""
    s = str(x)
    if s.startswith("C"):
        return s
    try:
        return "C" + str(int(float(s)))
    except ValueError:
        return s
 
def cluster_sort_key(x):
    s = str(x).replace("C", "")
    try:
        return int(float(s))
    except ValueError:
        return 9999
 
obs["cluster_label"] = obs[final_cluster_key].astype(str).apply(clean_cluster_label)
 
cluster_counts = (obs["cluster_label"].value_counts()
                  .rename_axis("Cluster").reset_index(name="n")
                  .sort_values("Cluster", key=lambda x: x.map(cluster_sort_key)))
 
# Row-wise percentages within each subcluster
tissue_pct = (pd.crosstab(obs["cluster_label"], obs[TISSUE_COL].astype(str),
                          normalize="index") * 100).round(1)
tissue_pct = tissue_pct.reset_index().rename(columns={"cluster_label": "Cluster"})
 
stage_pct = (pd.crosstab(obs["cluster_label"], obs[STAGE_COL].astype(str),
                         normalize="index") * 100).round(1)
stage_pct = stage_pct.reset_index().rename(columns={"cluster_label": "Cluster"})
 
table7 = (cluster_counts
          .merge(tissue_pct, on="Cluster", how="left")
          .merge(stage_pct,  on="Cluster", how="left"))
 
identity_map = {
    "C0":  "Classical CD14 monocytes",
    "C1":  "Lymphoid, non-myeloid (excluded from interpretation)",
    "C2":  "Lipid-/disease-associated macrophages",
    "C3":  "Inflammatory monocytes",
    "C4":  "Complement-high resident macrophages",
    "C5":  "Non-classical CD16 monocytes",
    "C6":  "MRC1+ macrophages",
    "C7":  "CD16 monocytes",
    "C8":  "Transitional monocytes",
    "C9":  "Erythroid, non-myeloid (excluded from interpretation)",
    "C10": "Kupffer-like resident macrophages",
    "C11": "Platelet, non-myeloid (excluded from interpretation)",
    "C12": "Lipid-associated macrophages (FABP4+)",
}
 
table7["Identity"] = table7["Cluster"].astype(str).map(identity_map)
table7["Identity"] = table7["Identity"].astype("object").fillna("")
table7 = table7.sort_values("Cluster", key=lambda x: x.map(cluster_sort_key))
table7 = table7.reset_index(drop=True)
 
table7.to_csv(OUTPUT_DIR / "Table7_subcluster_composition.csv", index=False)
print(table7.to_string(index=False))
 
 
# ============================================================================
# A3.3 Stage composition of C2 and C12 against the dataset baseline (Section 3.5)
# ============================================================================
 
# Dataset-wide F0_1 fraction, used as the baseline
total_f01 = (adata.obs[STAGE_COL] == "F0_1").mean()
print(f"Dataset-wide F0_1 fraction (baseline): {total_f01:.1%}\n")
 
# F0_1 fraction within C2 and C12
for c in ["2", "12"]:
    mask = adata.obs[CLUSTER_COL] == c
    c_f01 = (adata.obs.loc[mask, STAGE_COL] == "F0_1").mean()
    n = mask.sum()
    print(f"C{c}: F0_1 fraction {c_f01:.1%} (n={n}, baseline {total_f01:.1%})")
 
 
# ============================================================================
# A4.1 Broad cell-type markers (Section 3.1)
# ============================================================================
 
sc.tl.rank_genes_groups(adata, groupby=CELLTYPE_COL,
                        method="wilcoxon", use_raw=False, n_genes=50)
 
ranked_by_celltype = sc.get.rank_genes_groups_df(adata, group=None)
ranked_by_celltype.to_csv(OUTPUT_DIR / "ranked_marker_genes_by_CELLTYPE.csv",
                          index=False)
print(ranked_by_celltype.head(20))
 
 
# ============================================================================
# A4.2 Subcluster markers (Section 3.3, Figure 3B)
# ============================================================================
 
sc.tl.rank_genes_groups(adata, groupby=final_cluster_key,
                        method="wilcoxon", use_raw=False, n_genes=50)
 
ranked_markers_df = sc.get.rank_genes_groups_df(adata, group=None)
ranked_markers_df.to_csv(
    OUTPUT_DIR / "ranked_marker_genes_by_final_subcluster.csv", index=False)
 
# Top 10 markers per subcluster, for quick inspection
top10_markers = ranked_markers_df.groupby("group").head(10).reset_index(drop=True)
top10_markers.to_csv(
    OUTPUT_DIR / "top10_marker_genes_by_final_subcluster.csv", index=False)
 
# Ranked-marker dot plot, top 5 genes per subcluster (Figure 3B)
rg = sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby=final_cluster_key,
                                     standard_scale="var", show=False, return_fig=True)
rg.savefig(OUTPUT_DIR / "rank_genes_groups_dotplot_top5.png")
 
 
# ============================================================================
# A4.3 Subcluster identity checks (Section 3.3)
# ============================================================================
 
# Top 10 ranked genes for a single subcluster (example: C8)
result = adata.uns["rank_genes_groups"]
c8_markers = [result["names"]["8"][i] for i in range(10)]
print("C8 top 10 markers:", c8_markers)
 
# Composition of C8 by broad cell type
print(adata.obs[adata.obs[final_cluster_key] == "8"][CELLTYPE_COL].value_counts())
 
# Mapping of every Leiden subcluster onto the three annotated populations
ct = pd.crosstab(adata.obs[final_cluster_key], adata.obs[CELLTYPE_COL])
print("leiden_0.5 vs cell_type_with_cluster:")
print(ct)
 
