"""
02_subclustering.py

Graph-based subclustering of the monocyte/macrophage populations, followed by
marker-based interpretation of the resulting subclusters.

Clustering is performed on the pre-computed Harmony embedding (batch key =
sequencing run). Several Leiden resolutions are tested; resolution 0.5 is used
for the main analysis, giving 13 subclusters. Each subcluster is then
characterised using canonical markers and cluster-specific differential
expression (Wilcoxon rank-sum test).

Input : mono_macro.h5ad
Output: subcluster composition tables, marker plots, ranked marker genes,
        and an updated .h5ad containing the subcluster labels
"""

import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_PATH = "data/mono_macro.h5ad"
OUTPUT_DIR = "outputs/subclustering"
os.makedirs(OUTPUT_DIR, exist_ok=True)

adata = sc.read_h5ad(DATA_PATH)

# Drop empty categories so that 0-count groups do not appear in tables/plots
for col in ["cell_type_with_cluster", "Tissue", "StageGroup"]:
    if col in adata.obs.columns and hasattr(adata.obs[col], "cat"):
        adata.obs[col] = adata.obs[col].cat.remove_unused_categories()

# ---------------------------------------------------------------------------
# 1. Neighbour graph on the pre-computed Harmony embedding
# ---------------------------------------------------------------------------
# Harmony was computed upstream using the sample identifier (= sequencing run)
# as the batch key. We reuse these coordinates for consistency with the host
# laboratory's integration scheme.
use_rep = "X_harmony" if "X_harmony" in adata.obsm.keys() else "X_pca"

sc.pp.neighbors(
    adata,
    use_rep=use_rep,
    n_neighbors=15,
    key_added="neighbors_subcluster",
)

# ---------------------------------------------------------------------------
# 2. Leiden clustering across a range of resolutions
# ---------------------------------------------------------------------------
# Resolution controls cluster granularity. We test several values and use 0.5
# for the main analysis.
resolutions = [0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
for res in resolutions:
    sc.tl.leiden(
        adata,
        resolution=res,
        key_added=f"leiden_{res}",
        neighbors_key="neighbors_subcluster",
    )

final_cluster_key = "leiden_0.5"

# ---------------------------------------------------------------------------
# 3. Subcluster composition (by cell type, tissue and fibrosis-stage group)
# ---------------------------------------------------------------------------
# Cell counts per subcluster
adata.obs[final_cluster_key].value_counts().sort_index().to_csv(
    f"{OUTPUT_DIR}/subcluster_counts.csv"
)

# Row-normalised composition tables (percentage within each subcluster)
for var, name in [("cell_type_with_cluster", "cell_type"),
                  ("Tissue", "tissue"),
                  ("StageGroup", "stagegroup")]:
    (pd.crosstab(adata.obs[final_cluster_key], adata.obs[var], normalize="index") * 100
     ).to_csv(f"{OUTPUT_DIR}/subcluster_by_{name}_percent.csv")

# ---------------------------------------------------------------------------
# 4. Marker genes used to interpret the subclusters
# ---------------------------------------------------------------------------
marker_genes = [
    "LYZ", "LST1", "TYROBP", "AIF1",                              # general myeloid
    "CD14", "FCN1", "S100A8", "S100A9", "VCAN", "CCR2",           # CD14 monocyte
    "FCGR3A", "MS4A7", "LILRB1", "CX3CR1",                        # CD16 monocyte
    "C1QA", "C1QB", "C1QC", "APOE", "APOC1", "MRC1", "MARCO", "TIMD4",  # tissue macrophage
    "TREM2", "CD9", "SPP1", "GPNMB", "LGALS3",                    # disease-associated
]
present = [g for g in dict.fromkeys(marker_genes) if g in adata.var_names]

# Dot plot and matrix plot of markers across subclusters
sc.pl.dotplot(adata, var_names=present, groupby=final_cluster_key,
              standard_scale="var", use_raw=False, return_fig=True
              ).savefig(f"{OUTPUT_DIR}/dotplot_markers.png")
sc.pl.matrixplot(adata, var_names=present, groupby=final_cluster_key,
                 standard_scale="var", use_raw=False, return_fig=True
                 ).savefig(f"{OUTPUT_DIR}/matrixplot_markers.png")

# ---------------------------------------------------------------------------
# 5. Cluster-specific differential expression (Wilcoxon rank-sum)
# ---------------------------------------------------------------------------
# This confirms that each subcluster is enriched for its expected markers.
sc.tl.rank_genes_groups(
    adata,
    groupby=final_cluster_key,
    method="wilcoxon",
    use_raw=False,
    n_genes=50,
)
ranked = sc.get.rank_genes_groups_df(adata, group=None)
ranked.to_csv(f"{OUTPUT_DIR}/ranked_marker_genes.csv", index=False)

# ---------------------------------------------------------------------------
# 6. Save the object with subcluster labels for downstream steps
# ---------------------------------------------------------------------------
adata.write_h5ad(f"{OUTPUT_DIR}/mono_macro_with_subclusters.h5ad")

print(f"Subclustering outputs written to: {OUTPUT_DIR}")
