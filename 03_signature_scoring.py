"""
03_signature_scoring.py

Score four pre-defined transcriptional programmes in every cell using
scanpy.tl.score_genes.

The four signatures are canonical marker panels for CD14 monocytes, CD16
monocytes, complement-associated macrophages and disease-associated
macrophages. score_genes computes, for each cell, the mean expression of the
target gene set minus the mean expression of an expression-matched set of
control genes, so the scores are relative (negative values are expected).

Input : mono_macro_with_subclusters.h5ad  (from 02_subclustering.py)
Output: mean signature scores by subcluster, cell type, tissue and stage group,
        and an updated .h5ad containing the four score columns
"""

import os
import scanpy as sc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_PATH = "outputs/subclustering/mono_macro_with_subclusters.h5ad"
OUTPUT_DIR = "outputs/signatures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

adata = sc.read_h5ad(DATA_PATH)

# ---------------------------------------------------------------------------
# Gene sets for the four signatures
# ---------------------------------------------------------------------------
signatures = {
    "CD14_monocyte_score": ["CD14", "FCN1", "S100A8", "S100A9", "VCAN", "CCR2"],
    "CD16_monocyte_score": ["FCGR3A", "MS4A7", "LILRB1", "CX3CR1"],
    "Macrophage_complement_score": ["C1QA", "C1QB", "C1QC", "APOE", "APOC1"],
    "Disease_associated_macrophage_score": ["TREM2", "CD9", "SPP1", "GPNMB", "LGALS3"],
}

# ---------------------------------------------------------------------------
# Score each signature
# ---------------------------------------------------------------------------
# use_raw=False -> score on the log-normalised .X matrix.
# Default control parameters are used (ctrl_size=50, n_bins=25).
# A score is only computed if at least two genes of the set are present.
score_names = []
for score_name, genes in signatures.items():
    genes_present = [g for g in genes if g in adata.var_names]
    if len(genes_present) >= 2:
        sc.tl.score_genes(
            adata,
            gene_list=genes_present,
            score_name=score_name,
            use_raw=False,
        )
        score_names.append(score_name)

# ---------------------------------------------------------------------------
# Mean signature scores summarised by different groupings
# ---------------------------------------------------------------------------
final_cluster_key = "leiden_0.5"
for group, name in [(final_cluster_key, "subcluster"),
                    ("cell_type_with_cluster", "cell_type"),
                    ("Tissue", "tissue"),
                    ("StageGroup", "stagegroup")]:
    adata.obs.groupby(group, observed=True)[score_names].mean().to_csv(
        f"{OUTPUT_DIR}/signature_scores_by_{name}.csv"
    )

# ---------------------------------------------------------------------------
# Save the object with the four score columns for the statistical analysis
# ---------------------------------------------------------------------------
adata.write_h5ad(f"{OUTPUT_DIR}/mono_macro_with_scores.h5ad")

print(f"Signature scores written to: {OUTPUT_DIR}")
print("Scores computed:", score_names)
