"""
03_signature_scoring.py
 
Computes the four gene-signature scores with sc.tl.score_genes on the
log-normalised matrix, using the Scanpy default control parameters
(ctrl_size=50, n_bins=25), and summarises the mean score by subcluster, broad
cell type, tissue compartment and fibrosis-stage group.
 
Dissertation sections    : 2.6
Tables produced          : Table 5 (gene composition and scoring settings)
Supplementary Appendix   : A5
 
The Basis column of Table 5, which records the published evidence for each
signature, was compiled manually from the literature and has no corresponding
code.
 
Because scores are computed relative to expression-matched control genes,
negative values are expected and do not indicate absent expression.
 
Prerequisites
-------------
Depends on 02 (the subcluster labels are used for the per-subcluster summary).
Continues from the previous script in the same Python session; do not re-load
the object. Only the imports and the metadata column names are repeated below.
"""
 
import scanpy as sc
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
# A5. Gene-signature scoring (score_genes)
# ============================================================================
 
signatures = {
    "CD14_monocyte_score":                 ["CD14", "FCN1", "S100A8",
                                            "S100A9", "VCAN", "CCR2"],
    "CD16_monocyte_score":                 ["FCGR3A", "MS4A7", "LILRB1", "CX3CR1"],
    "Macrophage_complement_score":         ["C1QA", "C1QB", "C1QC", "APOE", "APOC1"],
    "Disease_associated_macrophage_score": ["TREM2", "CD9", "SPP1",
                                            "GPNMB", "LGALS3"],
}
 
score_names = []
for score_name, genes in signatures.items():
    genes_present = [g for g in genes if g in adata.var_names]
    if len(genes_present) >= 2:
        sc.tl.score_genes(adata, gene_list=genes_present,
                          score_name=score_name, use_raw=False)
        score_names.append(score_name)
 
# Mean signature scores by grouping variable
adata.obs.groupby(final_cluster_key)[score_names].mean().to_csv(
    OUTPUT_DIR / "signature_scores_by_final_subcluster.csv")
adata.obs.groupby(CELLTYPE_COL)[score_names].mean().to_csv(
    OUTPUT_DIR / "signature_scores_by_cell_type.csv")
adata.obs.groupby(TISSUE_COL)[score_names].mean().to_csv(
    OUTPUT_DIR / "signature_scores_by_tissue.csv")
adata.obs.groupby(STAGE_COL)[score_names].mean().to_csv(
    OUTPUT_DIR / "signature_scores_by_stagegroup.csv")
 
# Signature scores projected onto the UMAP
sc.pl.umap(adata, color=score_names, show=False)
plt.savefig(OUTPUT_DIR / "umap_signature_scores.png", dpi=300, bbox_inches="tight")
plt.close()
 
