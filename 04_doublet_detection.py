"""
04_doublet_detection.py

Objective doublet assessment using Scrublet.

No doublet score was stored in the original object, so Scrublet is applied to
the raw counts within each sequencing run (doublets form within a run). The
number of principal components is adapted to the size of each run, and runs
with fewer than 30 cells are too small to be scored and are skipped.

Input : mono_macro_with_subclusters.h5ad
Output: per-subcluster doublet summary table
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths  (requires scikit-image: pip install scikit-image)
# ---------------------------------------------------------------------------
DATA_PATH = "outputs/subclustering/mono_macro_with_subclusters.h5ad"
OUTPUT_DIR = "outputs/doublets"
os.makedirs(OUTPUT_DIR, exist_ok=True)

adata = sc.read_h5ad(DATA_PATH)

# Scrublet needs raw counts
adb = adata.copy()
adb.X = adb.layers["counts_RNA"].copy()

# ---------------------------------------------------------------------------
# Run Scrublet separately within each sequencing run
# ---------------------------------------------------------------------------
scores = pd.Series(np.nan, index=adata.obs_names, dtype=float)
preds = pd.Series(np.nan, index=adata.obs_names, dtype=float)
skipped = []

for s in adb.obs["Sample"].unique():
    sub = adb[adb.obs["Sample"] == s].copy()
    n = sub.n_obs

    # Runs with fewer than 30 cells are too small for the algorithm
    if n < 30:
        skipped.append((s, n))
        continue

    # Adapt the number of principal components to the run size
    n_pc = max(min(30, n // 3, sub.n_vars - 1), 2)
    try:
        sc.pp.scrublet(sub, n_prin_comps=n_pc, random_state=0, verbose=False)
        scores[sub.obs_names] = sub.obs["doublet_score"].values
        preds[sub.obs_names] = sub.obs["predicted_doublet"].astype(float).values
    except Exception as e:
        skipped.append((s, n))
        print(f"skipped {s} (n={n}): {type(e).__name__}")

adata.obs["doublet_score"] = scores.values
adata.obs["predicted_doublet"] = preds.values

print("Skipped runs:", skipped)
print("Cells scored: %d / %d" % (scores.notna().sum(), adata.n_obs))
print("Overall predicted-doublet rate: %.2f%%"
      % (100 * np.nanmean(adata.obs["predicted_doublet"])))

# ---------------------------------------------------------------------------
# Summarise doublet metrics per subcluster
# ---------------------------------------------------------------------------
summary = (
    adata.obs.dropna(subset=["doublet_score"])
    .groupby("leiden_0.5", observed=True)
    .agg(
        n_cells=("doublet_score", "size"),
        median_score=("doublet_score", "median"),
        mean_score=("doublet_score", "mean"),
        pct_predicted_doublet=("predicted_doublet", lambda x: 100 * np.mean(x)),
    )
    .round(3)
)
print(summary.to_string())
summary.to_csv(f"{OUTPUT_DIR}/doublet_by_subcluster.csv")

print(f"\nDoublet outputs written to: {OUTPUT_DIR}")
