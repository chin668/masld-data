"""
05_statistical_analysis.py

Statistical comparison of the four signature scores, using tests appropriate to
the study design.

Because all four compartments come from the same 20 patients, compartment
comparisons are treated as repeated measures (Friedman test + post-hoc Wilcoxon
signed-rank). Fibrosis-stage groups comprise different patients, so they are
compared with independent-samples tests (Kruskal-Wallis across stage groups;
Mann-Whitney U within each compartment). Effect sizes and bootstrap 95%
confidence intervals are reported alongside the p values, and Benjamini-Hochberg
correction is applied within each family of tests.

Input : mono_macro_with_scores.h5ad  (from 03_signature_scoring.py)
Output: Friedman, post-hoc, stage and within-compartment result tables
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon, kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import TTestIndPower

# ---------------------------------------------------------------------------
# Paths and settings
# ---------------------------------------------------------------------------
DATA_PATH = "outputs/signatures/mono_macro_with_scores.h5ad"
OUTPUT_DIR = "outputs/statistics"
os.makedirs(OUTPUT_DIR, exist_ok=True)

adata = sc.read_h5ad(DATA_PATH)

sig_cols = ["CD14_monocyte_score", "CD16_monocyte_score",
            "Macrophage_complement_score", "Disease_associated_macrophage_score"]
tissues = ["PBMC", "LIVER", "VAT", "SAT"]


def boot_ci(a, b, n=5000, seed=0):
    """Bootstrap 95% CI for the difference in means (a minus b)."""
    rng = np.random.default_rng(seed)
    diffs = [rng.choice(a, len(a), True).mean() - rng.choice(b, len(b), True).mean()
             for _ in range(n)]
    return np.percentile(diffs, [2.5, 97.5])


# ===========================================================================
# 1. COMPARTMENT COMPARISON  (paired -> Friedman + post-hoc Wilcoxon)
# ===========================================================================
# Aggregate to one value per patient per tissue.
pt = (adata.obs.groupby(["Patient_ID", "Tissue"], observed=True)[sig_cols]
      .mean().reset_index())

# Friedman test across the four compartments, per signature
rows = []
for s in sig_cols:
    wide = pt.pivot(index="Patient_ID", columns="Tissue", values=s)[tissues].dropna()
    chi2, p = friedmanchisquare(*[wide[t].values for t in tissues])
    rows.append({"signature": s, "n_patients": wide.shape[0],
                 "chi2": round(chi2, 2), "p_raw": p})
friedman = pd.DataFrame(rows)
friedman["FDR"] = multipletests(friedman["p_raw"], method="fdr_bh")[1]
friedman.to_csv(f"{OUTPUT_DIR}/friedman_tissue.csv", index=False)

# Post-hoc pairwise Wilcoxon signed-rank tests (24 comparisons, BH-corrected)
post = []
for s in sig_cols:
    wide = pt.pivot(index="Patient_ID", columns="Tissue", values=s)[tissues].dropna()
    for i in range(len(tissues)):
        for j in range(i + 1, len(tissues)):
            a, b = tissues[i], tissues[j]
            stat, p = wilcoxon(wide[a], wide[b])
            post.append({"signature": s, "comparison": f"{a} vs {b}",
                         "median_diff": round(np.median(wide[a] - wide[b]), 3),
                         "p_raw": p})
post = pd.DataFrame(post)
post["FDR"] = multipletests(post["p_raw"], method="fdr_bh")[1]
post.to_csv(f"{OUTPUT_DIR}/posthoc_wilcoxon_tissue.csv", index=False)

# ===========================================================================
# 2. STAGE COMPARISON  (independent -> patient-level Kruskal-Wallis)
# ===========================================================================
# Aggregate to ONE value per patient (mean across all compartments) to avoid
# pseudoreplication.
pat = (adata.obs.groupby("Patient_ID", observed=True)
       .agg(StageGroup=("StageGroup", "first"),
            **{c: (c, "mean") for c in sig_cols})
       .reset_index())

rows = []
for s in sig_cols:
    groups = [g[s].values for _, g in pat.groupby("StageGroup", observed=True)]
    h, p = kruskal(*groups)
    rows.append({"signature": s, "H": round(h, 2), "p_raw": p})
stage = pd.DataFrame(rows)
stage["FDR"] = multipletests(stage["p_raw"], method="fdr_bh")[1]
stage.to_csv(f"{OUTPUT_DIR}/stage_kruskal_patientlevel.csv", index=False)

# ===========================================================================
# 3. WITHIN-COMPARTMENT STAGE COMPARISON  (Mann-Whitney + effect size + CI)
# ===========================================================================
# F0_1 vs F2_3 inside each tissue, with rank-biserial effect size and a
# bootstrap 95% CI. The single healthy control is excluded.
pt_stage = (adata.obs.groupby(["Patient_ID", "Tissue"], observed=True)
            .agg(StageGroup=("StageGroup", "first"),
                 **{c: (c, "mean") for c in sig_cols})
            .reset_index())

rows = []
for t in tissues:
    sub = pt_stage[pt_stage["Tissue"] == t]
    for s in sig_cols:
        a = sub.loc[sub["StageGroup"] == "F0_1", s].values
        b = sub.loc[sub["StageGroup"] == "F2_3", s].values
        if len(a) < 2 or len(b) < 2:
            continue
        u, p = mannwhitneyu(a, b)
        rank_biserial = 1 - 2 * u / (len(a) * len(b))
        lo, hi = boot_ci(a, b)
        rows.append({"tissue": t, "signature": s,
                     "mean_F0_1": round(a.mean(), 3), "mean_F2_3": round(b.mean(), 3),
                     "diff": round(a.mean() - b.mean(), 3),
                     "CI_low": round(lo, 3), "CI_high": round(hi, 3),
                     "rank_biserial": round(rank_biserial, 2), "p_raw": round(p, 3)})
within = pd.DataFrame(rows)
within["FDR"] = multipletests(within["p_raw"], method="fdr_bh")[1]
within.to_csv(f"{OUTPUT_DIR}/within_compartment_stage.csv", index=False)

# ===========================================================================
# 4. POWER CHECK
# ===========================================================================
# Minimum detectable effect size (Cohen's d) at 80% power for 13 vs 6 patients.
mde = TTestIndPower().solve_power(nobs1=13, ratio=6 / 13, alpha=0.05, power=0.8)
print("Minimum detectable effect size (13 vs 6, 80%% power): %.2f SD" % mde)

print(f"\nStatistical outputs written to: {OUTPUT_DIR}")
