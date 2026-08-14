"""
05_statistical_analysis.py
 
All hypothesis testing reported in the dissertation:
 
  * cross-compartment comparison of the four signatures with Friedman tests on
    patient-level means, followed by pairwise Wilcoxon signed-rank post-hoc
    tests (Tables 8 upper half and 9, Figure 4);
  * the same comparison repeated within a single cell type under a minimum of
    10 cells per patient per compartment (Table 8, lower half);
  * fibrosis-stage comparisons with Mann-Whitney U tests, pooled across
    compartments and within each compartment (Tables 10 and 11, Figures 5
    and 7);
  * a sensitivity analysis excluding the three non-myeloid subclusters.
 
Aggregation levels differ within this file and should not be confused. The
Friedman and Mann-Whitney tests are performed on patient-level means. The
compartment means displayed in the upper half of Table 8 are cell-level means,
shown for description only. The boxplots in Figure 4 are at patient/tissue
level. Every significance conclusion reported in the dissertation comes from
the patient-level tests.
 
Benjamini-Hochberg correction is applied within each family of tests: across
the four signatures for the all-cell Friedman tests, across the 24 post-hoc
comparisons, across the three within-cell-type Friedman tests, separately
within each cell type for the within-cell-type post-hoc tests, across the four
signatures for the pooled stage comparison, and across the 16 within-
compartment stage comparisons.
 
Fibrosis-stage groups comprise different patients and are therefore treated as
independent; the single healthy control is excluded from all testing and shown
for reference only.
 
Dissertation sections    : 2.7, 3.3, 3.4, 3.5, 3.6
Tables produced          : Tables 8, 9, 10, 11
Figures produced         : Figures 4, 5, 7
Supplementary Appendix   : A6.1-A6.3, A7, A8.1, A8.2, A10
 
Prerequisites
-------------
Depends on 03 (signature scores) and, through it, on 02 (subcluster labels).
Continues from a previous script in the same Python session; do not re-load the
object. Only the imports and the metadata column names are repeated below.
"""
 
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
# A6.1 Friedman tests and post-hoc Wilcoxon signed-rank tests (Tables 8 and 9)
# ============================================================================
 
import pandas as pd, numpy as np
from scipy.stats import friedmanchisquare, wilcoxon
from statsmodels.stats.multitest import multipletests
 
sig_cols = ["CD14_monocyte_score", "CD16_monocyte_score",
            "Macrophage_complement_score", "Disease_associated_macrophage_score"]
tissues = ["PBMC", "LIVER", "VAT", "SAT"]
 
# Patient/tissue-level means: one value per patient per compartment
pt = (adata.obs.groupby([PATIENT_COL, TISSUE_COL], observed=True)[sig_cols]
      .mean().reset_index())
 
# --- Friedman: overall difference across the four compartments ----------------
rows = []
for s in sig_cols:
    wide = pt.pivot(index=PATIENT_COL, columns=TISSUE_COL, values=s)[tissues].dropna()
    chi2, p = friedmanchisquare(*[wide[t].values for t in tissues])
    rows.append({"signature": s, "n_patients": wide.shape[0],
                 "chi2": round(chi2, 2), "p_raw": p})
 
res = pd.DataFrame(rows)
res["FDR"] = multipletests(res["p_raw"], method="fdr_bh")[1]
print("=== Friedman (4 signatures x 4 tissues, paired) ===")
print(res.to_string(index=False))
 
# --- Post-hoc pairwise Wilcoxon signed-rank tests -----------------------------
post = []
for s in sig_cols:
    wide = pt.pivot(index=PATIENT_COL, columns=TISSUE_COL, values=s)[tissues].dropna()
    for i in range(len(tissues)):
        for j in range(i + 1, len(tissues)):
            a, b = tissues[i], tissues[j]
            stat, p = wilcoxon(wide[a], wide[b])
            post.append({"signature": s, "comparison": f"{a} vs {b}",
                         "median_diff": round(np.median(wide[a] - wide[b]), 3),
                         "p_raw": p})
 
post = pd.DataFrame(post)
post["FDR"] = multipletests(post["p_raw"], method="fdr_bh")[1]   # BH across 24
print("\n=== Post-hoc Wilcoxon signed-rank (BH across 24 comparisons) ===")
print(post.round(4).to_string(index=False))
 
res.to_csv(OUTPUT_DIR / "friedman_tissue.csv", index=False)
post.to_csv(OUTPUT_DIR / "posthoc_wilcoxon_tissue.csv", index=False)
 
 
# ============================================================================
# A6.2 Compartment means displayed in Table 8, upper half
# ============================================================================
 
all_cells = adata.obs.groupby(TISSUE_COL, observed=True)[sig_cols].mean()
print("Mean signature score per compartment (all cells):")
print(all_cells.round(2))
 
 
# ============================================================================
# A6.3 Patient/tissue-level boxplots (Figure 4)
# ============================================================================
 
patient_signature_scores = (adata.obs
                            .groupby([PATIENT_COL, TISSUE_COL, STAGE_COL],
                                     observed=True)[sig_cols]
                            .mean().reset_index())
patient_signature_scores.to_csv(
    OUTPUT_DIR / "patient_tissue_stage_signature_scores_mean.csv", index=False)
 
for score in sig_cols:
    plt.figure(figsize=(6, 4))
    patient_signature_scores.boxplot(column=score, by=TISSUE_COL, rot=45)
    plt.title(score + " by Tissue")
    plt.suptitle("")
    plt.ylabel(score)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / f"boxplot_{score}_by_tissue.png",
                dpi=300, bbox_inches="tight")
    plt.close()
 
 
# ============================================================================
# A7. Within-cell-type comparison (Table 8, lower half)
# ============================================================================
 
import pandas as pd
import numpy as np
from scipy.stats import friedmanchisquare, wilcoxon
from statsmodels.stats.multitest import multipletests
 
MIN_CELLS = 10  # supervisor-approved threshold
 
COMP_SIG = "Macrophage_complement_score"
DAM_SIG  = "Disease_associated_macrophage_score"
CD14_SIG = "CD14_monocyte_score"
 
obs = adata.obs.copy()
 
def build_matrix(obs, celltype_mask, score_col, compartments, min_cells=MIN_CELLS):
    """Patient x compartment mean scores, with and without the cell-count threshold."""
    sub = obs[celltype_mask].copy()
    counts = (sub.groupby([PATIENT_COL, TISSUE_COL], observed=True)
                 .size().unstack(TISSUE_COL))
    counts = counts.reindex(columns=compartments)
    means = (sub.groupby([PATIENT_COL, TISSUE_COL], observed=True)[score_col]
                .mean().unstack(TISSUE_COL))
    means = means.reindex(columns=compartments)
    wide_raw = means.dropna(subset=compartments)
    means_thr = means.copy()
    mask_low = counts < min_cells
    means_thr = means_thr.mask(mask_low)
    wide_thr = means_thr.dropna(subset=compartments)
    return wide_thr, wide_raw, counts
 
def paired_tests(wide, compartments, label):
    """Friedman across compartments, with post-hoc Wilcoxon signed-rank tests."""
    n = wide.shape[0]
    print(f"\n--- {label}  (n = {n} patients, threshold >= {MIN_CELLS} cells) ---")
    print("Mean score per compartment:")
    print(wide[compartments].mean().round(3).to_string())
    chi2, p = friedmanchisquare(*[wide[c].values for c in compartments])
    print(f"Friedman: chi2 = {chi2:.2f}, p = {p:.4g}")
    pairs, praw = [], []
    for i in range(len(compartments)):
        for j in range(i + 1, len(compartments)):
            a, b = compartments[i], compartments[j]
            stat, pp = wilcoxon(wide[a], wide[b])
            pairs.append(f"{a} vs {b}")
            praw.append(pp)
    fdr = multipletests(praw, method="fdr_bh")[1]
    print("Post-hoc (Wilcoxon signed-rank, BH):")
    for pr, raw, f in zip(pairs, praw, fdr):
        print(f"   {pr:16s} p = {raw:.4f}   FDR = {f:.4f}")
    return n, chi2, p
 
def highest_count(wide, compartments, focus="LIVER", label=""):
    """How many patients show the highest score in the focus compartment."""
    n = wide.shape[0]
    highest = wide[compartments].idxmax(axis=1)
    n_focus = (highest == focus).sum()
    print(f"   {label}: {focus} highest in {n_focus}/{n} patients")
    for c in compartments:
        if c == focus:
            continue
        n_gt = (wide[focus] > wide[c]).sum()
        print(f"      {focus} > {c}: {n_gt}/{n}")
    return n_focus, n
 
# ===== MACROPHAGES (LIVER / VAT / SAT; PBMC excluded) =====
mac_comp = ["LIVER", "VAT", "SAT"]
mac_mask = obs[CELLTYPE_COL].str.contains("Macrophage", na=False)
 
comp_thr, comp_raw, comp_counts = build_matrix(obs, mac_mask, COMP_SIG, mac_comp)
print("\nMacrophage counts per patient x compartment:")
print(comp_counts.fillna(0).astype(int).to_string())
 
print("\n### COMPLEMENT signature ###")
paired_tests(comp_thr, mac_comp, "Complement within macrophages (thresholded)")
print("Per-patient highest-compartment count:")
highest_count(comp_raw, mac_comp, focus="LIVER", label="UNTHRESHOLDED")
highest_count(comp_thr, mac_comp, focus="LIVER", label="THRESHOLDED ")
 
dam_thr, dam_raw, dam_counts = build_matrix(obs, mac_mask, DAM_SIG, mac_comp)
print("\n### DISEASE-ASSOCIATED signature ###")
paired_tests(dam_thr, mac_comp, "DAM within macrophages (thresholded)")
print("Per-patient highest-compartment count:")
highest_count(dam_raw, mac_comp, focus="SAT", label="UNTHRESHOLDED")
highest_count(dam_thr, mac_comp, focus="SAT", label="THRESHOLDED ")
 
# ===== CD14 MONOCYTES (all four compartments) =====
cd14_comp_full = ["PBMC", "LIVER", "VAT", "SAT"]
cd14_mask = obs[CELLTYPE_COL].str.contains("CD14", na=False)
 
cd14_thr_full, cd14_raw_full, cd14_counts = build_matrix(
    obs, cd14_mask, CD14_SIG, cd14_comp_full)
print("\nCD14 monocyte counts per patient x compartment:")
print(cd14_counts.fillna(0).astype(int).to_string())
 
print("\n### CD14 across 4 compartments ###")
if cd14_thr_full.shape[0] >= 3:
    paired_tests(cd14_thr_full, cd14_comp_full,
                 "CD14 within monocytes, 4 compartments (thresholded)")
else:
    print(f"   Only n={cd14_thr_full.shape[0]} complete across 4 compartments.")
 
print("\nRETAINED n:")
print(f"  Complement (LIVER/VAT/SAT): n = {comp_thr.shape[0]}")
print(f"  DAM        (LIVER/VAT/SAT): n = {dam_thr.shape[0]}")
print(f"  CD14 4-compartment        : n = {cd14_thr_full.shape[0]}")
 
 
# ============================================================================
# A7. Benjamini-Hochberg correction across the three within-cell-type Friedman tests
# ============================================================================
 
from statsmodels.stats.multitest import multipletests
 
# Raw Friedman p values from the three within-cell-type tests above.
# These values are copied from the printed Friedman output above; re-derive
# them if re-running with different data or a different cell-count threshold.
raw_p = [0.005248, 9.6e-5, 2.7e-5]        # complement, DAM, CD14
fdr = multipletests(raw_p, method="fdr_bh")[1]
for name, r, f in zip(["complement", "DAM", "CD14"], raw_p, fdr):
    print(f"{name}: raw p={r:.4g} -> FDR={f:.4g}")
 
 
# ============================================================================
# A8.1 Fibrosis stage, pooled across compartments (Table 10, Figure 5)
# ============================================================================
 
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
 
sig_cols = ["CD14_monocyte_score", "CD16_monocyte_score",
            "Macrophage_complement_score", "Disease_associated_macrophage_score"]
 
# --- Patient-level means across all four compartments, all three groups -------
pat = (adata.obs.groupby(PATIENT_COL, observed=True)
       .agg(StageGroup=(STAGE_COL, "first"),
            **{c: (c, "mean") for c in sig_cols})
       .reset_index())
 
# Group means reported in Table 10, including the healthy control (n = 1)
print("Patients per stage group:")
print(pat["StageGroup"].value_counts().to_string())
print("\nMean patient-level score per stage group:")
print(pat.groupby("StageGroup", observed=True)[sig_cols].mean().round(3).to_string())
 
# --- Mann-Whitney U, F0_1 vs F2_3, healthy control excluded -------------------
pat_two = pat[pat["StageGroup"].isin(["F0_1", "F2_3"])]
 
results = []
for col in sig_cols:
    g1 = pat_two.loc[pat_two["StageGroup"] == "F0_1", col].dropna()
    g2 = pat_two.loc[pat_two["StageGroup"] == "F2_3", col].dropna()
    U, p = mannwhitneyu(g1, g2, alternative="two-sided")
    results.append({"signature": col, "U": U, "p": p,
                    "mean_F0_1": g1.mean(), "mean_F2_3": g2.mean(),
                    "n_F0_1": len(g1), "n_F2_3": len(g2)})
 
df = pd.DataFrame(results)
df["FDR"] = multipletests(df["p"], method="fdr_bh")[1]   # BH across 4 signatures
 
print("\n=== Pooled stage comparison, Mann-Whitney U (healthy excluded) ===")
for _, r in df.iterrows():
    print(f"{r['signature']:40s} U={r['U']:.1f}  p={r['p']:.4f}  FDR={r['FDR']:.4f}  "
          f"(F0_1 {r['mean_F0_1']:.3f} vs F2_3 {r['mean_F2_3']:.3f}; "
          f"n={r['n_F0_1']}/{r['n_F2_3']})")
 
# --- Bootstrap 95% confidence intervals for the difference in means -----------
def boot_ci(a, b, n=5000, seed=0):
    rng = np.random.default_rng(seed)
    d = [rng.choice(a, len(a), True).mean() - rng.choice(b, len(b), True).mean()
         for _ in range(n)]
    return np.percentile(d, [2.5, 97.5])
 
for s in sig_cols:
    a = pat_two.loc[pat_two["StageGroup"] == "F0_1", s].values
    b = pat_two.loc[pat_two["StageGroup"] == "F2_3", s].values
    u, p = mannwhitneyu(a, b)
    rb = 1 - 2 * u / (len(a) * len(b))          # rank-biserial effect size
    lo, hi = boot_ci(a, b)
    print(f"{s:36} diff={a.mean() - b.mean():+.3f} "
          f"95%CI [{lo:+.3f}, {hi:+.3f}] rb={rb:+.2f}")
 
# --- Figure 5: patient-level boxplots by stage group --------------------------
groups = ["F0_1", "F2_3", "H"]
for score in sig_cols:
    fig, ax = plt.subplots(figsize=(6, 4))
    data = [pat.loc[pat["StageGroup"] == g, score].values for g in groups]
    ax.boxplot(data, tick_labels=groups)
    for i, d in enumerate(data, start=1):
        x = np.random.normal(i, 0.05, size=len(d))
        ax.scatter(x, d, alpha=0.7, s=25, color="steelblue", zorder=3)
    ax.set_ylabel(score.replace("_", " "))
    ax.set_xlabel("Fibrosis stage group")
    ax.set_title(score.replace("_", " ") + " by stage group")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / f"boxplot_{score}_by_stagegroup_patientlevel.png",
                dpi=300, bbox_inches="tight")
    plt.close()
 
 
# ============================================================================
# A8.2 Fibrosis stage within each compartment: tests (Table 11)
# ============================================================================
 
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
 
SIGNATURES = ["CD14_monocyte_score", "CD16_monocyte_score",
              "Macrophage_complement_score", "Disease_associated_macrophage_score"]
 
obs = adata.obs.copy()
 
def compare_stage_within(df, label):
    """Aggregate to patient level within a compartment, then Mann-Whitney U.
 
    F0_1 is compared with F2_3; the healthy control (n = 1 patient) is
    excluded from testing but its mean is retained for reference.
    """
    rows = []
    for sig in SIGNATURES:
        pat = (df.groupby([PATIENT_COL, STAGE_COL], observed=True)[sig]
                 .mean().reset_index())
        a = pat.loc[pat[STAGE_COL] == "F0_1", sig].dropna().values
        b = pat.loc[pat[STAGE_COL] == "F2_3", sig].dropna().values
        h = pat.loc[pat[STAGE_COL] == "H",    sig].dropna().values
 
        if len(a) >= 2 and len(b) >= 2:
            stat, p = mannwhitneyu(a, b, alternative="two-sided")
        else:
            stat, p = np.nan, np.nan
 
        rows.append({
            "group": label, "signature": sig,
            "n_patients_F0_1": len(a), "n_patients_F2_3": len(b),
            "n_patients_H": len(h),
            "mean_F0_1": np.mean(a) if len(a) else np.nan,
            "mean_F2_3": np.mean(b) if len(b) else np.nan,
            "mean_H":    np.mean(h) if len(h) else np.nan,
            "direction": ("F0_1 > F2_3" if (len(a) and len(b)
                                            and np.mean(a) > np.mean(b))
                          else "F2_3 > F0_1" if (len(a) and len(b)) else "NA"),
            "MannWhitney_U": stat, "p_raw": p,
        })
    return pd.DataFrame(rows)
 
res = []
for t in sorted(obs[TISSUE_COL].astype(str).unique()):
    sub = obs[obs[TISSUE_COL].astype(str) == t]
    res.append(compare_stage_within(sub, label=t))
res_all = pd.concat(res, ignore_index=True)
 
# BH correction across all 16 comparisons together
ok = res_all["p_raw"].notna()
res_all.loc[ok, "fdr_bh"] = multipletests(res_all.loc[ok, "p_raw"],
                                          method="fdr_bh")[1]
 
print(res_all[["group", "signature", "mean_F0_1", "mean_F2_3",
               "direction", "p_raw", "fdr_bh"]].to_string(index=False))
res_all.to_csv(OUTPUT_DIR / "within_tissue_stage_ALLCELLS.csv", index=False)
 
 
# ============================================================================
# A8.2 Fibrosis stage within each compartment: differences and bootstrap confidence intervals (Table 11)
# ============================================================================
 
# --- Differences and bootstrap 95% CIs for Table 11 ---------------------------
pt = (adata.obs.groupby([PATIENT_COL, TISSUE_COL], observed=True)
      .agg(StageGroup=(STAGE_COL, "first"),
           **{c: (c, "mean") for c in SIGNATURES})
      .reset_index())
 
def boot_ci(a, b, n=5000, seed=0):
    rng = np.random.default_rng(seed)
    d = [rng.choice(a, len(a), True).mean() - rng.choice(b, len(b), True).mean()
         for _ in range(n)]
    return np.percentile(d, [2.5, 97.5])
 
rows = []
for t in ["PBMC", "LIVER", "VAT", "SAT"]:
    sub = pt[pt[TISSUE_COL] == t]
    for s in SIGNATURES:
        a = sub.loc[sub["StageGroup"] == "F0_1", s].values
        b = sub.loc[sub["StageGroup"] == "F2_3", s].values
        if len(a) < 2 or len(b) < 2:
            continue
        u, p = mannwhitneyu(a, b)
        rb = 1 - 2 * u / (len(a) * len(b))       # rank-biserial effect size
        lo, hi = boot_ci(a, b)
        rows.append({"tissue": t, "signature": s,
                     "mean_F0_1": round(a.mean(), 3),
                     "mean_F2_3": round(b.mean(), 3),
                     "diff": round(a.mean() - b.mean(), 3),
                     "CI_low": round(lo, 3), "CI_high": round(hi, 3),
                     "rank_biserial": round(rb, 2),
                     "p": round(p, 3)})
 
effects = pd.DataFrame(rows)
print(effects.to_string(index=False))
effects.to_csv(OUTPUT_DIR / "effect_sizes_within_compartment.csv", index=False)
 
 
# ============================================================================
# A8.2 Fibrosis stage within each compartment: boxplots (Figure 7)
# ============================================================================
 
# --- Figure 7: signature scores by stage group, faceted by compartment --------
plot_df = (obs.groupby([PATIENT_COL, TISSUE_COL, STAGE_COL],
                       observed=True)[SIGNATURES].mean().reset_index())
 
tissue_order = [t for t in ["PBMC", "LIVER", "VAT", "SAT"]
                if t in plot_df[TISSUE_COL].astype(str).unique()]
stage_order = [s for s in ["F0_1", "F2_3", "H"]
               if s in plot_df[STAGE_COL].astype(str).unique()]
 
for sig in SIGNATURES:
    fig, axes = plt.subplots(1, len(tissue_order),
                             figsize=(4 * len(tissue_order), 4), sharey=True)
    if len(tissue_order) == 1:
        axes = [axes]
    for ax, t in zip(axes, tissue_order):
        d = plot_df[plot_df[TISSUE_COL].astype(str) == t]
        sns.boxplot(data=d, x=STAGE_COL, y=sig, order=stage_order, ax=ax,
                    showfliers=False)
        sns.stripplot(data=d, x=STAGE_COL, y=sig, order=stage_order, ax=ax,
                      color="black", size=4, alpha=0.6)
        ax.set_title(t)
        ax.set_xlabel("")
        if ax is not axes[0]:
            ax.set_ylabel("")
    fig.suptitle(f"{sig} by fibrosis stage, within each tissue", y=1.02)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / f"boxplot_{sig}_by_stage_within_tissue.png",
                dpi=200, bbox_inches="tight")
    plt.close()
 
 
# ============================================================================
# A10. Sensitivity analysis: exclusion of non-myeloid subclusters (Section 3.3)
# ============================================================================
 
sig_cols = ["CD14_monocyte_score", "CD16_monocyte_score",
            "Macrophage_complement_score", "Disease_associated_macrophage_score"]
contam = ["1", "9", "11"]     # C1 lymphoid, C9 erythroid, C11 platelet
 
# Same computation as Section A6.2, repeated here so that this sensitivity
# analysis is self-contained and gives the baseline for the comparison below.
all_cells = adata.obs.groupby(TISSUE_COL, observed=True)[sig_cols].mean()
 
clean = adata.obs[~adata.obs[CLUSTER_COL].isin(contam)]
clean_mean = clean.groupby(TISSUE_COL, observed=True)[sig_cols].mean()
 
print("After exclusion:")
print(clean_mean.round(3).to_string())
 
diff = (clean_mean - all_cells).abs()
print(f"Largest change in any compartment mean: {diff.values.max():.3f}")
 
for c in sig_cols:
    print(f"{c}: highest all cells = {all_cells[c].idxmax()}, "
          f"highest after exclusion = {clean_mean[c].idxmax()}")
 
