# MASLD Monocyte and Macrophage Single-Cell Analysis

Secondary single-cell analysis of monocyte and macrophage heterogeneity across liver, blood, visceral adipose and subcutaneous adipose tissue in MASLD-associated early fibrosis.

This repository contains the analysis code for the MSc dissertation *Dissecting Monocyte and Macrophage Heterogeneity Across Liver, Blood, and Adipose Tissues in MASLD-Associated Early Fibrosis* (MSc Gastroenterology, Barts and The London School of Medicine and Dentistry, Queen Mary University of London). The same code is reproduced in the Supplementary Appendix of the dissertation, which is the authoritative version; each script below is annotated with the appendix sections it corresponds to.

## Data

The dataset is a CITE-seq object generated and processed by the Alazawi Laboratory. It is not publicly available, for reasons of data governance and participant confidentiality, and is not included in this repository. This repository contains analysis code only.

The scripts expect a processed AnnData object in which `.X` holds log1p-transformed, library-size-normalised values (target sum 10,000), raw counts are retained in the `counts_RNA` layer, and a Harmony embedding and an existing cell-type annotation are already present. Normalisation, integration and the original annotation were performed upstream by the host laboratory and are not re-run here; `01_quality_control.py` only verifies them.

The path to the object is left as a placeholder (`INPUT_OBJECT`) at the top of `01_quality_control.py`, alongside `OUTPUT_DIR` for the tables and figures. Both must be supplied by the user; no path to the data is recorded in this repository.

## Environment

Python 3.11. Install dependencies with:

pip install -r requirements.txt


## How to run

These scripts are not independent programs. They operate on a single AnnData object held in memory and are meant to be executed one after another inside a single Python session, each continuing from the state left by the previous one. Later scripts depend on columns added by earlier ones — the Leiden labels from `02`, the signature scores from `03` — so no script other than `01` loads the object, and re-loading it midway would discard those columns.

One way to run the sequence:

python -i 01_quality_control.py


then, in the same interpreter:

```python
for f in ["02_subclustering.py",
          "03_signature_scoring.py",
          "04_doublet_detection.py",
          "05_statistical_analysis.py",
          "06_downsampling.py"]:
    exec(open(f).read())
```

Running them as separate `python 0X_....py` invocations will fail, because each new process starts without `adata`.

## Dependency order

01 (load, QC, annotation validation)
└── 02 (Leiden subclustering, marker genes)
├── 03 (signature scoring)
│ ├── 05 (statistical analysis)
│ └── 06 (downsampling sensitivity analysis)
└── 04 (doublet detection)


Run order: 01 → 02 → 03 → 04 → 05 → 06. `04` requires the subcluster labels from `02` and is independent of `03`, so it may be run at any point after `02`.

## Scripts

| Script | What it does | Dissertation output | Appendix |
|---|---|---|---|
| `01_quality_control.py` | Loads the object, defines the metadata column names used throughout, summarises dataset composition and the patient cohort, reports QC metrics, verifies the upstream normalisation, and validates the existing annotation with canonical markers and UMAP projections. | Tables 1, 2 (lower part), 3, 6; Figures 1A, 1B, 2, 6 | A0, A1.1–A1.5, A1.7, A2 |
| `02_subclustering.py` | Builds a 15-nearest-neighbour graph on the Harmony embedding, runs Leiden clustering at six resolutions (0.5 used downstream), assembles the subcluster composition table, and ranks marker genes at both the broad cell-type and subcluster level. | Table 7; Figures 3A, 3B | A3.1–A3.3, A4.1–A4.3 |
| `03_signature_scoring.py` | Scores the four gene signatures with `sc.tl.score_genes` and summarises them by subcluster, cell type, compartment and fibrosis-stage group. | Table 5 (composition and settings) | A5 |
| `04_doublet_detection.py` | Runs Scrublet per sequencing run and summarises doublet scores by subcluster. Requires the subcluster labels from `02`; independent of `03`. | Table 4 | A1.6 |
| `05_statistical_analysis.py` | Cross-compartment comparison (Friedman with post-hoc Wilcoxon signed-rank), the within-cell-type comparison under a minimum-cell threshold, the fibrosis-stage comparisons (Mann–Whitney U), and the sensitivity analysis excluding non-myeloid subclusters. | Tables 8, 9, 10, 11; Figures 4, 5, 7 | A6, A7, A8, A10 |
| `06_downsampling.py` | Sequencing-depth downsampling sensitivity analysis, repeating the cross-compartment tests at a common per-cell depth. | Section 4.5 | A9 |

## Statistical approach

Tissue compartments were sampled from the same patients and were therefore compared with Friedman tests on patient-level means, followed by pairwise Wilcoxon signed-rank post-hoc tests. Fibrosis-stage groups comprise different patients and were compared with Mann–Whitney U tests. The single healthy control was excluded from all testing and is shown for reference only.

Benjamini–Hochberg correction was applied within each family of tests:

- across the four signatures for the all-cell Friedman tests;
- across the 24 post-hoc compartment comparisons (4 signatures × 6 pairs);
- across the three within-cell-type Friedman tests, with their post-hoc comparisons corrected separately within each cell type;
- across the four signatures for the pooled fibrosis-stage comparison;
- across the 16 within-compartment fibrosis-stage comparisons (4 compartments × 4 signatures).

The within-cell-type analysis requires a minimum of 10 cells per patient per compartment, so that each patient contributes a stable per-cell estimate. Macrophages were compared across liver, VAT and SAT only, blood being excluded because it contained very few macrophages; CD14 monocytes were compared across all four compartments. Full details are given in the dissertation.

## Use of AI

An AI assistant was used as a coding assistant to help write, debug and refactor the Python code in this repository. The analysis itself was carried out by established deterministic tools (Scanpy, scikit-image, SciPy, statsmodels); no AI or machine-learning model was used to analyse the data, to generate or predict any value, or to build any predictive model. The analysis was designed, run and verified by the author. The full statement is given in the Supplementary Appendix of the dissertation.

## Author

Yeqin Huang — MSc dissertation, qmul
