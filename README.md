# MASLD Monocyte and Macrophage Single-Cell Analysis

Secondary analysis of a paired single-cell RNA-seq (CITE-seq) dataset examining
monocyte and macrophage heterogeneity across four compartments in metabolic
dysfunction-associated steatotic liver disease (MASLD).

## Overview

This repository contains the analysis code for an MSc dissertation investigating
how monocyte and macrophage transcriptional states vary across peripheral blood
(PBMC), liver, visceral adipose tissue (VAT) and subcutaneous adipose tissue
(SAT), and across fibrosis-stage groups, in 20 patients sampled across all four
compartments (a fully paired design).

The main finding is that tissue compartment, rather than fibrosis stage, is the
dominant determinant of monocyte and macrophage transcriptional state.

## Repository structure

The analysis is split into six scripts that run **in numerical order, 01 to 06**.
Scripts 02 and 03 write intermediate `.h5ad` files that later scripts read, so
they must be run before 04, 05 and 06.

| Script | Purpose | Reads | Writes |
|--------|---------|-------|--------|
| `01_quality_control.py` | Inspect existing QC metrics (counts, genes, mitochondrial and haemoglobin content) | `data/mono_macro.h5ad` | QC tables and plots |
| `02_subclustering.py` | Leiden subclustering on the Harmony embedding; marker validation and differential expression | `data/mono_macro.h5ad` | `mono_macro_with_subclusters.h5ad` |
| `03_signature_scoring.py` | Score four transcriptional programmes with `scanpy.tl.score_genes` | `mono_macro_with_subclusters.h5ad` | `mono_macro_with_scores.h5ad` |
| `04_doublet_detection.py` | Objective doublet assessment with Scrublet, per sequencing run | `mono_macro_with_subclusters.h5ad` | Doublet summary table |
| `05_statistical_analysis.py` | Friedman and post-hoc Wilcoxon (paired compartments); Kruskal-Wallis and Mann-Whitney U (stage groups); effect sizes and bootstrap confidence intervals | `mono_macro_with_scores.h5ad` | Statistical result tables |
| `06_downsampling.py` | Sequencing-depth sensitivity analysis | `mono_macro_with_scores.h5ad` | Depth-matched Friedman results |

## Expected input

The pipeline expects a processed AnnData object at `data/mono_macro.h5ad` with
the following structure:

**Matrices**

| Slot | Contents |
|------|----------|
| `.X` | log1p-transformed, library-size-normalised expression (target sum 10,000 per cell) |
| `.layers["counts_RNA"]` | raw integer counts (required by scripts 04 and 06) |
| `.obsm["X_harmony"]` | pre-computed 50-dimensional Harmony embedding (required by script 02) |

**Required `.obs` columns**

| Column | Contents |
|--------|----------|
| `Patient_ID` | patient identifier (20 patients) |
| `Sample` | sequencing-run identifier, used as the batch key and by Scrublet |
| `Tissue` | one of `PBMC`, `LIVER`, `VAT`, `SAT` |
| `StageGroup` | fibrosis-stage group: `F0_1`, `F2_3` or `H` |
| `cell_type_with_cluster` | existing annotation (CD14 monocytes, CD16 monocytes, macrophages) |
| `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `percent.hb` | quality-control metrics (used by script 01) |

Scripts 02 and 03 add further columns (`leiden_*` and the four `*_score`
columns) that downstream scripts rely on.

## Installation

```bash
pip install -r requirements.txt
```

Pinned versions are listed in `requirements.txt`. The analysis was run on
Python 3.11.

## How to run

1. Place the input object at `data/mono_macro.h5ad` (see **Expected input** above).
2. Run the scripts in numerical order:

```bash
python 01_quality_control.py
python 02_subclustering.py
python 03_signature_scoring.py
python 04_doublet_detection.py
python 05_statistical_analysis.py
python 06_downsampling.py
```

Scripts 02 and 03 must be run before 04, 05 and 06, as they write the
intermediate `.h5ad` files that those scripts read. Results are written to
`outputs/`.

## Data availability

The raw data are not included in this repository, as they are held by the host
laboratory and subject to data-governance and ethical restrictions. This
repository contains analysis code only. The `.obs` schema above documents the
input format so that the code can be inspected and adapted.

## Author

Yeqin Huang — MSc dissertation, [your institution/department]
