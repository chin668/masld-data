# MASLD Monocyte and Macrophage Single-Cell Analysis

Secondary analysis of a paired single-cell RNA-seq (CITE-seq) dataset examining
monocyte and macrophage heterogeneity across four compartments in metabolic
dysfunction-associated steatotic liver disease (MASLD).

## Overview

This repository contains the analysis code for an MSc dissertation investigating
how monocyte and macrophage transcriptional states vary across peripheral blood
(PBMC), liver, visceral adipose tissue (VAT) and subcutaneous adipose tissue (SAT),
and across fibrosis-stage groups, in 20 patients sampled across all four
compartments (a fully paired design).

The main finding is that tissue compartment, rather than fibrosis stage, is the
dominant determinant of monocyte and macrophage transcriptional state.

## Dataset

The analysis uses a processed AnnData object (`.h5ad`) provided by the host
laboratory, containing 22,983 monocytes and macrophages (36,601 genes) annotated
into CD14 monocytes, CD16 monocytes and macrophages.

**Note:** The raw data are not included in this repository, as they are held by
the host laboratory and subject to data-governance and ethical restrictions. This
repository contains analysis code only.

## Analysis workflow

The notebook performs the following steps:

1. **Quality control** — inspection of existing QC metrics (counts, genes,
   mitochondrial and haemoglobin content)
2. **Doublet detection** — Scrublet applied per sequencing run
3. **Subclustering** — Leiden clustering (resolution 0.5) on the pre-computed
   Harmony embedding, yielding 13 subclusters
4. **Annotation validation** — canonical marker expression and differential
   expression testing (Wilcoxon rank-sum)
5. **Signature scoring** — four transcriptional programmes scored with
   `scanpy.tl.score_genes` (CD14 monocyte, CD16 monocyte, complement-associated
   and disease-associated macrophage)
6. **Statistical analysis** — Friedman tests for paired compartment comparisons,
   post-hoc Wilcoxon signed-rank tests, Kruskal-Wallis and Mann-Whitney U tests
   for stage comparisons, with Benjamini-Hochberg correction and bootstrap
   confidence intervals
7. **Sensitivity analysis** — downsampling to a common sequencing depth

## Requirements

- Python 3.11
- scanpy
- pandas, numpy, scipy
- statsmodels
- scikit-image (required for Scrublet)
- matplotlib

## How to run

1. Place the `.h5ad` data object in the working directory.
2. Open the notebook in Jupyter.
3. Update the data file path at the top of the notebook.
4. Run the cells in order.

## Author

Yeqin Huang — MSc dissertation, [Blizard institute]
