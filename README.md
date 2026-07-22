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

## Dataset

The analysis uses a processed AnnData object (`.h5ad`) provided by the host
laboratory, containing 22,983 monocytes and macrophages (36,601 genes) annotated
into CD14 monocytes, CD16 monocytes and macrophages.

**Note:** The raw data are not included in this repository, as they are held by
the host laboratory and subject to data-governance and ethical restrictions.
This repository contains analysis code only.

## Repository structure

The analysis is split into scripts that run in order. Each script reads the
output of the previous one.

| Script | Purpose |
|--------|---------|
| `01_quality_control.py` | Inspect existing QC metrics (counts, genes, mitochondrial and haemoglobin content) |
| `02_subclustering.py` | Leiden subclustering on the Harmony embedding; marker validation and differential expression |
| `03_signature_scoring.py` | Score four transcriptional programmes with `scanpy.tl.score_genes` |
| `04_doublet_detection.py` | Objective doublet assessment with Scrublet (per sequencing run) |
| `05_statistical_analysis.py` | Friedman + post-hoc Wilcoxon (paired compartments), Kruskal-Wallis and Mann-Whitney U (stage), effect sizes and bootstrap CIs |
| `06_downsampling.py` | Sequencing-depth sensitivity analysis |

## Requirements

- Python 3.11
- scanpy
- pandas, numpy, scipy
- statsmodels
- scikit-image (required for Scrublet in `04_doublet_detection.py`)
- matplotlib
- leidenalg, igraph (required for Leiden clustering)

## How to run

1. Place the `.h5ad` data object in a `data/` folder.
2. Update the paths at the top of each script if needed.
3. Run the scripts in numerical order (`01` to `06`). Scripts `02` and `03`
   write intermediate `.h5ad` files that later scripts read.

## Author

Yeqin Huang — MSc dissertation, [your institution/department]
