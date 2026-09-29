# Single-cell RNA-seq analysis for “Targeting a Treg-driven multicellular relay limits pathological immunosuppression”

## Overview

This repository contains code for two complementary single-cell analyses investigating the role of Treg-derived, GARP-mediated TGF-β1 activity in lung cancer progression:

1. Single-cell RNA sequencing (scRNA-seq): Analysis of lung immune cells from control and Lrrc32^fl/fl Foxp3-Cre⁺ mice inoculated with Lewis lung carcinoma (LLC) cells.

2. Single-cell multiome sequencing (scRNA-seq and scATAC-seq): Multiomic analysis of gene expression and chromatin accessibility in immune cells from the lungs of wild-type (WT) mice inoculated with LLC cells with receiving control treatment or MERIT treatment.

## Data availability

Sequencing data are available at DDBJ under accession number xxxxxxxx.

## Analysis scripts

Run the scripts in the following order:

1. `01_preprocessing.R` — Quality control and normalization.
2. `02_clustering.R` — Dimensionality reduction and clustering.
3. `03_annotation.R` — Cell-type annotation.
4. `04_downstream_analysis.R` — Differential expression, cell-cell interaction and visualization.

## Usage

The analyses were performed using R , Seurat , and [other key packages and versions].

## Contact
    
For questions, please contact gou-im@m.u-tokyo.ac.jp.
