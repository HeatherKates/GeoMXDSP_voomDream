HuBMAP Pancreas Spatial Transcriptomics Analysis
This R-based workflow analyzes GeoMX DSP spatial transcriptomics data from the HuBMAP Pancreas project, focusing on differences between islet-present and islet-absent regions across different pancreatic cell types.

Input Data
File: P1-P4.QC.v3.xlsx
Exported from GeoMX DSP Analysis Suite, containing:

BioProbeCountMatrix: Raw target-level probe counts

SegmentProperties: Segment-level sample metadata

Folder: ../data/
Contains the .xlsx input and any additional needed files.

Required R Packages

library(standR)
library(SpatialExperiment)
library(edgeR)
library(limma)
library(variancePartition)
library(BiocParallel)
library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(msigdbr)

Workflow Summary
1. Data Import & Preprocessing
Import raw counts, metadata, and feature annotations using a custom geomx_import_fun_HRK() helper.

Create a SpatialExperiment object.

Filter for pancreas samples and remove low-quality segments and genes.

2. Normalization & Log Transformation
Apply log2(count + 1) transformation.

Batch correction with limma::removeBatchEffect() (using Case and SlideName).

3. PCA Visualization
PCA plots (corrected and uncorrected) colored by:

ROI_type (islet_present / islet_absent)

organ_region

Split by Case and AOI_target

4. Differential Expression (DE) Analysis
For each AOI_target:

Use edgeR for TMM normalization.

Apply voomWithDreamWeights() and dream() from variancePartition.

Fit mixed-effects model:
~ ROI_type + (1 | Case/SlideName)

5. DE Filtering & Annotation
Filter for |logFC| > 0.5 and adj.P.Val < 0.05

Annotate genes using MSigDB C8 cell-type signatures

Categorize genes into: endocrine, exocrine, endothelial, or unknown

6. Visualization
Generate high-res heatmaps per AOI_target (e.g., endothelial_cells)

Expression plots for genes of interest across ROI types

Output Files
PCA_by_*.png: PCA plots (corrected and uncorrected)

heatmap_*.png: Heatmaps of differentially expressed genes for each AOI target

geomx_pancreas_analysis_workspace.RData: Saved R session

Example Plot Function

plot_gene_by_ROItype(spe, "SCG5")

Notes
The "NegProbe-WTX" probe is retained but flagged during preprocessing.

AOIs are grouped into islet-present and islet-absent based on manual annotation and segment-level metadata.

Cell type assignment is based on MSigDB C8 collections with custom parsing and grouping logic.
