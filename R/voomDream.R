# Load required packages
library(standR)
library(SpatialExperiment)
library(ggplot2)
library(dplyr)
library(limma)
library(edgeR)
library(readxl)
library(dplyr)
source("helpers.R")
select <- dplyr::select
data_path="../data/"
data_file="P1-P4.QC.v3.xlsx"

# Get the sheet names
sheet_names <- excel_sheets(paste(data_path,data_file,sep="/"))

# Read each sheet and assign it to a data frame named after the sheet's name
for (sheet in sheet_names) {
  assign(sheet, read_excel(paste(data_path,data_file,sep="/"), sheet = sheet))
}

# Create the count data frame. Samples in columns and features/genes in rows. The first column is the gene names/ids

# Create the sample annotation data frame.
sampleAnno <- SegmentProperties

# 1. Create a unique ID for each probe row (including duplicates like "NegProbe-WTX")
BioProbeCountMatrix <- BioProbeCountMatrix %>%
  mutate(unique_id = make.unique(TargetName))

# 2. Extract feature annotation — keep all columns *including* unique_id
featureAnno <- BioProbeCountMatrix %>%
  select(-all_of(SegmentProperties$SegmentDisplayName))  # drop sample columns
rownames(featureAnno) <- featureAnno$unique_id

# 3. Extract count matrix — keep unique_id + sample columns only
counts <- BioProbeCountMatrix %>%
  select(unique_id, all_of(SegmentProperties$SegmentDisplayName))

# 4. Sample metadata
sampleAnno <- SegmentProperties

# 5. Run the function
spe <- geomx_import_fun_HRK(
  countFile = counts,
  sampleAnnoFile = sampleAnno,
  featureAnnoFile = featureAnno,
  rmNegProbe = FALSE,
  NegProbeName = "NegProbe-WTX",
  colnames.as.rownames = c("unique_id", "SegmentDisplayName", "unique_id"),
  coord.colnames = c("ROICoordinateX", "ROICoordinateY")
)
# Flag low quality ROIs and probes (optional but helps)


# 1. Filter for Pancreas organ
spe <- spe[, colData(spe)$organ == "Pancreas"]
spe <- spe[, !colData(spe)$Case == "P1"]

# 2. Quality Control

## Segment QC: Remove segments with low reads or other QC metrics
# Define your QC thresholds based on your data characteristics
min_reads <- 1000  # Example threshold
spe <- spe[, colData(spe)$RawReads >= min_reads]

## Gene QC: Remove genes with low expression across samples
# Define your QC thresholds based on your data characteristics
min_counts <- 10  # Example threshold
min_samples <- 5  # Example threshold
keep_genes <- rowSums(assay(spe) >= min_counts) >= min_samples
spe <- spe[keep_genes, ]
spe <- spe[, colData(spe)$AOI_target != "beta_cells"]
# 3. PCA Visualization

plot_pca <- function(spe_object, 
                     title = "", 
                     color_by = "ROI_type", 
                     split_by1 = NULL, 
                     split_by2 = NULL, 
                     logcounts_assay = "logcounts") {
  library(ggplot2)
  library(SummarizedExperiment)
  
  # Extract the logcounts matrix from the specified assay
  log_mat <- assay(spe_object, logcounts_assay)
  
  # Run PCA
  pca <- prcomp(t(log_mat), center = TRUE, scale. = TRUE)
  
  # Calculate % variance explained
  percentVar <- (pca$sdev)^2 / sum(pca$sdev^2)
  x_lab <- paste0("PC1 (", round(100 * percentVar[1], 1), "%)")
  y_lab <- paste0("PC2 (", round(100 * percentVar[2], 1), "%)")
  
  # Extract metadata
  metadata <- as.data.frame(colData(spe_object))
  
  # Build PCA dataframe
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    color = metadata[[color_by]]
  )
  
  # Add split variables if provided
  if (!is.null(split_by1)) pca_df$split1 <- metadata[[split_by1]]
  if (!is.null(split_by2)) pca_df$split2 <- metadata[[split_by2]]
  
  # Base plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = color)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(title = title, color = color_by, x = x_lab, y = y_lab) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "grey50", fill = NA, size = 0.5),  # visible borders
      panel.spacing = unit(0.5, "lines"),  # spacing between panels
      strip.background = element_rect(fill = "#f0f0f0", color = NA)
    )
  
  # Apply faceting
  if (!is.null(split_by1) & !is.null(split_by2)) {
    p <- p + facet_grid(split1 ~ split2)
  } else if (!is.null(split_by1)) {
    p <- p + facet_wrap(~ split1)
  }
  
  return(p)
}
# Define combinations of settings
pca_settings <- list(
  list(filename = "PCA_by_ROItype_corrected.png",
       assay = "logcounts_corrected",
       color_by = "ROI_type",
       title = "PCA by ROI_type (corrected)"),
  
  list(filename = "PCA_by_organregion_corrected.png",
       assay = "logcounts_corrected",
       color_by = "organ_region",
       title = "PCA by organ_region (corrected)"),
  
  list(filename = "PCA_by_ROItype_uncorrected.png",
       assay = "logcounts",
       color_by = "ROI_type",
       title = "PCA by ROI_type (uncorrected)"),
  
  list(filename = "PCA_by_organregion_uncorrected.png",
       assay = "logcounts",
       color_by = "organ_region",
       title = "PCA by organ_region (uncorrected)")
)

# Loop through and save each plot
for (setting in pca_settings) {
  png(filename = setting$filename, width = 3000, height = 2500, res = 350)
  
  print(plot_pca(spe,
                 title = setting$title,
                 color_by = setting$color_by,
                 split_by1 = "Case",
                 split_by2 = "AOI_target",
                 logcounts_assay = setting$assay))
  
  dev.off()
}

# Log-transform counts if not already
logcounts(spe) <- log2(counts(spe) + 1)

# Run PCA
pca <- prcomp(t(logcounts(spe)), scale. = TRUE)

# Extract metadata
metadata <- as.data.frame(colData(spe))

# Build PCA dataframe
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],  # You can include more PCs if needed
  ROI_type = metadata$ROI_type,
  AOI_target = metadata$AOI_target,
  organ_region = metadata$organ_region,
  SlideName = metadata$SlideName,
  Case = metadata$Case
)
summary(lm(PC2 ~ ROI_type, data = pca_df))

# 4. Differential Expression Analysis with dream (mixed model)
library(variancePartition)
library(edgeR)
library(limma)
library(BiocParallel)
register(MulticoreParam(workers = 8))

#subset_genes <- rownames(counts(spe))[1:100]  # or sample()
#spe_subset_small <- spe[subset_genes, ]

# Ensure logcounts exists
logcounts(spe) <- log2(counts(spe) + 1)

# Initialize list to hold DE results
de_results <- list()

# Loop over AOI_targets
unique_targets <- unique(colData(spe)$AOI_target)

for (target in unique_targets) {
  
  # Subset for current AOI_target
  spe_subset <- spe[, colData(spe)$AOI_target == target]
  meta <- as.data.frame(colData(spe_subset))
  
  # Skip target if only one level of ROI_type (e.g. beta_cells)
  if (length(unique(meta$ROI_type)) < 2) next
  
  # Create DGEList and calculate normalization
  dge <- DGEList(counts = counts(spe_subset))
  dge <- calcNormFactors(dge)
  
  # Define formula: fixed = ROI_type, random = (1|Case/SlideName)
  formula <- ~ ROI_type + (1 | Case/SlideName)
  
  # Voom with dream
  vobj <- voomWithDreamWeights(dge, formula, meta)
  
  # Fit mixed model
  fit <- dream(vobj, formula, meta)
  
  # Extract top results for ROI_type
  top <- topTable(fit, coef = "ROI_typeislet_present", number = Inf)
  
  # Store
  de_results[[target]] <- top
}

# Ensure logcounts are present
if (is.null(logcounts(spe))) {
  logcounts(spe) <- log2(counts(spe) + 1)
}

plot_gene_by_ROItype <- function(spe_object, gene_symbol) {
  # Check if logcounts are present, compute if needed
  if (is.null(logcounts(spe_object))) {
    logcounts(spe_object) <- log2(counts(spe_object) + 1)
  }
  
  # Check if gene exists
  if (!gene_symbol %in% rownames(spe_object)) {
    stop(paste("Gene", gene_symbol, "not found in the object"))
  }
  
  # Extract expression and metadata
  df <- data.frame(
    expression = as.numeric(logcounts(spe_object)[gene_symbol, ]),
    ROI_type = colData(spe_object)$ROI_type,
    AOI_target = colData(spe_object)$AOI_target
  )
  
  # Plot
  ggplot(df, aes(x = ROI_type, y = expression, fill = ROI_type)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
    facet_wrap(~ AOI_target) +
    labs(
      title = paste(gene_symbol, "expression by ROI type"),
      x = "ROI Type",
      y = paste0("log2(", gene_symbol, " expression + 1)")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
}
plot_gene_by_ROItype(spe, "SCG5")

#plot_gene_by_ROItype(spe, "SIPA1L3")

save.image(file = "geomx_pancreas_analysis_workspace.RData")
load("geomx_pancreas_analysis_workspace.RData")

# Filter and plot
# --- STEP 1: Batch-correct logcounts ---
meta <- as.data.frame(colData(spe))
log_expr <- logcounts(spe)

log_expr_corrected <- limma::removeBatchEffect(
  log_expr,
  batch = meta$Case,
  batch2 = meta$SlideName
)

# Store corrected matrix in spe for future use
assay(spe, "logcounts_corrected") <- log_expr_corrected

filter_de_results <- function(de_results, spe, 
                               logcounts_assay = "logcounts_corrected",
                               logfc_threshold = 0.5,
                               pval_threshold = 0.05) {
  results_filtered <- list()
  
  for (target in names(de_results)) {
    res <- de_results[[target]]
    
    # Skip if result is NULL (e.g., skipped during DE)
    if (is.null(res)) next
    
    # Apply filters
    res <- res[abs(res$logFC) > logfc_threshold & res$adj.P.Val < pval_threshold, , drop = FALSE]
    if (nrow(res) == 0) next  # skip if nothing passes
    
    # Get relevant samples for that AOI_target
    idx <- colData(spe)$AOI_target == target
    roi <- colData(spe)$ROI_type[idx]
    expr_mat <- assay(spe, logcounts_assay)[rownames(res), idx]
    
    # Average expression in each group
    res$expr_islet_present <- rowMeans(expr_mat[, roi == "islet_present", drop = FALSE])
    res$expr_islet_absent  <- rowMeans(expr_mat[, roi == "islet_absent", drop = FALSE])
    
    # Store
    results_filtered[[target]] <- res
  }
  
  return(results_filtered)
}


filtered_de_results <- filter_de_results(de_results, spe)

library(msigdbr)
library(dplyr)
library(stringr)

annotate_de_genes_with_celltype <- function(de_result, species = "Homo sapiens", category = "C8") {
  # 1. Get MSigDB C8 Cell Type signatures
  msigdb_genes <- msigdbr(species = species, category = category)
  
  # 2. Simplify cell type names from gs_name
  msigdb_filtered <- msigdb_genes %>%
    filter(str_detect(gs_name, "PANCREAS|PANCREATIC")) %>%
    mutate(celltype_simple = str_replace(gs_name, ".*?(PANCREAS.*|PANCREATIC.*)", "\\1")) %>%
    select(gene_symbol, celltype_simple)
  
  # 3. Summarize mapping: some genes may map to multiple cell types
  msigdb_map <- msigdb_filtered %>%
    group_by(gene_symbol) %>%
    summarize(msigdb_celltype = paste(unique(celltype_simple), collapse = "; "), .groups = "drop")
  
  # 4. Join to your DE result by gene symbol (rownames)
  annotated_de_result <- de_result %>%
    tibble::rownames_to_column("gene") %>%
    left_join(msigdb_map, by = c("gene" = "gene_symbol")) %>%
    tibble::column_to_rownames("gene")  # Re-assign rownames for heatmap compatibility
  
  return(annotated_de_result)
}

# Example use:
filtered_de_results <- lapply(filtered_de_results, annotate_de_genes_with_celltype)
categorize_celltypes <- function(de_annotated_df) {
  de_annotated_df <- de_annotated_df %>%
    mutate(
      celltype_group = case_when(
        # If endocrine-related terms exist but no exocrine
        str_detect(msigdb_celltype, "ISLET|ALPHA|BETA|DELTA|POLYPEPTIDE") & 
          !str_detect(msigdb_celltype, "DUCTAL|MESENCHYMAL|STROMAL") ~ "endocrine",
        
        # If exocrine terms only
        str_detect(msigdb_celltype, "DUCTAL|MESENCHYMAL|STROMAL") & 
          !str_detect(msigdb_celltype, "ISLET|ALPHA|BETA|DELTA|POLYPEPTIDE") ~ "exocrine",
        
        # If both present
        str_detect(msigdb_celltype, "ISLET|ALPHA|BETA|DELTA|POLYPEPTIDE") & 
          str_detect(msigdb_celltype, "DUCTAL|MESENCHYMAL|STROMAL") ~ "both",
        
        # If endothelial
        str_detect(msigdb_celltype, "ENDOTHELIAL") ~ "endothelial",
        
        # Otherwise unknown
        TRUE ~ "unknown"
      )
    )
  
  return(de_annotated_df)
}

filtered_de_results <- lapply(filtered_de_results, categorize_celltypes)

 library(pheatmap)
 library(SpatialExperiment)
 
library(pheatmap)

plot_de_heatmap <- function(spe, de_result, target_name, 
                            logcounts_assay = "logcounts_corrected") {
  # Subset to AOI target
  target_idx <- colData(spe)$AOI_target == target_name
  spe_target <- spe[, target_idx]
  
  # Get DE genes for this target
  de_genes <- rownames(de_result[[target_name]])
  
  # Subset expression matrix
  expr_mat <- assay(spe_target, logcounts_assay)
  expr_mat <- expr_mat[rownames(expr_mat) %in% de_genes, ]
  
  # Scale rows
  expr_scaled <- t(scale(t(expr_mat)))
  
  # Row annotation
  de_info <- de_result[[target_name]]
  row_annot <- data.frame(CellType = de_info$celltype_group)
  rownames(row_annot) <- rownames(de_info)
  row_annot <- row_annot[rownames(expr_scaled), , drop = FALSE]
  
  # Column annotation
  col_annot <- data.frame(
    ROI_type = colData(spe_target)$ROI_type,
    Case = colData(spe_target)$Case
  )
  rownames(col_annot) <- colnames(expr_scaled)
  
  # Define clear color schemes
  annotation_colors <- list(
    ROI_type = c(
      "islet_present" = "#377eb8",  # Blue
      "islet_absent" = "#e41a1c"    # Red
    ),
    Case = c(
      "P2" = "#4daf4a",  # Green
      "P3" = "#984ea3",  # Purple
      "P4" = "#ff7f00"   # Orange
    ),
    CellType = c(
      "endocrine" = "#f781bf",   # Pink
      "exocrine" = "#a65628",    # Brown
      "endothelial" = "#999999", # Gray
      "unknown" = "#cccccc"      # Light gray
    )
  )
  
  # Plot
  pheatmap::pheatmap(
    mat = expr_scaled,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_row = row_annot,
    annotation_col = col_annot,
    annotation_colors = annotation_colors,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    fontsize_col = 6,
    main = paste0(target_name,": DE Genes by Proximity to Islet"),
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100)  # Better heatmap gradient
  )
}

# Loop over each AOI_target and generate high-res heatmap
for (target_name in names(filtered_de_results)) {
  # Define output filename
  outfile <- paste0("heatmap_", target_name, ".png")
  
  # Save to file
  png(filename = outfile, width = 3000, height = 2500, res = 350)
  plot_de_heatmap(
    spe = spe,
    de_result = filtered_de_results,
    target_name = target_name,
    logcounts_assay = "logcounts_corrected"
  )
  dev.off()
}

