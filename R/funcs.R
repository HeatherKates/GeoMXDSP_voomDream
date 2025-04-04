library(SpatialExperiment)
library(ggplot2)
library(gridExtra)

# Define a function to run drawPCA with up to two subset variables and a color variable
run_drawPCA <- function(spe, assay, color, subset_var1 = NULL, subset_val1 = NULL, subset_var2 = NULL, subset_val2 = NULL) {
  # Subset the SpatialExperiment object if subset variables are provided
  if (!is.null(subset_var1)) {
    spe <- spe[, colData(spe)[[subset_var1]] == subset_val1]
  }
  if (!is.null(subset_var2)) {
    spe <- spe[, colData(spe)[[subset_var2]] == subset_val2]
  }
  
  # Run drawPCA
  drawPCA(spe, assay = assay, color = color)
}
