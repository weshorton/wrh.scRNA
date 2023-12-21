getProjectRPlottingData <- function(obj, batch_v, pop_v, sigName_v = "sncScore", binaryName_v = "SnC", popCol_v = "collapsePop", metaCol_v = "Treatment") {
  #' Get ProjectR Plotting Data
  #' @description extract meta.data from seurat object for easier plotting with ggplot
  #' @param obj seurat object with new columns created by prepProjectRPlotting
  #' @param batch_v name of batch (required to get appropriate treatments)
  #' @param pop_v name of population (required to get appropriate population order)
  #' @param sigName_v must be same as what was used in prepProjectRPlotting
  #' @param binaryName_v must be same as what was used in prepProjectRPlotting
  #' @param popCol_v name of meta.data columbn containing population information. Can usually be 'mPop', 'sPop', or 'collapsedPop'
  #' @param metaCol_v name of the meta.data column containing grouping variable. Can be more than one, Treatment is required.
  #' @return data.frame
  #' @export
  
  # Extract data
  cols_v <- c(popCol_v, metaCol_v, grep(paste(sigName_v, binaryName_v, sep = "|"), colnames(obj@meta.data), value = T))
  plot_df <- obj@meta.data[,cols_v]
  
  # Factorize treatments and populations
  if (!"Treatment" %in% metaCol_v) stop("Must have 'Treatment' in getCols_v.\n")
  if (batch_v == "batch3") {
    plot_df$Treatment <- factor(plot_df$Treatment, levels = rev(b3_treats_v))
  } else if (batch_v == "batch12") {
    plot_df$Treatment <- factor(plot_df$Treatment, levels = rev(b12_treats_v))
  } else if (batch_v == "orig") {
    plot_df$Treatment <- factor(plot_df$Treatment, levels = rev(b12_treats_v))
  }
  
  # Factor spops if not (lymphoid and myeloid should be) Warning! I'm allowing length(popCol_v) > 1, which works above b/c just using c()
  # to get colnames. If neoplastic, length(popCol_v) must equal 1, but it should always do so b/c we didn't ID the pops.
  if (pop_v == "neoplastic") {
    plot_df[[popCol_v]] <- factor(plot_df[[popCol_v]], 
                                  levels = unique(plot_df[[popCol_v]])[order(as.numeric(gsub("neo.c", "", unique(plot_df[[popCol_v]]))))])
  }
  
  # Return
  return(plot_df)
  
} # getProjectRPlottingData
