leadingEdgeVolcano <- function(data_dt, gsea_res, pathway_v, treat_v, lfc_v = 0.5, pval_v = 0.05, title_v = NULL) {
  #' Leading Edge Volcano
  #' @description
  #' Make a volcano plot of Leading Edge genes from GSEA results
  #' @param data_dt table of DEG. Must be same one that created gsea_res
  #' @param gsea_res output of GSEA runFGSEA (must be wrangled...change this)
  #' @param pathway_v vector of the pathway name to make plot for
  #' @details
  #' For a given pathway in gsea_res, grab the leading edge genes and then plot
  #' a volcano plot using the supplied DEG results.
  #' @return ggplot
  #' @export 
  
  # Wrangle data, if needed
  if (!"diffExp" %in% colnames(data_dt)) {
    data_dt <- wrh.scRNA::volcanoWrangleMarkers(data_dt = data_dt, lfc_v = lfc_v, pval_v = pval_v)
  }
  
  # Update labels to only be pathway genes
  leGenes_v <- gsea_res[pathway == pathway_v, leadingEdge][[1]]
  data_dt$DElabel <- ""
  data_dt[Gene %in% leGenes_v, DElabel := Gene]
  
  # Update title
  if (!is.null(title_v)) title_v <- paste0(title_v, "\n", pathway_v, " Genes")
  
  # Make plot
  plot_gg <- plotVolcano(data_dt = data_dt, splitVar_v = NULL, runNames_v = '', geneCol_v = "Gene", lfc_v = lfc_v, pval_v = pval_v,
                         colorCol_v = "diffExp", ident1_v = treat_v, title_v = title_v)
  
  
} # leadingEdgeVolcano