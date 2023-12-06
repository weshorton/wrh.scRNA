getHeatGenes <- function(data_dt, kind_v, geneCol_v = "Gene", lfcCol_v = "avg_log2FC", lfc_v = 0.25, pvalCol_v = "p_val_adj",
                         pval_v = 0.05, labelAll_v = T) {
  #' Get Heatmap Genes - NOT UPDATED!!
  #' @description filter find marker results to get genes for heatmap plotting
  #' @param data_dt data.table output from FindAllMarkers
  #' @param kind_v either "sig" to grab significant genes or an integer to grab top # genes.
  #' @param geneCol_v name of column containing gene names. Default is Gene
  #' @param lfcCol_v name of column containing log fold change values. Default is avg_log2FC (used to subset top genes, optionally used to determine sig genes)
  #' @param lfc_v log fold change cut-off for significance
  #' @param pvalCol_v name of column containing adjusted pvalues. Default is p_val_adj (used to subset top and sig genes)
  #' @param pval_v pvalue cut-off to use
  #' @param labelAll_v passed to volcanoWrangleMarkers. TRUE - all sig genes are labeled. FALSE - only top 25 sig genes labeled
  #' @param lfcLabelCutOff_v passed to volcanoWrangleMarkers. TRUE - use lfc_v threshold to determine significant genes in addition to pval_v
  #' @details Need to update the label options to mirror volcanowrangleMarkers. Essentially need to add labelGenes_v argument.
  #' @return list with two vectors. First, 'heat' are the genes to be displayed in the heatmap. Second, 'label' contains the genes that are significant. 
  #' If kind_v == "sig", then these will be the same. If kind_v = #, then it will be the significant genes and is used to add a * to their name in heatmap
  #' @export
  
  ### Wrangle with volcano function
  wrangle_dt <- volcanoWrangleMarkers(data_dt, geneCol_v = geneCol_v, lfcCol_v = lfcCol_v, lfc_v = lfc_v, pvalCol_v = pvalCol_v,
                                      pval_v = pval_v, labelAll_v = labelAll_v)
  
  ### Get all significant genes determined by volcano wrangle
  sigGenes_v <- wrangle_dt[diffExp %in% c("UP", "DOWN"),get(geneCol_v)]
  
  ### If kind_v is an integer, grab those
  if (kind_v == "sig") {
    out_ls <- list("heat" = sigGenes_v, "label" = sigGenes_v)
  } else if (length(kind_v == 1) & class(kind_v) == "numeric") {
    topUp_v <- wrangle_dt[get(lfcCol_v) > 0, get(geneCol_v)][1:kind_v]
    topDown_v <- wrangle_dt[get(lfcCol_v) < 0, get(geneCol_v)][1:kind_v]
    out_ls <- list("heat" = c(topUp_v, topDown_v), "label" = sigGenes_v)
  } else {
    stop("kind_v must be either 'sig' or of length 1 and numeric.")
  }
  
  return(out_ls)
  
} # getHeatGenes
