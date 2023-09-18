volcanoWrangleMarkers <- function(data_dt, geneCol_v = "Gene", 
                                     lfcCol_v = "avg_log2FC", lfc_v = 0.25, 
                                     pvalCol_v = "p_val_adj", pval_v = 0.05,
                                     labelAll_v = T, lfcLabelCutOff_v = T) {
  #' Wrangle Marker Data for Volcano Plot
  #' @description Add information for volcano plot
  #' @param data_dt data.table output from FindAllMarkers
  #' @param geneCol_v name of column containing gene names. Default is Gene
  #' @param lfcCol_v name of column containing log fold change values. Default is avg_log2FC
  #' @param lfc_v log fold change cut-off for labels
  #' @param pvalCol_v name of column contianing adjusted pvalues. Default is p_val_adj
  #' @param pval_v pvalue cut off for labels
  #' @param labelAll_v logical. TRUE - label all significant genes. FALSE - label top 25 only
  #' @param lfcLabelCutOff_v logical. TRUE - only label genes that are significant and above the lfc_v threshold.
  #' @return original data_dt with updated metadata info for plotting
  #' @export
  
  ### Make a new column called "diffExp" that binarizes direction of change
  data_dt$diffExp <- "NO"
  data_dt[get(lfcCol_v) > lfc_v & get(pvalCol_v) < pval_v, diffExp := "UP"]
  data_dt[get(lfcCol_v) < -lfc_v & get(pvalCol_v) < pval_v, diffExp := "DOWN"]
  
  ### Pull out new tables of only up/down and sort them by pvalue
  up_dt <- data_dt[diffExp == "UP",][order(get(pvalCol_v))]
  dn_dt <- data_dt[diffExp == "DOWN",][order(get(pvalCol_v))]
  
  ### Extract genes for labels
  if (labelAll_v) {
    upGenes_v <- up_dt[[geneCol_v]]
    dnGenes_v <- dn_dt[[geneCol_v]]
  } else {
    upGenes_v <- up_dt[1:min(up_dt[,.N], 25), get(geneCol_v)]
    dnGenes_v <- up_dt[1:min(dn_dt[,.N], 25), get(geneCol_v)]
  } # fi
  
  ### Add new column for labels. If in upGenes_v or dnGenes_v, will get gene name
  ### if not, will be blank ('')
  data_dt$DElabel <- ""
  data_dt[get(geneCol_v) %in% upGenes_v, DElabel := get(geneCol_v)]
  data_dt[get(geneCol_v) %in% dnGenes_v, DElabel := get(geneCol_v)]
  
  ### Make diffExp facctor
  data_dt$diffExp <- factor(data_dt$diffExp, levels = c("NO", "DOWN", "UP"))
  
  ### Return
  return(data_dt)
  
} # volcanoWrangleAllMarkers