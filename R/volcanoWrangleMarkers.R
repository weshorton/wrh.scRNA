volcanoWrangleMarkers <- function(data_dt, geneCol_v = "Gene", 
                                  lfcCol_v = "avg_log2FC", lfc_v = 0.25, 
                                  pvalCol_v = "p_val_adj", pval_v = 0.05,
                                  labelGenes_v = NULL, labelAll_v = T, labelTop_v = 20, labelDir_v = "both"
                                  ) {
  #' Wrangle Marker Data for Volcano Plot
  #' @description Add information for volcano plot
  #' @param data_dt data.table output from FindMarkers
  #' @param geneCol_v name of column containing gene names. Default is Gene
  #' @param lfcCol_v name of column containing log fold change values. Default is avg_log2FC
  #' @param lfc_v log fold change cut-off for labels
  #' @param pvalCol_v name of column contianing adjusted pvalues. Default is p_val_adj
  #' @param pval_v p-value cut off for labels
  #' @param labelGenes_v character vector of gene names to label. See details.
  #' @param labelAll_v logical. TRUE - label all significant genes. See details.
  #' @param labelTop_v numeric. if labelAll_v == F, number of significant genes to label
  #' @details
    #' Reformat the default FindMarkers() DEG output to contain information we want for volcano plots.
    #' 1. Binarize differential expression by assigning all significant genes (i.e. p-value less than pval_v) to either "UP"
    #' or "DOWN" depending of if lfcCol_v is greater than or less than lfc_v, respectively.
    #' 2. Assign labels to specific genes in the plot (3 options). Default is option b.
    #'  a. If labelGenes_v is set (i.e. not NULL), then these genes will be labeled, if present. This takes precedence 
    #'  b. If labelAll_v == T (and is.null(labelGenes_v)), then label all genes that pass the pvalue and lfc threshold
    #'  c. If labelAll_v == F (and is.null(labelGenes_v)), then label labelTop_v genes. If labelTop_v == 0, no labels.
  #' @return original data_dt with updated metadata info for plotting
  #' @export
  
  ### Set local variables
  dirVarCheck_lsv <- list("up" = c("up", "Up", "u", "U"), "down" = c("down", "Down", "D", "d"), "both" = c("both", "Both", "B", "b"))
  
  ### Make a new column called "diffExp" that binarizes direction of change
  data_dt$diffExp <- "NO"
  data_dt[get(lfcCol_v) > lfc_v & get(pvalCol_v) < pval_v, diffExp := "UP"]
  data_dt[get(lfcCol_v) < -lfc_v & get(pvalCol_v) < pval_v, diffExp := "DOWN"]
  
  ### Pull out new tables of only up/down and sort them by pvalue
  up_dt <- data_dt[diffExp == "UP",][order(get(pvalCol_v))]
  dn_dt <- data_dt[diffExp == "DOWN",][order(get(pvalCol_v))]
  
  ### Assign labels
  if (labelAll_v) {
    
    upGenes_v <- up_dt[[geneCol_v]]
    dnGenes_v <- dn_dt[[geneCol_v]]
    labelUpGenes_v <- labelDnGenes_v <- labelNonGenes_v <- NULL
    
  } else {
    
    upGenes_v <- dnGenes_v <- NULL
    labelUpGenes_v <- labelDnGenes_v <- labelNonGenes_v <- NULL
    
    if (!is.null(labelTop_v)) {
      upGenes_v <- c(upGenes_v, up_dt[1:min(up_dt[,.N], labelTop_v), get(geneCol_v)])
      dnGenes_v <- c(dnGenes_v, dn_dt[1:min(dn_dt[,.N], labelTop_v), get(geneCol_v)])
    } # fi labelTop_v
    
    if (!is.null(labelGenes_v)) {
      labelUpGenes_v <- intersect(up_dt[[geneCol_v]], labelGenes_v)
      labelDnGenes_v <- intersect(dn_dt[[geneCol_v]], labelGenes_v)
      labelNonGenes_v <- intersect(data_dt[diffExp == "NO",get(geneCol_v)], labelGenes_v)
      
      ### Remove from label top
      if (!is.null(labelTop_v)) {
        upGenes_v <- setdiff(upGenes_v, labelUpGenes_v)
        dnGenes_v <- setdiff(dnGenes_v, labelDnGenes_v)
      } # fi labelTop_v
    } # fi labelGenes_v
    
  } # fi labelAll_v
  
  ### Now have to figure out labels.
  data_dt$DElabel <- ""
  
  if (is.null(upGenes_v) & is.null(labelUpGenes_v)) {
    
    warning("Both upGenes_v and labelUpGenes_v are NULL, indicating:
    
            labelAll_v == F
            labelTop_v == NULL
            labelGenes_v == NULL\n\nat least one of these needs to be set.")
    
    return(data_dt)
    
  } else {
    
    if (!is.null(upGenes_v)) {
      
      if (labelDir_v %in% dirVarCheck_lsv$up | labelDir_v %in% dirVarCheck_lsv$both) {
        data_dt[get(geneCol_v) %in% upGenes_v, DElabel := get(geneCol_v)]
      } # fi up
      
      if (labelDir_v %in% dirVarCheck_lsv$dn | labelDir_v %in% dirVarCheck_lsv$both) {
        data_dt[get(geneCol_v) %in% dnGenes_v, DElabel := get(geneCol_v)]
      } # fi dn
      
      
    } # fi !is.null(upGenes_v)
    
    if (!is.null(labelUpGenes_v)) {
      
      data_dt$setLabel <- ""
      
      if (labelDir_v %in% dirVarCheck_lsv$up | labelDir_v %in% dirVarCheck_lsv$both) {
        data_dt[get(geneCol_v) %in% labelUpGenes_v, setLabel := get(geneCol_v)]
      } # fi up
      
      if (labelDir_v %in% dirVarCheck_lsv$dn | labelDir_v %in% dirVarCheck_lsv$both) {
        data_dt[get(geneCol_v) %in% labelDnGenes_v, setLabel := get(geneCol_v)]
      } # fi dn
      
      data_dt[get(geneCol_v) %in% labelNonGenes_v, setLabel := get(geneCol_v)]
      
    } # fi !is.null(labelUpGenes_v)
    
  } # fi double NULL
  
  ### Make diffExp factor
  data_dt$diffExp <- factor(data_dt$diffExp, levels = c("NO", "DOWN", "UP"))
  if (!is.null(labelUpGenes_v)) data_dt$setColor <- data_dt$diffExp
  
  ### Return
  return(data_dt)
  
} # volcanoWrangleAllMarkers