checkForGene <- function(deg_lslslslsdt, gene_v, pval_v = NULL) {
  #' Check For Gene
  #' @description
    #' Check all of the deg tables in either the vsEach or global outputs for a given gene
  #' @param deg_lslslslsdt either a vsEach or a global DEG output.
  #' @param gene_v gene to check
  #' @param pval_v optional numeric value to subset the deg tables for significant hits of gene_v only
  #' @return data.table with one row per hit of gene_v in deg_lslslslsdt
  #' @export
  #' 
  
  sigResults_mat <- NULL
  
  for (i in 1:length(deg_lslslslsdt)) {
    
    currI_v <- names(deg_lslslslsdt)[i]
    deg_lslslsdt <- deg_lslslslsdt[[currI_v]]
    if (length(deg_lslslsdt) == 0) next
    
    for (j in 1:length(deg_lslslsdt)) {
      
      currJ_v <- names(deg_lslslsdt)[j]
      deg_lslsdt <- deg_lslslsdt[[currJ_v]]
      if (length(deg_lslsdt) == 0) next
      
      for (k in 1:length(deg_lslsdt)) {
        
        currK_v <- names(deg_lslsdt)[k]
        deg_lsdt <- deg_lslsdt[[currK_v]]
        if (is.null(deg_lsdt)) next
        
        if (is.data.table(deg_lsdt)) {
          
          padjCol_v <- ifelse("p_val_adj" %in% colnames(deg_lsdt), "p_val_adj", ifelse("padj" %in% colnames(deg_lsdt), "padj", stop("Must have 'padj' or 'p_val_adj' in data.table")))
          pvalCol_v <- ifelse("p_val" %in% colnames(deg_lsdt), "p_val", ifelse("pval" %in% colnames(deg_lsdt), "pval", stop("Must have 'pval' or 'p_val' in data.table")))
          
          if (!is.null(pval_v)) {
            subDEG_dt <- deg_lsdt[Gene == gene_v & get(padjCol_v) < pval_v,]
          } else {
            subDEG_dt <- deg_lsdt[Gene == gene_v,]
          }
          if (nrow(subDEG_dt) > 0) {
            row_v <- c(currI_v, currJ_v, currK_v, subDEG_dt$avg_log2FC, subDEG_dt[[pvalCol_v]], subDEG_dt[[padjCol_v]])
            sigResults_mat <- rbind(sigResults_mat, row_v)
          } # fi
          
        } else {
          
          for (l in 1:length(deg_lsdt)) {
            
            currL_v <- names(deg_lsdt)[l]
            deg_dt <- deg_lsdt[[currL_v]]
            
            padjCol_v <- ifelse("p_val_adj" %in% colnames(deg_dt), "p_val_adj", ifelse("padj" %in% colnames(deg_dt), "padj", stop("Must have 'padj' or 'p_val_adj' in data.table")))
            pvalCol_v <- ifelse("p_val" %in% colnames(deg_dt), "p_val", ifelse("pval" %in% colnames(deg_dt), "pval", stop("Must have 'pval' or 'p_val' in data.table")))
            
            if (is.null(pval_v)) {
              subDEG_dt <- deg_dt[Gene == gene_v]
            } else {
              subDEG_dt <- deg_dt[Gene == gene_v & get(padjCol_v) < pval_v]
            }
            if (nrow(subDEG_dt) > 0) {
              row_v <- c(currI_v, currJ_v, currK_v, currL_v, subDEG_dt$avg_log2FC, subDEG_dt[[pvalCol_v]], subDEG_dt[[padjCol_v]])
              sigResults_mat <- rbind(sigResults_mat, row_v)
            } # fi
            
          } # for l
        } # fi
      } # for k
    } # for j
  } # for i
  
  if (!is.null(sigResults_mat)) {
    if (ncol(sigResults_mat) == 7) {
      colnames(sigResults_mat) <- c("Object", "Cluster", "Treat1", "Treat2", "l2fc", "pval", "padj")
    } else {
      colnames(sigResults_mat) <- c("Object", "Treat1", "Treat2", "l2fc", "pval", "padj")
    }
    sigResults_dt <- as.data.table(sigResults_mat)
    for (col_v in c("l2fc", "pval", "padj")) sigResults_dt[, (col_v) := as.numeric(get(col_v))]
  } else {
    sigResults_dt <- NULL
  }
  
  return(sigResults_dt)
  
}