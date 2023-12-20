prepProjectRPlotting <- function(sig_lssig, obj, slot_v = "data", pvals_v = c(0.01, 0.05), sigName_v = "sncScore", binaryName_v = "SnC") {
  #' Prep ProjectR Plotting
  #' @description Add projection output to seurat object metadata
  #' @param sig_lssig projection list object output by projectR (contains 'pval' and 'projection', 1-row arrays of same number of columns)
  #' @param obj seurat object used to generate projections
  #' @param slot_v slot used in projectR. Used for naming columns
  #' @param pval_v numeric vector of p-values to filter by
  #' @param sigName_v name to call the signature. Used for naming columns
  #' @param binaryName_v character vector used to name cells that pass the filter. Cells will be 'binaryName_v' and 'not binaryName_v'
  #' @return updated seurat object with new columns.
  #' ProjectR results:
  #'  1. sig_lssig$pval will be an exact copy, called paste(slot_v, sigName_v, "pval", sep = "_")
  #'  2. sig_lssig$projection will be an exact copy, called paste(slot_v, sigName_v, "weight", sep = "_")
  #' Wrangled ProjectR results (for each pvalue provided)
  #' 1. column called paste(slot_v, sigName_v, "weight_sig", pval_v, sep = "_") is a copy of the projection weights, but only for significant ones
  #' 2. column called paste(slot_v, pval_v, "sig", binaryName_v, sep = "_") binarizes the weights into "pos" and "neg" if they are significant
  #' 3. column called paste(slot_v, pval_v, binaryName_v, sep = "_") binarizes ALL the weights into "pos" and "neg"
  #' For each pvalue provided in pval_v, two columns will
  #' @export
  
  ### Make column names for new results
  weightCol_v <- paste(slot_v, sigName_v, "weight", sep = "_")
  pvalCol_v <- paste(slot_v, sigName_v, "pval", sep = "_")
  
  ### Add signature score info to metadata
  obj <- AddMetaData(object = obj, metadata = t(sig_lssig$projection), col.name = weightCol_v)
  obj <- AddMetaData(object = obj, metadata = t(sig_lssig$pval), col.name = pvalCol_v)
  
  ### Add another set of columns
  ### One set will hold significant weights only and binarize from that
  ### The other will be binary of all weights
  for (i in 1:length(pvals_v)) {
    
    # New colum name
    pval_v <- pvals_v[i]
    sigWeightCol_v <- paste(slot_v, sigName_v, "weight_sig", pval_v, sep = "_")

    # Fill only with weights of significant ones
    obj[[sigWeightCol_v]] <- ifelse((obj@meta.data[[pvalCol_v]] < pval_v), obj@meta.data[[weightCol_v]], NA)
    
    # Two new binary columns (one for all genes, and one for significant)
    sigBinaryCol_v <- paste(slot_v, pval_v, "sig", binaryName_v, sep = "_")
    binaryCol_v <- paste(slot_v, pval_v, binaryName_v, sep = "_")
    
    # Add first value
    obj[[sigBinaryCol_v]] <- paste0("pos_", binaryName_v)
    obj[[binaryCol_v]] <- paste0("pos_", binaryName_v)
    
    # Add second value (negative)
    obj@meta.data[(!is.na(obj[[sigWeightCol_v]]) & weightCol_v < 0), (sigBinaryCol_v)] <- paste0("neg_", binaryName_v)
    obj@meta.data[(weightCol_v < 0), (binaryCol_v)] <- paste0("neg_", binaryName_v)
    
    # Add third value (non-scoring)
    obj@meta.data[is.na(obj[[sigWeightCol_v]]), (binaryCol_v)] <- paste0("not_", binaryName_v)
    obj@meta.data[(weightCol_v == 0), (sigBinaryCol_v)] <- paste0("not_", binaryName_v)
    
  } # for i in pvals
  
  ### Return object
  return(obj)
  
} # prepProjectRPlotting
