prepProjectR <- function(obj, assay_v = "SCT", slot_v = "data", filter_dt = NULL, filterCol_v = "gene_name") {
  #' Prep data for ProjectR Run
  #' @description Exctract correct assay/slot counts and filter
  #' @param obj seurat object
  #' @param assay_v name of assay to get data from
  #' @param filter_dt optional table containing at least one column of gene names that will be checked against genes in obj
  #' @param filterCol_v column name from filter_dt used to identify genes.
  #' @return a data.frame of counts, filtered for provided genes.
  #' @export
  
  ### Check size of scale.data
  if (slot_v == "scale.data") {
    dim_v <- dim(obj@assays[[assay_v]][slot_v])
    dim2_v <- dim(obj@assays[[assay_v]]$counts)
    if (!is.logical(all.equal(dim_v, dim2_v))) {
      warning(sprintf("sale.data slot provided, but full object wasn't return in SCT run.
              Scale.data only has %s genes, while full object has %s.
              ProjectR will only be run on the scale.data genes, so any signature genes missing will be omitted.
              Be sure to run SCTransform with return.only.var.genes = F if you want to use them all.\n"))
    } # fi dim_v == dim2_v
  } # fi slot_v == "scale.data
  
  ### Extract data
  data_mat <- obj@assays[[assay_v]][slot_v]
  
  ### Filter
  if (!is.null(filter_dt)) {
    data_mat <- as.matrix(data_mat[rownames(data_mat) %in% filter_dt[[filterCol_v]],])
    missingGenes_v <- setdiff(filter_dt[[filterCol_v]], rownames(data_mat))
    if (length(missingGenes_v) > 0) cat(sprintf("\n%s signature genes not found.\n", length(missingGenes_v)))
  } # fi
  
  ### Output
  return(data_mat)
} # prepProjectR
