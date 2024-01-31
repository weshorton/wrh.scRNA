avgConservedMarkers <- function(conserved_dt, geneCol_v = "Gene", extraCols_v = c("gene_biotype", "description")) {
  #' Average Conserved Markers
  #' @description take the output of FindConservedMarkers and average the log2FC and adj p-value of all the comparisons
  #' @param conserved_dt output of FindConservedMarkers, converted to data.table.
  #' @param geneCol_v column name of gene names
  #' @param extraCols_v any columns that aren't in default FCM() output. Default here is gene_biotype which is added when I merge output with annotations.
  #' @return data.table with same number of rows as conserved_dt, but only 'Gene', 'avg_padj', 'avg_log2FC'
  #' @export
  
  ### Get columns
  pctCols_v <- grep("_pct", colnames(conserved_dt), value = T)
  minMaxP_v <- c("max_pval", "minimum_p_val")
  pvalCols_v <- setdiff(grep("_p_val$", colnames(conserved_dt), value = T), minMaxP_v)
  padjCols_v <- grep("_adj", colnames(conserved_dt), value = T)
  lfcCols_v <- grep("_log2FC", colnames(conserved_dt), value = T)
  
  ### Calculate averages
  avgPadj_v <- rowSums(convertDFT(conserved_dt[,mget(c(geneCol_v, padjCols_v))], col_v = geneCol_v))
  avgPadj_dt <- data.table("V1" = names(avgPadj_v), "avg_padj" = avgPadj_v); colnames(avgPadj_dt)[1] <- geneCol_v
  
  avgLFC_v <- rowSums(convertDFT(conserved_dt[,mget(c(geneCol_v, lfcCols_v))], col_v = geneCol_v))
  avgLFC_dt <- data.table("V1" = names(avgLFC_v), "avg_log2FC" = avgLFC_v); colnames(avgLFC_dt)[1] <- geneCol_v
  
  ### Output
  out_dt <- merge(avgPadj_dt, avgLFC_dt, by = geneCol_v, sort = F)
  return(out_dt)
  
} # avgConservedMarkers

