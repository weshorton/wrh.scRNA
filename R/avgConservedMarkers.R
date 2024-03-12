avgConservedMarkers <- function(conserved_dt, geneCol_v = "Gene", extraCols_v = c("gene_biotype", "description")) {
  #' Average Conserved Markers
  #' @description take the output of FindConservedMarkers and average the log2FC and adj p-value of all the comparisons
  #' @param conserved_dt output of FindConservedMarkers, converted to data.table.
  #' @param geneCol_v column name of gene names
  #' @param extraCols_v any columns that aren't in default FCM() output. Default here is gene_biotype which is added when I merge output with annotations.
  #' @return data.table with same number of rows as conserved_dt, but only 'Gene', 'max_pval', 'minimump_p_val', 'avg_padj' (calculated in script), and 'avg_log2FC' (calculated in script)
  #' @export
  
  ### Get columns
  pctCols_v <- grep("_pct", colnames(conserved_dt), value = T)
  minMaxP_v <- c("max_pval", "minimump_p_val")
  pvalCols_v <- setdiff(grep("_p_val$", colnames(conserved_dt), value = T), minMaxP_v)
  padjCols_v <- grep("_adj", colnames(conserved_dt), value = T)
  lfcCols_v <- grep("_log2FC", colnames(conserved_dt), value = T)
  
  ###
  ### Calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Calculate min and max p-values (to compare to theirs)
  myMinP_v <- apply(conserved_dt[,mget(pvalCols_v)], 1, min)
  myMaxP_v <- apply(conserved_dt[,mget(pvalCols_v)], 1, max)
  myMinPadj_v <- apply(conserved_dt[,mget(padjCols_v)], 1, min)
  myMaxPadj_v <- apply(conserved_dt[,mget(padjCols_v)], 1, max)
  
  ### Calculate averages
  temp_df <- convertDFT(conserved_dt[,mget(c(geneCol_v, padjCols_v))], col_v = geneCol_v)
  avgPadj_v <- rowMeans(convertDFT(conserved_dt[,mget(c(geneCol_v, padjCols_v))], col_v = geneCol_v))
  avgPadj_dt <- data.table("V1" = names(avgPadj_v), "avg_padj" = avgPadj_v); colnames(avgPadj_dt)[1] <- geneCol_v
  
  avgLFC_v <- rowMeans(convertDFT(conserved_dt[,mget(c(geneCol_v, lfcCols_v))], col_v = geneCol_v))
  avgLFC_dt <- data.table("V1" = names(avgLFC_v), "avg_log2FC" = avgLFC_v); colnames(avgLFC_dt)[1] <- geneCol_v
  
  ###
  ### Combine Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Combine my calculations, min/max, and averages.
  myCalcs_dt <- data.table("V1" = conserved_dt[[geneCol_v]], 
                           "myMax_p" = myMaxP_v, "myMin_p" = myMinP_v,
                           "myMax_padj" = myMaxPadj_v, "MyMin_padj" = myMinPadj_v)
  colnames(myCalcs_dt)[1] <- geneCol_v
  minMax_dt <- conserved_dt[,mget(c(geneCol_v, minMaxP_v))]
  avg_dt <- merge(avgPadj_dt, avgLFC_dt, by = geneCol_v, sort = F)
  
  ### Combine them all
  out_dt <- merge(minMax_dt, myCalcs_dt, by = geneCol_v, sort = F)
  out_dt <- merge(out_dt, avg_dt, by = geneCol_v, sort = F)
  return(out_dt)
  
} # avgConservedMarkers

