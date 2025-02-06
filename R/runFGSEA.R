### GSEA wrapper
runFGSEA <- function(data_df, geneCol_v = "Gene", rankCol_v = "avg_log2FC", pathways_v, minSize_v = 1, 
                     maxSize_v = nrow(data_df), seed_v = 42) {
  #' FGSEA Wrapper
  #' @description Run fgsea and wrangle output
  #' @param data_df table of data to run. must have log2FC column and column for genes
  #' @param geneCol_v name of gene column. Default is Gene
  #' @param rankCol_v column name of data_df to use to rank. Default is avg_log2FC
  #' @param pathways_v list of gene pathways to test against
  #' @param minSize_v passed to fgsea. minimum size of gene set to test
  #' @param maxSize_v passed to fgsea. Maximum size of gene set to test
  #' @param seed_v set for reproducibility
  #' @return data.frame with gsea results
  #' @export
  
  ### Make sure sorted by rank col
  setorderv(data_df, rankCol_v, order = -1)
  
  ### Grab ranks and name them with genes
  rank <- data_df[[rankCol_v]]
  names(rank) <- data_df[[geneCol_v]]
  
  ### Run function
  set.seed(seed_v)
  fgsea_res <- fgsea::fgsea(pathways = pathways_v,
                     stats = rank,
                     minSize = minSize_v,
                     maxSize = 500)
  
  ### Wrangle output
  results <- fgsea::collapsePathways(fgseaRes = fgsea_res, pathways = pathways_v, stats = rank)
  df <- as.data.frame(results$parentPathways)
  df <- tibble::rownames_to_column(df, var="pathway")
  names(df)[2] <- 'Parent'
  df <- merge(fgsea_res, df, by='pathway', sort=FALSE)
  df <- (df[order(-NES), ])
  #df <- subset(df, is.na(df$Parent))
  
  ### Return
  return(df)
}
