mutationCounter <- function(mutation_dt, idCols_v = c("Hugo_Symbol", "Tumor_Sample_Barcode"), 
                            specificSummary_lsv = list("4x" = c("I4", "I5", "I6")),
                            specificCol_v = "Tumor_Sample_Barcode") {
  #' Mutation Counter
  #' @description summarize WES mutation table in multiple ways
  #' @param mutation_dt WES mutation table
  #' @param idCols_v column names of mutation_dt to use for mutation summary tables. Should generally be a gene ID and a sample ID
  #' @param specificSummary_lsv list of vectors of sample IDs to use for a specific mutation count
  #' @param specificCol_v column name of mutation_dt that is used for specific summary
  #' @details Summary columns
  #' 'nSamp' simple count of how many samples gene is mutated in
  #' 'nSampSpecific' simple count of how many specific samples gene is mutated in
  #' 'total' total count of all mutations in all samples of this gene
  #' @return summary data.table with above columns, sort by decreasing nSamp
  #' @export
  
  # Create a count summary table
  counts_df <- as.data.frame.matrix(table(mutation_dt[,mget(idCols_v)]))
  
  # Get binary count and total count for full table
  nSamp_v <- apply(counts_df, 1, function(x) length(which(x != 0)))
  total_v <- rowSums(counts_df)
  
  # Get binary counts and total counts for specific samples
  if (!is.null(specificSummary_lsv)) {
    
    specific_nSamp_lsv <- specific_total_lsv <- list()
    for (i in 1:length(specificSummary_lsv)) {
      currCount_v <- apply(counts_df[,colnames(counts_df) %in% specificSummary_lsv[[i]]], 1, function(x) length(which(x != 0)))
      currTotal_v <- rowSums(counts_df[,colnames(counts_df) %in% specificSummary_lsv[[i]]])
      specific_nSamp_lsv[[names(specificSummary_lsv)[i]]] <- currCount_v
      specific_total_lsv[[names(specificSummary_lsv[i])[i]]] <- currTotal_v
    } # for i
    
  } # fi !is.null(specificSummary_lsv)
  
  # Add full dataset summary
  counts_df$nSamp <- nSamp_v
  counts_df$total <- total_v
  
  # Add specific summary
  if (!is.null(specificSummary_lsv)) {
    
    for (i in 1:length(specificSummary_lsv)) {
      name_v <- names(specificSummary_lsv)[i]
      sCol_v <- paste0("nSamp", name_v)
      tCol_v <- gsub("nSamp", "total", sCol_v)
      counts_df[,sCol_v] <- specific_nSamp_lsv[[name_v]]
      counts_df[,tCol_v] <- specific_total_lsv[[name_v]]
    } # for i
    
  } # fi !is.null(specificSummary_lsv)
  
  # Convert to data.table and outpu8t
  out_dt <- convertDFT(counts_df, newName_v = idCols_v[1])
  out_dt <- out_dt[,lapply(.SD, function(x) if (class(x) == "integer") {return(as.numeric(x))}else{return(x)})]
  setorder(out_dt, -nSamp)
  return(out_dt)
  
} # mutationCounter