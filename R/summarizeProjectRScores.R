summarizeProjectRScores <- function(plot_df, cols_v = NULL) {
  #' Summarize ProjectR Scores
  #' @description Get a summary of scores
  #' @param plot_df plot table output gy getProjectRPlottingData
  #' @param cols_v columns to summarize by
  #' @return summary table
  #' @export
  
  # Get all the weights and tabulate
  weights_v <- grep("weight", colnames(plot_df), value = T)
  weightCounts_lstab <- sapply(weights_v, function(x) tabulateScores(plot_df = plot_df, weight_v = x, cols_v = cols_v), simplify = F)
  
  # Handle differently if cols
  if (is.null(cols_v)) {
    
    # Check
    if (length(unique(sapply(weightCounts_lstab, sum))) > 1) stop("Unequal tabulation of weights...\n")
    
    # Wrangle
    summary_dt <- do.call(rbind, weightCounts_lstab)
    rownames(summary_dt) <- paste0("weight", gsub("^.*weight", "", rownames(summary_dt)))
    summary_lsdt <- list("summary" = wrh.rUtils::convertDFT(summary_dt, newName_v = "Projection"))
    
  } else {
    
    # Check
    if (length(unique(sapply(weightCounts_lstab, function(x) sum(x[,mget(setdiff(colnames(x), cols_v))], na.rm = T)))) > 1) stop("Unequal tabulation of weights...\n")
    summary_lsdt <- sapply(weightCounts_lstab, function(x) {
      if (!is.logical(all.equal(class(x), c("data.table", "data.frame")))) {
        return(convertDFT(x, newName_v = "Rownames"), simplify = F, USE.NAMES = T)
      } else {
        return(x)
      } 
    }, simplify = F, USE.NAMES = T)
  } # fi
  
  # Build Pos/NegSummary
  summary_lsdt <- sapply(summary_lsdt, function(x) {
    x[,pctNeg := round(neg/(neg+pos)*100, digits = 2)]
    x[,pctPos := round(pos/(neg+pos)*100, digits = 2)]
    x[,pctNotSig := round(zero/(neg+pos+zero)*100, digits = 2)]
    return(x)}, simplify = F, USE.NAMES = T)
  
  # Return
  return(summary_lsdt)
  
} # summarizeProjectRScores


tabulateScores <- function(plot_df, weight_v, cols_v = NULL) {
  #' Tabulate Signature Scores
  #' @description
    #' Summarize input table into Pos/Neg/Zero by provided variables
  #' @param plot_df data.frame containing plotting information
  #' @param weight_v column name of weight that is getting summarized
  #' @param cols_v other columns in plot_df that should be used to stratify the summary of weight_v
  #' @return table with counts of each group 
  #' @export
  
  ### Convert to data.table
  if (!is.logical(all.equal(class(plot_df), c("data.table", "data.frame")))) {
    plot_df <- convertDFT(plot_df, newName_v = "Rownames")
  } # fi data.table
  
  ### Subset data
  subData_lsdf <- list("pos" = plot_df[get(weight_v) > 0 & !is.na(plot_df[[weight_v]]),],
                       "neg" = plot_df[get(weight_v) < 0 & !is.na(plot_df[[weight_v]]),],
                       "zero" = plot_df[get(weight_v) == 0 | is.na(plot_df[[weight_v]]),])
  
  ### Tabulate
  counts_v <- sapply(subData_lsdf, function(x) {
    y <- x[, .(count = .N), by = cols_v]
    return(y)
  }, simplify = F, USE.NAMES = T)
  
  if (is.null(cols_v)) {
    counts_v <- unlist(counts_v)
    names(counts_v) <- gsub("\\.count", "", names(counts_v))
  } else {
    counts_v <- wrh.rUtils::mergeDTs(data_lsdt = counts_v, mergeCol_v = cols_v)
  }
  
  ### Return
  return(counts_v)
  
} # tabulateScores