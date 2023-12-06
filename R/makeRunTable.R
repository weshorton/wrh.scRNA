makeRunTable <- function(list_ls, treats_v, ...) {
  #' Make run Table
  #' @description
  #' Use names from input list to determine what should be run.
  #' @param list_ls input list
  #' @param treats_v treatment vector, used for ordering
  #' @param ... passed to orderTreatments
  #' @return data.table with one row per unique combination of columns
  #' @export
  
  names_lsv <- orderTreatments(names_lsv = wrh.rUtils::getAllListNames(list_ls),
                               order_v = treats_v, ...)
  
  runTable_dt <- as.data.table(expand.grid(names_lsv, stringsAsFactors = F))
  setkeyv(runTable_dt, cols = paste0("Var", 1:ncol(runTable_dt)))
  
  counts_v <- apply(runTable_dt, 1, function(x) {
    y <- table(as.character(x))
    if (any(y>1)) {return(2)} else {return(1)}
  })
  
  runTable_dt <- runTable_dt[which(counts_v < 2),]
  
  return(runTable_dt)
  
} # makeRunTable