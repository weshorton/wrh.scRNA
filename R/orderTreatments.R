orderTreatments <- function(names_lsv, order_v, skip_v = "batch|lymphoid|myeloid|neoplastic") {
  #' Order treatments
  #' @description
  #' Order a vector of treatments either standalone or as part of a list.
  #' @param names_lsv list of vectors or just a single vector
  #' @param order_v vector of matching treatments in the correct order
  #' @param skip_v regex expression of list elements that should be skipped b/c they're not treatments.
  #' @return same object, but ordered
  #' @export
  
  ### Handle class
  if (!is.logical(all.equal(class(names_lsv), "list"))) names_lsv <- list(names_lsv)
  
  ### Order
  out_lsv <- sapply(names_lsv, function(x) {
    if (FALSE %in% grepl(skip_v, x)) {
      missing_v <- setdiff(x, order_v)
      x <- as.character(sort(factor(x, levels = order_v)))
      if (length(missing_v) > 0) { 
        warning(sprintf("Input names had values not found in treatments.
                      Appending to front. Please check! Treatments: %s\n",
                        paste(missing_v, collapse = "; ")))
        x <- c(missing_v, x)
      } # fi length()
    } # fi treat
    return(x)}, simplify = F, USE.NAMES = T)
  
  ### Return
  return(out_lsv)
  
} # orderTreatments
