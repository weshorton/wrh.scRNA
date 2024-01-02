getTreatmentComparisons <- function(treats_v) {
  #' Get Treatment Comparisons
  #' @description
  #' Create a list of pairwise treatment comparisons from an input vector of treatments
  #' @param treats_v vector of treatments to run comparisons for
  #' @description Returns a non-repeating list of pairs, so order of treats_v matters.
  #' @return list. One element per "treat1" comparison. Each "treat1" element has a vector
  #' containing each "treat2" to be compared against it.
  #' @export
  
  treatComparisons_lsv <- list()
  for (i in 1:length(treats_v)) {
    if (i == length(treats_v)) next
    treatComparisons_lsv[[treats_v[i]]] <- treats_v[(i+1):length(treats_v)]
  } # for i
  
  return(treatComparisons_lsv)
  
}