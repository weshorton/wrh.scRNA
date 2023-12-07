makePairs <- function(vector_v, pairs_lsv = NULL, counter_v = 1) {
  #' Make Pairs
  #' @description
    #' Get all unique pairs of values from a vector.
    #' (order doesn't matter)
  #' @param vector_v vector of elements to combine pairwise
  #' @param paris_lsv list of pairs fed back in recursively
  #' @param counter_v counter for recursion
  #' @return list of pairs. Each element of vector_v will be a list element whose values are the other elements of vector_v.
  #' @export
  
  if (is.null(pairs_lsv)) pairs_lsv <- list()
  if (counter_v != length(vector_v)) {
    pairs_lsv[[vector_v[counter_v]]] <- vector_v[(counter_v+1):length(vector_v)]
    makePairs(vector_v = vector_v, pairs_lsv = pairs_lsv, counter_v = counter_v + 1)
  } else {
    return(pairs_lsv)
  } # fi
} # makePairs

