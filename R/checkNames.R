checkNames <- function(obj, deg) {
  #' Check Object and DEG Names
  #' @description Make sure that seurat object and DEG table have same pop names
  #' @param obj seurat object
  #' @param deg deg table
  #' @export
  
  ### Get object names
  if ("collapsePop" %in% colnames(obj@meta.data)) {
    # objNames_v <- union(unique(obj$sPop), unique(obj$collapsePop))
    objNames_v <- unique(obj$collapsePop)
  } else {
    objNames_v <- unique(obj$sPop)
  } # fi
  
  ### Get deg names
  degNames_v <- names(deg)
  
  ### Run Checks
  diff_ls <- list("objOnly" = setdiff(objNames_v, degNames_v),
                  "degOnly" = setdiff(degNames_v, objNames_v))
  
  ### Update
  if (length(diff_ls$degOnly) > 0) {
    degOnly_v <- paste(diff_ls$degOnly, collapse = ", ")
    warning(sprintf("Population(s) with DEG results but not in object:\t\t%s\nThis shouldn't happen. Check your inputs!\n", degOnly_v))
  }
  
  if (length(diff_ls$objOnly) > 0) {
    objOnly_v <- paste(diff_ls$objOnly, collapse = ", ")
    warning(sprintf("Population(s) in object, but no DEG results :\t\t%s\nThis shouldn't happen. Check your inputs!\n", objOnly_v))
  }
  
  return(diff_ls)
} # checkNames