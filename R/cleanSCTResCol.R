cleanSCTResCol <- function(obj, grep_v = "SCT_snn.*", which_v = 2) {
  #' Clean SCT Res Col
  #' @description remove extra resolution columns in obj@meta.data
  #' @param obj seurat object with multiple clustering resolutions
  #' @param grep_v regex used to grab the columns to be cleaned.
  #' @param which_v which of the results returned by grep_v should be kept?
  #' @details When processing the data, there are two runs of SCT with and without doublets.
  #' We want the second one, without doublets and here will remove the first. This fxn can be
  #' used with other parameter values, but the defaults are what are needed for this project.
  #' @return same input seurat object, but only one resolution column
  #' @export
  
  ### Get columns
  grepRes_v <- grep(grep_v, colnames(obj@meta.data), value = T)
  
  ### If needed, select 1 and remove the rest
  if (length(grepRes_v) > 1) {
    toRm_v <- grepRes_v[-which_v]
    for (col_v in toRm_v) obj[[col_v]] <- NULL
    grepRes_v <- grepRes_v[which_v]
  } # fi
  
  ### Output
  return(obj)
} # cleanSCTResCol