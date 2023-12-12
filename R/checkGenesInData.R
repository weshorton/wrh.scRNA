checkGenesInData <- function(obj, objName_v, geneList_v, geneListName_v = "", verbose_v = T) {
  #' Check Genes In Data
  #' @description See if provided genes are present in the dataset
  #' @param obj seurat object to check
  #' @param objName_v name of seurat object for update
  #' @param geneList_v vector of gene names to check. OR named list of vectors of gene names
  #' @param geneListName_v name of gene list for output. Only used if geneList_v is a vector.
  #' @param verbose_v logical. true - print updates. false - no printing
  #' @return same gene list, but subset for only those that are in obj.
  #' @export
  
  ### Get intersect
  if (is.logical(all.equal(class(geneList_v), "list"))) {
    isList_v <- T
    intersect_v <- sapply(geneList_v, function(x) intersect(x, rownames(obj)), simplify = F, USE.NAMES = T)
  } else {
    isList_v <- F
    intersect_v <- intersect(geneList_v, rownames(obj))
  }
  
  ### Update
  if (verbose_v) {
    if (isList_v) {
      invisible(sapply(names(intersect_v), function(x) cat(sprintf("Found %s of the %s %s module gene list in %s (%s%%)\n",
                                                                   length(intersect_v[[x]]), length(geneList_v[[x]]), x,
                                                                   objName_v, round(length(intersect_v[[x]])/length(geneList_v[[x]])*100,digits = 2)))))
    } else {
      cat(sprintf("Found %s of the %s %s module gene list in %s (%s%%)\n",
                  length(intersect_v), length(geneList_v), geneListName_v,
                  objName_v, round(length(intersect_v)/length(geneList_v)*100,digits = 2)))
    } # fi isList_v
  } # fi verbose
  
  ### Remove genes
  return(intersect_v)
  
} # checkGenesInData

