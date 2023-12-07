wrangleCollapsePops <- function(obj, name_v, batch_v, col_v = "sPop", newCol_v = NULL) {
  #' Wrangle Collapse Populations
  #' @description
    #' Use population map tables to collapse original sPop names
    #' in seurat object meta.data
  #' @param obj seurat object
  #' @param name_v name of seurat object. Usually b12_lymphoid or similar. Used to get map_dt!
  #' @param batch_v name of batch. Used if name_v doesn't have batch in the name
  #' @param col_v what column of meta.data are we updating? sPop is only one so far.
  #' @param newCol_v name of column to put collapsed populations in. Default NULL will replace values in col_v
  #' @return Either: (1) same exact seurat object, but col_v is now updated with collapsed populations
  #' (2) same seurat object with one more metadata column containing the collapsed populations
  #' @export
  
  # Get name of map table
  if (grepl("batch", name_v)) {
    mapName_v <- paste0(gsub("atch", "", name_v), "PopMap_dt")
  } else {
    mapName_v <- paste0(gsub("atch", "", batch_v), "_", name_v, "PopMap_dt")
  }
  
  # Only run if map table exists
  if (!exists(mapName_v)) {
    
    cat(sprintf("Map table %s doesn't exist. This object won't be changed.
                Fine if this is neo or other, not good if lymphoid or myeloid.\n", mapName_v))
    return(obj)
    
  } else {
    
    # Wrangle map
    map_dt <- eval(as.name(mapName_v))
    collapse_v <- unique(map_dt$collapsePop)
    updateCol_v <- ifelse(is.null(newCol_v), col_v, newCol_v)
    
    # Run update
    obj@meta.data[[updateCol_v]] <- as.character(obj@meta.data[[col_v]])
    for (c_v in collapse_v) {
      obj@meta.data[obj@meta.data[[col_v]] %in% map_dt[collapsePop == c_v,sPop], updateCol_v] <- unique(map_dt[collapsePop == c_v,collapsePop])
    }
    obj@meta.data[[updateCol_v]] <- factor(obj@meta.data[[updateCol_v]], levels = unique(obj@meta.data[[updateCol_v]]))
    
    # Return
    return(obj)
    
  } # fi !exists(mapName_v)
  
} # wrangleCollapsePops