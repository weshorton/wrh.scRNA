wrangleCollapsePops <- function(obj, name_v, batch_v, col_v = "sPop") {
  #' Wrangle Collapse Populations
  #' @description
    #' Use population map tables to collapse original sPop names
    #' in seurat object meta.data
  #' @param obj seurat object
  #' @param name_v name of seurat object. Usually b12_lymphoid or similar. Used to get map_dt!
  #' @param batch_v name of batch. Used if name_v doesn't have batch in the name
  #' @param col_v what column of meta.data are we updating? sPop is only one so far.
  #' @return same exact seurat object, but col_v is now updated with fewer populations
  #' @export
  
  if (grepl("batch", name_v)) {
    mapName_v <- paste0(gsub("atch", "", name_v), "PopMap_dt")
  } else {
    mapName_v <- paste0(gsub("atch", "", batch_v), "_", name_v, "PopMap_dt")
  }
  map_dt <- eval(as.name(mapName_v))
  
  collapse_v <- unique(map_dt$collapsePop)
  
  obj@meta.data[[col_v]] <- as.character(obj@meta.data[[col_v]])
  
  for (c_v in collapse_v) {
    obj@meta.data[obj@meta.data[[col_v]] %in% map_dt[collapsePop == c_v,sPop], col_v] <- unique(map_dt[collapsePop == c_v,collapsePop])
  }
  
  obj@meta.data[[col_v]] <- factor(obj@meta.data[[col_v]], levels = unique(obj@meta.data[[col_v]]))
  
  return(obj)
  
} # wrangleCollapsePops