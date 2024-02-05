wrangleCollapsePops <- function(obj, name_v, batch_v, col_v = "sPop", newCol_v = NULL, refCol_v = "sPop") {
  #' Wrangle Collapse Populations
  #' @description
    #' Use population map tables to collapse original sPop names
    #' in seurat object meta.data
  #' @param obj seurat object
  #' @param name_v name of seurat object. Usually b12_lymphoid or similar. Used to get map_dt!
  #' @param batch_v name of batch. Used if name_v doesn't have batch in the name
  #' @param col_v what column of meta.data are we updating? sPop is only one so far.
  #' @param newCol_v name of column to put collapsed populations in. Default NULL will replace values in col_v
  #' @param refCol_v name of column to use as reference in map_dt. Should always be sPop I think.
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
    
    # Check map
    inMapNotInData_v <- setdiff(map_dt[[col_v]], obj@meta.data[[col_v]])
    if (length(inMapNotInData_v) > 0) {
      if (col_v == "collapsePop") {
        cat(sprintf("At least one population in map_dt is not in object.
                    Update column (col_v) is 'collapsePop', so likely issue is wrangleCollapsePops has been run
                    previously and col_v is shortened now (no spaces or slashes). Going to apply same shortening to 
                    map_dt and then map."))
        map_dt <- as.data.table(apply(map_dt, 2, function(x) gsub(" ", "\\.", gsub("\\/", "-", x))))
        collapse_v <- unique(map_dt$collapsePop)
      } else if (col_v == "sPop") {
        warning("At least one population in map_dt is not in the object, but the update column is sPop.
                map_dt was specifically made with sPop and this should not happen. Be sure to check output.")
      } else {
        stop(sprintf("col_v argument may be 'sPop' or 'collapsePop'. %s was provided.\n", col_v))
      }
    } # fi 
    
    # Run update
    obj@meta.data[[updateCol_v]] <- as.character(obj@meta.data[[col_v]])
    for (c_v in collapse_v) {
      obj@meta.data[obj@meta.data[[col_v]] %in% map_dt[collapsePop == c_v,get(refCol_v)], updateCol_v] <- unique(map_dt[collapsePop == c_v,collapsePop])
    }
    
    # Shorten
    obj@meta.data[[col_v]] <- gsub(" ", "\\.", gsub("\\/", "-", obj@meta.data[[col_v]]))
    obj@meta.data[[updateCol_v]] <- gsub(" ", "\\.", gsub("\\/", "-", obj@meta.data[[updateCol_v]]))
    
    # Factor
    obj@meta.data[[col_v]] <- factor(obj@meta.data[[col_v]], levels = unique(obj@meta.data[[col_v]]))
    obj@meta.data[[updateCol_v]] <- factor(obj@meta.data[[updateCol_v]], levels = unique(obj@meta.data[[updateCol_v]]))
    
    # Return
    return(obj)
    
  } # fi !exists(mapName_v)
  
} # wrangleCollapsePops