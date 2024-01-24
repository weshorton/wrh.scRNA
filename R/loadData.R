loadData <- function(markerDir_v = NULL, markerClass_v = "global", metaDir_v = NULL, rdsDir_v = NULL,
                     batch_v, pop_v, toRead_v = "_2x_|_3xNR_|_3xR_|_4x_|_d80_|_Ent_|_NT_|_PTX_|_3x_|_3xR\\+4x|_Unt_") {
  #' Load Data
  #' @description
    #' Load DEG, meta, and seurat data
  #' @param markerDir_v path to where marker results are held
  #' @param markerClass_v either 'global' or 'vsEach', indicating which read in function to use.
  #' @param metaDir_v path to metadata directory
  #' @param rdsDir_v path to seurat objects
  #' @param batch_v batch name that identifies files to read in
  #' @param pop_v population names that identify files to read in.
  #' @param toRead_v used to grep what to input.
  #' @return list with marker results, metadata, and seurat objects. Read in objects with provided directories. If dir == NULL, not read in.
  #' @export
  
  out_ls <- list()
  
  ### DEG Results
  if (!is.null(markerDir_v)) {
    
    if (markerClass_v == "global") {
      fxn <- wrh.scRNA::readGlobalMarkers
    } else if (markerClass_v == "vsEach") {
      fxn <- wrh.scRNA::readVsEachMarkers
    } else {
      stop(sprintf("markerClass_v must be 'global' or 'vsEach'. %s provided.\n", markerClass_v))
    } # fi markerClass
    
    if (!is.null(toRead_v)) {
      out_ls[["deg"]] <- fxn(inDir_v = markerDir_v, dirs_v = grep(paste(pop_v, collapse = "|"), list.files(markerDir_v, pattern = batch_v), value = T), toRead_v = toRead_v)
    } else {
      out_ls[["deg"]] <- fxn(inDir_v = markerDir_v, dirs_v = grep(paste(pop_v, collapse = "|"), list.files(markerDir_v, pattern = batch_v), value = T))
    }
    
  } # fi !is.null(markerDir_v)
  
  ### Metadata
  if (!is.null(metaDir_v)) {
    
    out_ls[["meta"]] <- fread(file.path(metaDir_v, paste0(batch_v, "_meta.txt")))
    
  } # fi !is.null(metaDir_v)
  
  ### Seurat objects
  if (!is.null(rdsDir_v)) {
    
    files_v <- grep(paste(paste(pop_v, collapse = "|"), paste(simpleCap(pop_v), collapse = "|"), sep = "|"), list.files(rdsDir_v, pattern = batch_v), value = T)
    out_ls[["seurat"]] <- sapply(files_v, function(x) readRDS(file.path(rdsDir_v, x)), USE.NAMES = T)
    names(out_ls[["seurat"]]) <- gsub("\\.rds", "", files_v)
  } # fi !is.null(rdsDir_v)
  
  ### Output
  return(out_ls)
  
} # loadData