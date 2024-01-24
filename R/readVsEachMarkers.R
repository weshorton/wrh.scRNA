readVsEachMarkers <- function(inDir_v, dirs_v, grep_v = "2x|3xNR|3xR|4x|d80|Ent|NT|PTX|3x|Unt", toRead_v = NULL,
                              toReadSpop_v = NULL) {
  #' Read VS Each Markers
  #' @description
  #' Read in directory of vs each find markers results. For each meta population tested (e.g. lymphoid, myeloid, neoplastic),
  #' there is a sub-directory for each subpopulation. For each of those, there is a FindMarkers excel workbook for each treatment. 
  #' Each workbook has a sheet for every treatment it's compared to.
  #' Want to read these in so we have a list of list of lists of all of these results
  #' @param inDir_v path to the main data directory that contains dirs_v
  #' @param dirs_v vector of directories to read data from. Should be batch3_lymphoid, batch3_myeloid, batch3_neoplastic (or batch12)
  #' @param grep_v character string to extract names. Should be "2x|3xNR|3xR|4x|d80|Ent|NT|PTX"
  #' @param toRead_v used to grep which files to read in within an spop directory
  #' @param toReadSpop_v used to grep which spop directories to read in
  #' @return list of length(dirs_v), each element is a list whose length is the number of treatments run, and each of those
  #' is a list whose length is the number of treatments compared to base treatment. Each element of the list is the data.table
  #' of differential expression results.
  #' @export
  
  ### Grab all of the subDirs in each main directory
  subDirs_lsv <- sapply(dirs_v, function(x) {
    y <- list.files(file.path(inDir_v, x))
    if (!is.null(toReadSpop_v)) y <- intersect(y, toReadSpop_v)
    return(y)}, simplify = F, USE.NAMES = T)
  
  ### Read in each workbook within each of those
  deg_lslslslsdt <- sapply(dirs_v, function(x) {
    subDirs_v <- subDirs_lsv[[x]]
    deg_lslslsdt <- sapply(subDirs_v, function(y) {
      files_v <- list.files(file.path(inDir_v, x, y))
      if (!is.null(toRead_v)) files_v <- files_v[grep(toRead_v, files_v)]
      data_lsdt <- sapply(files_v, function(z) {
        readAllExcel(file.path(inDir_v, x, y, z))
      }, simplify = F, USE.NAMES = T)
      names(data_lsdt) <- sapply(names(data_lsdt), function(xx) grep(grep_v, unlist(strsplit(xx, split = "_")), value = T))
      return(data_lsdt)
    }, simplify = F, USE.NAMES = T)
  }, simplify = F, USE.NAMES = T)
  
  ### Output
  return(deg_lslslslsdt)
  
} # readVsEachMarkers