readGlobalMarkers <- function(inDir_v, dirs_v, grep_v = "2x|3xNR|3xR|4x|d80|Ent|NT|PTX") {
  #' Read Global Markers
  #' @description
  #' Read in directory of global find markers results. For each meta population tested (e.g. lymphoid, myeloid, neoplastic),
  #' there is a FindMarkers excel workbook for each treatment. Each workbook has a sheet for every treatment it's compared to.
  #' Want to read these in so we have a list of list of lists of all of these results
  #' @param inDir_v path to the main data directory that contains dirs_v
  #' @param dirs_v vector of directories to read data from. Should be batch3_lymphoid, batch3_myeloid, batch3_neoplastic (or batch12)
  #' @param grep_v character string to extract names. Should be "2x|3xNR|3xR|4x|d80|Ent|NT|PTX"
  #' @return list of length(dirs_v), each element is a list whose length is the number of treatments run, and each of those
  #' is a list whose length is the number of treatments compared to base treatment. Each element of the list is the data.table
  #' of differential expression results.
  #' @export
  
  ### Grab all of the workbooks in each directory
  subFiles_lsv <- sapply(dirs_v, function(x) list.files(file.path(inDir_v, x)), simplify = F, USE.NAMES = T)
  
  ### Read in each workbook
  deg_lslslsdt <- sapply(dirs_v, function(x) {
    subFiles_v <- subFiles_lsv[[x]]
    deg_lslsdt <- sapply(subFiles_v, function(y) readAllExcel(file.path(inDir_v, x, y)), simplify = F, USE.NAMES = T)
    names(deg_lslsdt) <- sapply(names(deg_lslsdt), function(y) grep(grep_v, unlist(strsplit(y, split = "_")), value = T))
    return(deg_lslsdt)
  }, simplify = F, USE.NAMES = T)
  
  ### Output
  return(deg_lslslsdt)
  
} # readGlobalMarkers