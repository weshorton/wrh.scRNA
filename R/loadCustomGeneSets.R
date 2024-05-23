loadCustomGeneSets <- function(file_v, keep_v = NULL, omit_v = NULL) {
  #' Load Custom Gene Sets
  #' @description
  #' Load custom gene sets
  #' @param file_v path to excel file containing one or more sheets with gene lists. See details.
  #' @param keep_v vector of sheet names to keep
  #' @param omit_v vector of sheet names to omit.
  #' @details
  #' Only one of keep_v or omit_v may be set.
  #' Gene list format is:
  #' 1. one sheet per gene set
  #' 2. sheet name is gene set name
  #' 3. Must at least contain a column GENE with gene names.
  #' @return list of tables
  #' @export
  
  ### Read in gene sets and edit if necessary
  if (!is.null(keep_v) & !is.null(omit_v)) stop("Only one of keep_v or omit_v may be set.\n")
  geneSets_lsdt <- wrh.rUtils::readAllExcel(file_v)
  if (!is.null(keep_v)) geneSets_lsdt <- geneSets_lsdt[which(names(geneSets_lsdt) %in% keep_v)]
  if (!is.null(omit_v)) geneSets_lsdt <- geneSets_lsdt[which(!(names(geneSets_lsdt) %in% omit_v))]
  
  ### Change gene names
  ### So far, any gene set can be coerced to mouse gene format using the simpleCap()
  ### function, except for the MHC genes.
  geneSets_lsdt <- sapply(names(geneSets_lsdt), function(xx) {
    
    x <- geneSets_lsdt[[xx]]
    
    if ("Gene" %in% colnames(x)) {
      
      return(x)
      
    } else if ("GENE" %in% colnames(x)) {
      
      if (xx == "MHC") { 
        x$Gene <- x$GENE
      } else { 
        x$Gene <- convertMouseHumanGenes(x$GENE, species_v = "hg")
      } 
      
      return(x)
      
    } else {
      stop("Must have column 'Gene' or 'GENE'.")
    }
  })
  
  ### Return
  return(geneSets_lsdt)
  
} # loadCustomGeneSets