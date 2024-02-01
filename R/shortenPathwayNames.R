shortenPathwayNames <- function(gsea_res) {
  #' Shorten Pathway Names
  #' @description shorten gsea pathway names to make output easier to read
  #' @param gsea_res input gsea results
  #' @return gsea_res, except 'pathway' column is shorter
  #' @export
  
  for (i in 1:nrow(gsea_res)) {
    
    ### Get pathway name into words
    currPathWords_v <- strsplit(gsea_res$pathway[i], split = "_")[[1]]
    
    ### Shorten if longer than 4 words
    if (length(currPathWords_v) < 4) {
      
      currPath_v <- paste(currPathWords_v, collapse = "_")
      
    } else {
      
      currFirst_v <- currPathWords_v[1:3]
      currSecond_v <- currPathWords_v[4:length(currPathWords_v)]
      currSecond_v <- sapply(currSecond_v, function(x) strtrim(x, 4))
      currPath_v <- paste(c(currFirst_v, currSecond_v), collapse = "_")
      
    } # fi
    
    gsea_res$pathway[i] <- currPath_v
    
  } # for i
  
  return(gsea_res)
  
} # shortenPathwayNames

