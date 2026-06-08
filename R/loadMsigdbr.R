loadMsigdbr <- function(species_v = "MM", category_v = c("H", "C2", "C2", "C2", "C5", "C5"), 
                        subcategory_v = c("", "CP:KEGG_LEGACY", "CP:BIOCARTA", "CP:REACTOME", "GO:BP", ""), 
                        name_v = c("HALLMARK", "KEGG", "BIOCARTA", "REACTOME", "GOBP", "GO")) {
  #' Load msigdbr
  #' @description
    #' Load various reference sets from msigdbr database
  #' @param species_v what species to use. 
  #' For mouse, use "MM". For human, use "HS"
  #' @param category_v vector of categories to grab data for
  #' @param subcategory_v vector of subcategories. Must be same length as category_v. Use "" if no subcat
  #' @param name_v character vector of names to label outputs.
  #' @return list of lists of msigdbr results
  #' @export
  
  ### Handle species
  species_v <- sapply(species_v, function(x) {
    if (x %in% c("Mouse", "mouse", "mmu", "Mmu", "Mus musculus")) {
      x <- "Mus musculus"
    } else if (x %in% c("Human", "human", "hg", "Hg", "Homo sapiens")) {
      x <- "Homo sapiens"
    } # fi x
    return(x)})
  if (length(species_v) == 1) species_v <- rep(species_v, length(category_v))
  
  ### Handle db species
  db_species_v <- sapply(species_v, function(x) {
    if (x %in% c("Mouse", "mouse", "mmu", "Mmu", "Mus musculus", "MM")) {
      x <- "MM"
    } else if (x %in% c("Human", "human", "hg", "Hg", "Homo sapiens")) {
      x <- "HS"
    } # fi x
    return(x)})
  if (length(db_species_v) == 1) db_species_v <- rep(db_species_v, length(category_v))
  
  ### Run table
  run_dt <- data.table("DB_species" = db_species_v, "Species" = species_v, "Cat" = category_v, "SubCat" = subcategory_v, "Name" = name_v)
  run_dt[run_dt == ""] <- NA
  
  ### Fix collections if mouse
  run_dt[DB_species == "MM", Cat := gsub("^H$", "MH", Cat)]
  run_dt[DB_species == "MM", Cat := gsub("^C", "M", Cat)]
  
  ### Download
  out_ls <- sapply(1:nrow(run_dt), function(x) {
    ### Get sub-category from run table
    if (is.na(run_dt[x,SubCat])) {subcat_v <- NULL} else {subcat_v <- run_dt[x,SubCat]}
    # New in msigdbr 10.0.0
    y <- msigdbr::msigdbr(db_species = run_dt[x,DB_species], species = run_dt[x,Species], collection = run_dt[x,Cat], subcollection = subcat_v)
    ### format output
    out <- split(x = y$gene_symbol, f = y$gs_name)
    return(out)}, simplify = F)
  names(out_ls) <- run_dt$Name
  
  ### Output
  return(out_ls)
  
} # loadMsigdbr
