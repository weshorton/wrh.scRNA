getSummaryTables <- function(obj, subset_v = NULL, subCol_v = "mPop", treatCol_v = "Treatment", popCol_v = "mPop",
                             summary_lsv = list("treat" = "Treatment",
                                                "pop" = "mPop",
                                                "combo" = c("Treatment", "mPop"))) {
  #' Get Tables of Cell Counts
  #' @description Get summary tables of selected columns
  #' @param obj seurat object to summarize. Must have the specified columns
  #' @param subset_v subset data by matching provided value with entries in subCol_v (e.g. myeloid, neoplastic
  #' @param subCol_v specify which column to subset using subset_v
  #' @param treatCol_v column specifying where treatment is recorded in meta.data
  #' @param popCol_v column name for population annotations
  #' @param summary_lsv list of vectors of columns to summarize.
  #' @details Currently make 3 different tables: treatment, population, and treatment x population.  
  #' Treatment should almost always be 'treat = Treatment'.  
  #' Population should be pop='mPop' if subset_v = NULL and sPop' if not (or sPop2 or another version).  
  #' treatment x population (combo) must agree with the others.  
  #' 
  #' @export
  
  ### Extract meta data
  meta_dt <- as.data.table(obj@meta.data)
  
  ### Get treatment and population order
  if (is(obj@meta.data[[popCol_v]], "character")) {
    warning("Assigning population level order. Should only happen for neoplastic.\n")
    levels_v <- unique(obj@meta.data[[popCol_v]])
    levels_v <- levels_v[order(as.numeric(gsub("neo\\.c", "", levels_v)))]
    obj@meta.data[[subCol_v]] <- factor(obj@meta.data[[popCol_v]], levels = levels_v)
  } else {
    levels_v <- levels(obj@meta.data[[popCol_v]])
  } # fi
  
  treatLevels_v <- intersect(names(treatColors_v), unique(meta_dt[[treatCol_v]]))
  
  ### Subset if desired
  if (!is.null(subset_v)) meta_dt <- meta_dt[get(subCol_v) %in% subset_v,]
  
  ### Output list
  out_ls <- list()
  
  ### Have to do it slightly differently if b3 or b12
  batches_v <- unique(meta_dt$batchID)
  
  if (length(batches_v) == 1) {
  
    ### Treatments
    if ("treat" %in% names(summary_lsv)) {
      treat <- as.data.table(table(as.character(meta_dt[[summary_lsv$treat]]))); colnames(treat) <- c("Treat", "nCells")
    } # fi treat
    
    ### Major or minor populations
    if ("pop" %in% names(summary_lsv)) {
      pop <- as.data.table(table(as.character(meta_dt[[summary_lsv$pop]]))); colnames(pop) <- c("Pop", "nCells")
    } # fi pop
    
    ### Combo
    if ("combo" %in% names(summary_lsv)) {
      combo_tab <- t(as.data.frame.matrix(table(meta_dt[,mget(summary_lsv$combo)])))
      combo_tab <- combo_tab[unname(apply(combo_tab, 1, function(x) length(which(x == 0)) != length(x))),]
      combo <- convertDFT(combo_tab, newName_v = "Pop")
      combo$Pop <- factor(combo$Pop, levels = levels_v); setorder(combo, Pop)
      combo <- combo[,mget(c("Pop", treatLevels_v))]
      out_ls[["combo"]] <- combo
    } # fi combo
    
  } else {
    
    ### Treatments
    if ("treat" %in% names(summary_lsv)) {
      treat <- convertDFT(as.data.frame.matrix(table(meta_dt[,mget(c(summary_lsv$treat, "batchID"))])), newName_v = "Treat")
      treat[,nCells := rowSums(.SD), .SDcols = batches_v]
    } # fi treat
    
    ### Major or minor populations
    if ("pop" %in% names(summary_lsv)) {
      pop <- convertDFT(as.data.frame.matrix(table(meta_dt[,mget(c(summary_lsv$pop, "batchID"))])), newName_v = "Pop")
      pop[,nCells := rowSums(.SD), .SDcols = batches_v]
    } # fi pop
    
    ### Combo
    if ("combo" %in% names(summary_lsv)) {
      combo_ls <- sapply(batches_v, function(x) {
        y <- convertDFT(t(as.data.frame.matrix(table(meta_dt[batchID == x, mget(summary_lsv$combo)]))), newName_v = "Pop")
        y$Pop <- factor(y$Pop, levels = levels_v); setorder(y, Pop)
        y <- y[,mget(c("Pop", intersect(treatLevels_v, colnames(y))))]
      }, simplify = F)
      merge_dt <- convertDFT(t(as.data.frame.matrix(table(meta_dt[,mget(summary_lsv$combo)]))), newName_v = "Pop")
      merge_dt$Pop <- factor(merge_dt$Pop, levels = levels_v); setorder(merge_dt, Pop)
      merge_dt <- merge_dt[,mget(c("Pop", intersect(treatLevels_v, colnames(merge_dt))))]
      combo_ls[["Total"]] <- merge_dt
      out_ls[["combo"]] <- combo_ls
    } # fi combo
    
  } # if (length(batches))
  
  ### Set Orders and add to output
  if ("treat" %in% names(summary_lsv)) {
    treat$Treat <- factor(treat$Treat, levels = treatLevels_v)
    setorder(treat, Treat)
    out_ls[["treat"]] <- treat
  } # fi pop
  
  if ("pop" %in% names(summary_lsv)) {
    pop$Pop <- factor(pop$Pop, levels = levels_v)
    setorder(pop, Pop)
    out_ls[["pop"]] <- pop
  } # fi pop
  
  return(out_ls)
  
} # getSummaryTables

