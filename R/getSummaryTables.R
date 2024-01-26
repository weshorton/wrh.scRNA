getSummaryTables <- function(obj, subset_v = NULL, subCol_v = "mPop",
                             summary_lsv = list("treat" = "Treatment",
                                                "pop" = "mPop",
                                                "combo" = c("Treatment", "mPop"))) {
  #' Get Tables
  #' @description Get summary tables of selected columns
  #' @param obj seurat object to summarize. Must have specified columns
  #' @param subset_v subset data on provided value
  #' @param subCol_v specify which column to subset using subset_v
  #' @param summary_lsv list of vectors of columns to summarize.
  #' @description Currently make 3 different tables: treatment, population, and treatment x population
  #' Treatment should almost always be 'treat = Treatment'; population should be 'mPop' if subset_v = NULL and
  #' 'sPop' if not (or sPop2 or another version); treatment x population must agree with the others.
  #' @export
  
  ### Extract meta data
  meta_dt <- as.data.table(obj@meta.data)
  
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
      #treat <- rbind(treat, list("Total", sum(treat$nCells)))
      out_ls[["treat"]] <- treat
    } # fi
    
    ### Major or minor populations
    if ("pop" %in% names(summary_lsv)) {
      pop <- as.data.table(table(as.character(meta_dt[[summary_lsv$pop]]))); colnames(pop) <- c("Pop", "nCells")
      out_ls[["pop"]] <- pop
    }
    
    ### Combo
    if ("combo" %in% names(summary_lsv)) {
      #newName_v <- ifelse(is.null(subset_v), paste(subset_v, "Pop", collapse = ", "))
      newName_v <- "Pop"
      combo_tab <- t(as.data.frame.matrix(table(meta_dt[,mget(summary_lsv$combo)])))
      combo_tab <- combo_tab[unname(apply(combo_tab, 1, function(x) length(which(x == 0)) != length(x))),]
      combo <- convertDFT(combo_tab, newName_v = newName_v)
      out_ls[["combo"]] <- combo
    }
    
  } else {
    
    ### Treatments
    if ("treat" %in% names(summary_lsv)) {
      treat <- convertDFT(as.data.frame.matrix(table(meta_dt[,mget(c(summary_lsv$treat, "batchID"))])), newName_v = "Treat")
      treat[,Total := rowSums(.SD), .SDcols = batches_v]
      out_ls[["treat"]] <- treat
    }
    
    ### Major or minor populations
    if ("pop" %in% names(summary_lsv)) {
      pop <- convertDFT(as.data.frame.matrix(table(meta_dt[,mget(c(summary_lsv$pop, "batchID"))])), newName_v = "Pop")
      pop[,Total := rowSums(.SD), .SDcols = batches_v]
      out_ls[["pop"]] <- pop
    }
    
    ### Combo
    if ("combo" %in% names(summary_lsv)) {
      combo_ls <- sapply(batches_v, function(x) {
        convertDFT(t(as.data.frame.matrix(table(meta_dt[batchID == x, mget(summary_lsv$combo)]))), newName_v = "Pop")
      }, simplify = F)
      merge_dt <- convertDFT(t(as.data.frame.matrix(table(meta_dt[,mget(summary_lsv$combo)]))), newName_v = "Pop")
      combo_ls[["Total"]] <- merge_dt
      out_ls[["combo"]] <- combo_ls
    }
    
  } # if (length(batches))
  
  return(out_ls)
  
} # getSummaryTables

