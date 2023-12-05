compareMarkerResults <- function(data1_dt, data2_dt, name_v = NULL, ident1_v, ident2_v, identCol_v, meta_dt = NULL, geneCol_v = "Gene",
                                    runNames_v = c("sct", "rna"), lfc_v = 0.5, pval_v = 0.05, toPlot_v = c("scatter", "volcano"), compareName_v, verbose_v = F) {
  #' Compare Marker Results Log-Fold-Change
  #' @description plot lfc results of different findMarkers runs
  #' @param data1_dt findMarker results
  #' @param data2_dt findmarker results
  #' @param name_v character of global name of what is presented. Usually lymphoid, neoplastic, etc.
  #' @param ident1_v first ident in findmarker call
  #' @param ident2_v second ident in findMarker call
  #' @param identCol_v column name of meta.data that corresponds to idents
  #' @param meta_dt optional meta.data. If not provided, will take cell counts from objects.
  #' @param geneCol_v column name in data1_dt and data2_dt that correspond to genes. Default is "Gene"
  #' @param runNames_v vector of length two specifying the names for data1_dt and data2_dt.
  #' @param lfc_v log fold change to use for cut-offs
  #' @param pval_v adjusted p-value to use for cut-offs
  #' @param toPlot_v vector of plot types to output. Current options are 'scatter' and 'volcano'
  #' @param compareName_v passed to splitVar_v in plotVolcano. Name to describe what's different between data1_dt/data2_dt (i.e. why we're comparing them)
  #' @param verbose_v passed to plotVolcano. Print messages.
  #' @details Take two different FindMarkers() results and compare the log-fold change outputs.
  #' Written to compare runs between different assays (SCT vs RNA), but could compare other things, I think.
  #' @export
  
  data1_dt = currSCTDEG_dt
  data2_dt = currRNADEG_dt
  name_v = currmPop_v
  ident1_v = currTreat_v
  ident2_v = currOtherTreat_v
  identCol_v = "Treatment"
  compareName_v = "assay"
  meta_dt = meta_dt
  geneCol_v = "Gene"
  runNames_v = c("sct", "rna")
  lfc_v = 0.5
  pval_v = 0.05
  toPlot_v = c("scatter", "volcano")
  verbose_v = F
  
  ###
  ### WRANGLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Get number of cells
  cells1_v <- nrow(meta_dt[get(identCol_v) == ident1_v,])
  cells2_v <- nrow(meta_dt[get(identCol_v) == ident2_v,])
  
  ### Base title name
  baseTitle_v <- paste0(paste(cells1_v, ident1_v, "Cells"), " vs ",
                        paste(cells2_v, ident2_v, "Cells"))
  if (!is.null(name_v)) baseTitle_v <- paste(name_v, baseTitle_v, sep = "\n")
  
  ### Update percentage column names
  data_lsdt <- list(data1_dt, data2_dt); names(data_lsdt) <- runNames_v
  data_lsdt <- sapply(data_lsdt, function(x) {
    colnames(x)[colnames(x) == "pct.1"] <- paste0("pct.", ident1_v)
    colnames(x)[colnames(x) == "pct.2"] <- paste0("pct.", ident2_v)
    return(x)
  }, simplify = F)
  
  ### Extract necessary columns
  data_lsdt <- sapply(data_lsdt, function(x) {
    return(x[,mget(grep(paste0(geneCol_v, "|p_val|log2FC"), colnames(x), value = T))])
  }, simplify = F)
  
  ### Merge
  merge_dt <- mergeDTs(data_lsdt = data_lsdt, mergeCol_v = "Gene", suffixes = paste0(".", runNames_v))
  
  ### Split into shared and individual
  cols_lsv <- sapply(c("p_val_adj", "avg_log2FC"), function(x) paste(runNames_v, x, sep = "_"), simplify = F)
  split_lsdt <- list("shared" = merge_dt[!is.na(get(cols_lsv$avg_log2FC[1])) & !is.na(get(cols_lsv$avg_log2FC[2])),],
                     "data1" = merge_dt[is.na(get(cols_lsv$avg_log2FC[2])),],
                     "data2" = merge_dt[is.na(get(cols_lsv$avg_log2FC[1])),])
  
  ### Subset for significant by p-value 
  split_sigP_lsdt <- list("sharedAnd" = split_lsdt$shared[get(cols_lsv$p_val_adj[1]) < pval_v & get(cols_lsv$p_val_adj[2]) < pval_v,],
                          "sharedOr" = split_lsdt$shared[get(cols_lsv$p_val_adj[1]) < pval_v | get(cols_lsv$p_val_adj[2]) < pval_v,],
                          "data1" = split_lsdt$data1[get(cols_lsv$p_val_adj[1]) < pval_v,],
                          "data2" = split_lsdt$data2[get(cols_lsv$p_val_adj[2]) < pval_v,])
  
  ### Subset for higher lfc
  if (!is.null(lfc_v)) {
    split_sigP_filterLFC_lsdt <- list("sharedAnd" = split_sigP_lsdt$sharedAnd[abs(get(cols_lsv$avg_log2FC[1])) > lfc_v & abs(get(cols_lsv$avg_log2FC[2])) > lfc_v,],
                           "sharedOr" = split_sigP_lsdt$sharedOr[abs(get(cols_lsv$avg_log2FC[1])) > lfc_v | abs(get(cols_lsv$avg_log2FC[2])) > lfc_v,],
                           "data1" = split_sigP_lsdt$data1[abs(get(cols_lsv$avg_log2FC[1])) > lfc_v,],
                           "data2" = split_sigP_lsdt$data2[abs(get(cols_lsv$avg_log2FC[2])) > lfc_v,])
  } # fi is.null(lfc_v)
  
  ###
  ### PLOTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  out_lsgg <- list()
  
  ###
  ### SCATTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  if ("scatter" %in% toPlot_v) {
    
    ### Get data
    scatterData_lsdt <- list("allShared" = split_lsdt$shared,
                             "andSharedSigP" = split_sigP_lsdt$sharedAnd,
                             "orSharedSigP" = split_sigP_lsdt$sharedOr)
    
    ### Add others
    if (!is.null(lfc_v)) {
      scatterData_lsdt[["andSharedSigP_filterLFC"]] <- split_sigP_filterLFC_lsdt$sharedAnd
      scatterData_lsdt[["orSharedSigP_filterLFC"]] <- split_sigP_filterLFC_lsdt$sharedOr
    } # fi !is.null(lfc_v)
    
    ### Make plots
    scatter_lsgg <- list()
    
    for (i in 1:length(scatterData_lsdt)) {

      ### Build title
      currTitle_v <- paste(baseTitle_v, paste(runNames_v[1], " vs ", runNames_v[2]), 
                           paste(nrow(scatterData_lsdt[[i]]), names(scatterData_lsdt)[i], "Genes"), sep = "\n")
      
      if (nrow(scatterData_lsdt[[i]]) == 0) {
        if (verbose_v) cat(sprintf("No results for scatter %s.\n", names(scatterData_lsdt)[i]))
        next
      }
      
      ### Make plot
      scatter_lsgg[[names(scatterData_lsdt)[i]]] <- plotXYScatter(data_dt = scatterData_lsdt[[i]],
                                                         x_v = paste0(runNames_v[1], "_avg_log2FC"),
                                                         y_v = paste0(runNames_v[2], "_avg_log2FC"),
                                                         title_v = currTitle_v,
                                                         cor_v = T)
    } # for i
    
    out_lsgg[["scatter"]] <- scatter_lsgg
    
  } # fi scatter_v
  
  ###
  ### VOLCANO ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  if ("volcano" %in% toPlot_v) {
    
    ### Get data
    volcanoData_lsdt <- list("allShared" = split_lsdt$shared)
    volcanoData_lsdt[[runNames_v[1]]] <- split_lsdt$data1
    volcanoData_lsdt[[runNames_v[2]]] <- split_lsdt$data2
    volcanoData_lsdt[["andSharedSigP"]] <- split_sigP_lsdt$sharedAnd
    volcanoData_lsdt[["orSharedSigP"]] <- split_sigP_lsdt$sharedOr
    
    volcanoColorCol_v <- c(rep("diffExp", 3), rep(paste0(compareName_v, "DE"), 2))
    
    ### Add others
    if (!is.null(lfc_v)) {
      volcanoData_lsdt[["andSharedSigP_filterLFC"]] <- split_sigP_filterLFC_lsdt$sharedAnd
      volcanoData_lsdt[["orSharedSigP_filterLFC"]] <- split_sigP_filterLFC_lsdt$sharedOr
      volcanoColorCol_v <- c(volcanoColorCol_v, rep(paste0(compareName_v, "DE"), 2))
    } # fi !is.null(lfc_v)
    
    ### Make plots
    volcano_lsgg <- list()
    
    for (i in 1:length(volcanoData_lsdt)) {
      
      if (nrow(volcanoData_lsdt[[i]]) == 0) {
        if (verbose_v) cat(sprintf("No results for volcano %s.\n", names(volcanoData_lsdt)[i]))
        next
      }
      
      ### Build title
      currTitle_v <- paste(baseTitle_v, paste(runNames_v[1], " vs ", runNames_v[2]), 
                           paste(nrow(volcanoData_lsdt[[i]]), names(volcanoData_lsdt)[i], "Genes"), sep = "\n")
      
      ### Make plot
      volcano_lsgg[[names(volcanoData_lsdt)[i]]] <- plotVolcano(data_dt = volcanoData_lsdt[[i]],
                                                                splitVar_v = compareName_v,
                                                                colorCol_v = volcanoColorCol_v[i],
                                                                title_v = currTitle_v,
                                                                runNames_v = runNames_v,
                                                                geneCol_v = geneCol_v,
                                                                lfc_v = lfc_v,
                                                                pval_v = pval_v,
                                                                ident1_v = ident1_v,
                                                                verbose_v = verbose_v)
    } # for i
    
    out_lsgg[["volcano"]] <- volcano_lsgg
    
  } # fi volcano
  
  return(out_lsgg)
  
} # compareMarkerResults




