plotVolcano <- function(data_dt, splitVar_v = NULL, runNames_v = '', geneCol_v, lfc_v, pval_v, ident1_v, colorCol_v, title_v = NULL, verbose_v = F,
                        labelGenes_v = NULL, labelAll_v = F, labelTop_v = 20, labelDir_v = "both") {
  #' Plot Volcano
  #' @description
  #' Make a volcano plot of DEG results
  #' @param data_dt data.table of gene expression. Can be wrangled for volcano or not. Must have ident1_v and another ident.
  #' @param splitVar_v name to call the variable that the data is split on. If provided. If null, no split. Have only used 'assay' so far.
  #' @param runNames_v values that splitVar_v is split on. If provided.
  #' @param geneCol_v name of column name containing genes. Should be Gene and passed from parent funtion
  #' @param lfc_v log-fold change to include in plot as cut-off
  #' @param pval_v adjusted p-value to include in plot as cut-off
  #' @param colorCol_v column to use for colors. Should be either diffExp for standard, or paste0(splitVar_v, "DE") for comparison ones.
  #' @param title_v optional plot title.
  #' @param verbose_v logical indicating to print messages.
  #' @param labelGenes_v passed to volcanoWrangleMarkers() optional set of genes to label (regardless of significance)
  #' @param labelAll_v passed to volcanoWrangleMarkers() (logical indicating to plot all sig genes.)
  #' @param labelTop_v passed to volcanoWrangleMarkers() (if labelAll_v == F, this will come into play and set number of sig genes to label.)
  #' @param labelDir_v passed to volcanoWrangleMarkers() (can either be 'both', 'up', or 'down', indicating with DEG direction(s) to take top genes from)
  #' @details
  #' Make a volcano plot comparing the DEGs of two different groups. Can be one findmarker result, or can plot two sets of results using splitVar_v
  #' @return volcano plot
  #' @export
  
  ### Prefix/Suffixes to search for in column names
  ixes_v <- paste(unlist(sapply(runNames_v, function(y) paste(c(paste0(c("\\.", "_"), y), paste0(y, c("\\.", "_"))), collapse = "|"), simplify = F)), collapse = "|")
  
  ### Determine if data is wrangled and wide or long
  if (!("diffExp" %in% colnames(data_dt))) {
    
    if (verbose_v) cat(sprintf("diffExp not found in data_dt, indicating it hasn't been wrangled for volcano.\nWrangling now.\n"))
    doubleCols_v <- as.data.table(table(gsub(ixes_v, "", colnames(data_dt))))[N > 1,.N]
    
    if (doubleCols_v > 0) {
      
      if (verbose_v) cat(sprintf("Non-unique columns when splitVar prefix/suffix removed, indicating wide format.\nConverting to long.\n"))
      data_lsdt <- sapply(runNames_v, function(x) data_dt[,mget(c(geneCol_v, grep(x, colnames(data_dt), value = T)))], simplify = F)
      data_lsdt <- sapply(seq_along(data_lsdt), function(i) {
        dat <- data_lsdt[[i]]
        dat <- dat[,mget(grep("pct\\.", colnames(dat), invert = T, value = T))]
        if (!is.null(splitVar_v)) dat[[splitVar_v]] <- runNames_v[i]
        colnames(dat) <- gsub(ixes_v, "", colnames(dat))
        return(dat)}, simplify = F, USE.NAMES = T)
      data_dt <- do.call(rbind, data_lsdt)
      
    } # fi doubleCols_v
    
    data_dt <- wrh.scRNA::volcanoWrangleMarkers(data_dt = data_dt, lfc_v = lfc_v, pval_v = pval_v,
                                                labelGenes_v = labelGenes_v, labelAll_v = labelAll_v, labelTop_v = labelTop_v, labelDir_v = labelDir_v)
    
    if (!is.null(splitVar_v)) {
      data_dt[[paste0(splitVar_v, "DE")]] <- paste(data_dt[[splitVar_v]], data_dt$diffExp)
    }
    
  } # fi !diffExp
  
  ### Get colors
  if (colorCol_v == "diffExp") {
    colors_v <- wrh.scRNA::oneGroupVolcanoColors_v
    labels_v <- names(colors_v)
  } else if (is.null(splitVar_v)) {
    stop("Colors column is not diffExp, but splitVar_v is NULL.\n")
  } else if (colorCol_v == paste0(splitVar_v, "DE")) {
    colors_v <- wrh.scRNA::twoGroupVolcanoColors_v
    labels_v <- gsub("_2", runNames_v[2], gsub("_1", runNames_v[1], names(colors_v)))
  } else {
    stop(sprintf("Colors column must be diffExp or paste0(splitVar_v, 'DE'). %s provided.\n", colorCol_v))
  }
  
  ### Make plot
  plot_gg <- ggplot2::ggplot(data = data_dt, aes(x = avg_log2FC, y = -log10(p_val_adj), col = !!sym(colorCol_v), label = DElabel)) +
    geom_point() + my_theme() +
    theme(legend.position = "right") +
    geom_label_repel(show.legend = F, max.overlaps = Inf, label.padding = 0.1) +
    geom_vline(xintercept=c(-lfc_v, lfc_v), col="black", linetype = "dashed") +
    geom_hline(yintercept=-log10(pval_v), col="black", linetype = "dashed") +
    scale_color_manual(values=colors_v, labels = c("NO", "DOWN", "UP")) +
    guides(color = guide_legend(title = paste0(ident1_v, " Dir")))
  
  ### Add titile
  if (!is.null(title_v)) plot_gg <- plot_gg + ggtitle(title_v)
  
  ### Return
  return(plot_gg)
  
} # plotVolcano
