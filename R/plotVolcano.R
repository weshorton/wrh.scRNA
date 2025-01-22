plotVolcano <- function(data_dt, splitVar_v = NULL, runNames_v = '', geneCol_v = "Gene", theme_v = "my", mult_v = 1,
                        lfc_v = 0.5, lfcCol_v = "avg_log2FC", pval_v = 0.05, pvalCol_v = "p_val_adj", force_v = 1, thinAxis_v = T,
                        ident1_v, colorCol_v = "diffExp", title_v = NULL, verbose_v = F, labelGenes_v = NULL, hideNonSigLabels_v = F,
                        labelAll_v = F, labelTop_v = NULL, labelDir_v = "both", labelSize_v = 1, labelBox_v = F) {
  #' Plot Volcano
  #' @description
  #' Make a volcano plot of DEG results
  #' @param data_dt data.table of gene expression. Can be wrangled for volcano or not. Must have ident1_v and another ident.
  #' @param splitVar_v name to call the variable that the data is split on. If provided. If null, no split. Have only used 'assay' so far.
  #' @param runNames_v values that splitVar_v is split on. If provided.
  #' @param geneCol_v name of column name containing genes. Should be Gene and passed from parent funtion
  #' @param lfc_v log-fold change to include in plot as cut-off
  #' @param lfcCol_v name of column containing log fold change values. Default is avg_log2FC
  #' @param pval_v adjusted p-value to include in plot as cut-off
  #' @param pvalCol_v name of column contianing adjusted pvalues. Default is p_val_adj
  #' @param force_v fed to force argument in geom_label_repel to push labels apart
  #' @param thinAxis_v logical indicating to override theme and make axes thin
  #' @param colorCol_v column to use for colors. Should be either diffExp for standard, or paste0(splitVar_v, "DE") for comparison ones.
  #' @param title_v optional plot title.
  #' @param verbose_v logical indicating to print messages.
  #' @param labelGenes_v passed to volcanoWrangleMarkers() optional set of genes to label (optionally, regardless of significance)
  #' @param hideNonSigLabels_v option to hide any genes in labelGenes_v that aren't significant. Default is to show all.
  #' @param labelAll_v passed to volcanoWrangleMarkers() (logical indicating to plot all sig genes.)
  #' @param labelTop_v passed to volcanoWrangleMarkers() (if labelAll_v == F, this will come into play and set number of sig genes to label.)
  #' @param labelDir_v passed to volcanoWrangleMarkers() (can either be 'both', 'up', or 'down', indicating with DEG direction(s) to take top genes from)
  #' @param labelSize_v value to increase text size of labels on plot
  #' @param labelBox_v logical indicating to draw the boxes around the label or not
  #' @details
  #' Make a volcano plot comparing the DEGs of two different groups. Can be one findmarker result, or can plot two sets of results using splitVar_v.
  #' The different labeling options are a little complex:
  #' 1. labelAll_v takes precendence over everything. If this is set to T, then all sig genes will be labeled.
  #' 1. labelTop_v is used if labelAll_v == F. This is a numeric value setting the top N significant genes to label.
  #' 1. labelDir_v is used in conjunction with labelTop_v. This determines if the top N up-regulated, down-regulated, or both genes should be labeled.
  #' 1. labelGenes_v can be used independently or in conjunction with labelTop_v. If both are set, both will be labeled. If you only want labelGenes_v to be labeled, then labelTop_v must be NULL
  #' @return volcano plot
  #' @export
  
  ### Volcano colors
  oneGroupVolcanoColors_v <- c("NO" = "grey", "DOWN" = "blue", "UP" = "red")
  twoGroupVolcanoColors_v <- c("NO_1" = "grey", "DOWN_1" = "blue", "UP_1" = "red", "NO_2" = "darkgrey", "DOWN_2" = "darkblue", "UP_2" = "darkred")
  darkVolcanoColors_v <- c("NO" = "darkgrey", "DOWN" = "darkblue", "UP" = "darkred")
  
  ### Prefix/Suffixes to search for in column names
  ixes_v <- paste(unlist(sapply(runNames_v, function(y) paste(c(paste0(c("\\.", "_"), y), paste0(y, c("\\.", "_"))), collapse = "|"), simplify = F)), collapse = "|")
  
  ### Handle theme
  if (theme_v == "my") {
    theme_v <- my_theme()
  } else if (theme_v == "big") {
    theme_v <- big_label()
  } else if (theme_v == "umap") {
    theme_v <- umapFigureTheme()
  } else if (theme_v == "massive") {
    theme_v <- massive_label(multiplier_v = mult_v)
  }
  
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
                                                geneCol_v = geneCol_v, lfcCol_v = lfcCol_v, pvalCol_v = pvalCol_v,
                                                labelGenes_v = labelGenes_v, labelAll_v = labelAll_v, labelTop_v = labelTop_v, labelDir_v = labelDir_v)
    
    if (!is.null(splitVar_v)) {
      data_dt[[paste0(splitVar_v, "DE")]] <- paste(data_dt[[splitVar_v]], data_dt$diffExp)
    }
    
  } # fi !diffExp
  
  ### Get colors
  if (colorCol_v == "diffExp") {
    colors_v <- rev(wrh.scRNA::oneGroupVolcanoColors_v)
    labels_v <- names(colors_v)
  } else if (is.null(splitVar_v)) {
    stop("Colors column is not diffExp, but splitVar_v is NULL.\n")
  } else if (colorCol_v == paste0(splitVar_v, "DE")) {
    colors_v <- wrh.scRNA::twoGroupVolcanoColors_v
    labels_v <- gsub("_2", runNames_v[2], gsub("_1", runNames_v[1], names(colors_v)))
  } else {
    stop(sprintf("Colors column must be diffExp or paste0(splitVar_v, 'DE'). %s provided.\n", colorCol_v))
  }
  
  ## Get label box
  if (labelBox_v) {
    boxLabelSize_v <- 0.25
  } else {
    boxLabelSize_v <- NA
  } # fi
  
  ### Make plot
  plot_gg <- ggplot2::ggplot(data = data_dt, aes(x = avg_log2FC, y = -log10(p_val_adj), col = !!sym(colorCol_v), label = DElabel)) +
    geom_point() + theme_v +
    theme(legend.position = "right") +
    ggrepel::geom_label_repel(size = labelSize_v, show.legend = F, max.overlaps = Inf, label.padding = 0.1, 
                              force = force_v, label.size = boxLabelSize_v) +
    geom_vline(xintercept=c(-lfc_v, lfc_v), col="black", linetype = "dashed") +
    geom_hline(yintercept=-log10(pval_v), col="black", linetype = "dashed") +
    scale_color_manual(values=colors_v, labels = c("NO", "DOWN", "UP")) +
    guides(color = guide_legend(title = paste0(ident1_v, " Dir")))
  
  ### Add titile
  if (!is.null(title_v)) plot_gg <- plot_gg + ggtitle(title_v)
  
  ### Add gene set labels, if required
  if (!is.null(labelGenes_v)) {
    setLabels_v <- setdiff(unique(data_dt$setLabel), "")
    nLabelsFound_v <- length(setLabels_v)
    sigSetLabels_v <- setLabels_v[which(as.character(data_dt[setLabel %in% setLabels_v, setColor]) != "NO")]
    notSigSetLabels_v <- setLabels_v[which(as.character(data_dt[setLabel %in% setLabels_v, setColor]) == "NO")]
    if (hideNonSigLabels_v) {
      data_dt[Gene %in% notSigSetLabels_v, setLabel := ""]
    }
    if (length(labelGenes_v) > 0) {
      
      if (length(setLabels_v) > 0) {
        
        ### Update plot
        plot_gg <- plot_gg + 
          ggnewscale::new_scale_colour() +
          ggrepel::geom_label_repel(data = data_dt, inherit.aes = F, aes(x = avg_log2FC, y = -log10(p_val_adj), col = setColor, label = setLabel),
                           max.overlaps = Inf, label.padding = 0.1, size = labelSize_v, label.size = boxLabelSize_v) +
          scale_colour_manual(values = darkVolcanoColors_v, labels = names(darkVolcanoColors_v)) +
          guides(colour = guide_legend(title = "Gene Set"))
        
      } # fi
        
      ### Determine Significance
      #nSigSetLabels_v <- length(which(as.character(data_dt[setLabel %in% setLabels_v, setColor]) != "NO"))
      
      ### Change theme, if required
      if (thinAxis_v) {
        plot_gg <- plot_gg + theme(axis.line = element_line(linewidth = .25), axis.ticks = element_line(linewidth = .25))
      }
      
      ### Add label
      plot_gg <- ggpubr::annotate_figure(p = plot_gg, 
                                         bottom = text_grob(label = paste0("Found ", nLabelsFound_v, " of ", length(labelGenes_v), 
                                                                           " genes from list (", length(sigSetLabels_v), " significant)"), size = 8))
    } # fi
  } # fi
  
  ### Return
  return(plot_gg)
  
} # plotVolcano
