displaySummary <- function(obj, subObj = NULL, name_v = NULL, subCol_v = "sPop", batch_v, cPops_v = F, popCol_v = "mPop",
                           groupBy_v = list("Population", "Treatment", "Cluster"), alpha_v = 0.5,
                           summary_lsv = list(treat = "Treatment", pop = "mPop", combo = c("Treatment", "mPop")), 
                           outDir_v = NULL, outName_v = NULL, print_v = T, displayOrder_v = "decreasing", cutOff_v = cutOff_v) {
  #' Display Summary Info
  #' @description
    #' Display summary tables and umaps of provided object.
  #' @param obj seurat object to summarize
  #' @param subObj subset object. need to update getSummaryTables to work on subObj, so i don't need both of these
  #' @param name_v population name. will remove "batch#", if present. See details.
  #' @param subCol_v Column to use in conjunction with name_v to subset data to summarize.
  #' @param batch_v name of batch
  #' @param cPops_v logical indicating if the sub-populations are the "collapsed" pops or not. Used to get appropriate colors.
  #' @param popCol_v column name for population annotations
  #' @param groupBy_v vector containing any or all of 'Population', 'Treatment', and 'Cluster', indicating which UMAPs to make.
  #' @param alpha_v alpha parameter only used for by-treatment plot.
  #' @param summary_lsv list of vectors of treatments to summarize. See details.
  #' @param outDir_v optional output directory to save plots and data
  #' @param outName_v optional name for output file.
  #' @param print_v logical indicating to print tables and plots to console.
  #' @param displayOrder_v either 'decreasing' to indicate that identity with most number of cells goes on the bottom, or 'factor' to go by factor levels.
  #' @param cutOff_v numeric value indicating minimum number of cells needed in each group.
  #' @details
    #' This function calls getSummaryTables(), formats those tables, and makes umaps
    #' of seurat clusters and population identities.
    #' 
    #' summary_lsv determines how the data are summarized. It is a named list where
    #' each element's name is user defined and should describe what the output is summarizing.
    #' Refer to the defaults for help. The values of each element correspond to the column or
    #' columns that will be subset from obj@meta.data in order to generate that summary.
    #' 
    #' The default will summarize Treatment alone, major population alone, and treat x pop
    #' together. Right now can only handle 'treat', 'pop' and 'combo'
  #' @return returns list of summary tables as well as two optional outputs: 
  #' 1. save plots and data to file 
  #' and/or
  #' 2. print plots to console.
  #' @export
  
  ###
  ### Wrangle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Get colors, treatments, cluster name, umap name, and subset object
  
  name_v <- gsub("^.*_", "", name_v)
  subBatch_v <- gsub("atch", "", batch_v)
  if (name_v == 'full') {
    subset_v <- NULL
  } else if (name_v == "lymphoid") {
    subset_v <- "Lymphoid/NK"
  } else {
    subset_v <- simpleCap(name_v)
  }
  
  clustName_v <- paste0(subBatch_v, "_clusterNames_v")
  umapName_v <- paste0(subBatch_v, "_umapNames_v")
  if (name_v == 'full') {
    colorName_v <- "mPopColors_v"
  } else {
    if (cPops_v) {
      colorName_v <- paste0(subBatch_v, "_", name_v, "2Colors_v")
    } else {
      colorName_v <- paste0(subBatch_v, "_", name_v, "Colors_v")
    } # fi
  } # fi
  
  clust_v <- eval(as.name(clustName_v)); clust_v <- clust_v[which(names(clust_v) == name_v)]
  umap_v <- eval(as.name(umapName_v)); umap_v <- umap_v[which(names(umap_v) == name_v)]
  
  if (!exists(colorName_v)) {
    colors_v <- NULL
  } else {
    colors_v <- eval(as.name(colorName_v))
  }
  
  treats_v <- intersect(eval(as.name(paste0(subBatch_v, "_treats_v"))), unique(obj$Treatment))
  if (class(obj@meta.data$Treatment) != "factor") obj$Treatment <- factor(obj$Treatment, levels = treats_v)
  
  if (displayOrder_v == "decreasing") {
    counts_dt <- as.data.table(table(obj$Treatment)); setorder(counts_dt, N)
    order_v <- counts_dt$V1
    # obj$displayTreat <- as.character(obj$Treatment)
    # obj$displayTreat <- factor(obj$displayTreat, levels = order_v)
  } else if (displayOrder_v == "factor") {
    # obj$displayTreat <- obj$Treatment
    order_v <- treats_v
  } else {
    stop(sprintf("Either 'factor' or 'decreasing' must be supplied for displayTreat_v. Provided: %s\n", displayOrder_v))
  } # fi displayOrder_v
  
  if (is.null(subObj)) subObj <- obj
  
  if (length(setdiff(names(colors_v), unique(subObj@meta.data[[popCol_v]]))) > 0) {
    names(colors_v) <- gsub("\\/", "-", gsub(" ", ".", names(colors_v)))
  }
  
  ###
  ### Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Call summary fxn
  summary_lsdt <- getSummaryTables(obj = obj, subset_v = subset_v, subCol_v = subCol_v,
                                   summary_lsv = summary_lsv)
  
  ### Make plots
  plots_lsgg <- list()
  
  if ("Population" %in% groupBy_v) {
    plots_lsgg[["Population"]] <- DimPlot(subObj, reduction = umap_v, group.by = popCol_v, pt.size = 0.2, label = F, cols = colors_v) + coord_equal() +
      ggtitle(paste0(batch_v, " ", name_v, " Populations")) + 
      umapFigureTheme()
      
  } # fi pop
  
  if ("Treatment" %in% groupBy_v) {
    treats_gg <- DimPlot(subObj, reduction = umap_v, group.by = "Treatment", pt.size = 0.2, order = order_v, label = F, cols = wrh.scRNA::treatColors_v, alpha = alpha_v) + coord_equal() +
      ggtitle(paste0(batch_v, " ", name_v, " Populations")) +
      umapFigureTheme() + theme(legend.position = "none")
    
    newLegend_gg <- g_legend(DimPlot(subObj, reduction = umap_v, group.by = "Treatment", pt.size = 0.2, cols = wrh.scRNA::treatColors_v) + umapFigureTheme())
    
    outTreats_gg <- ggpubr::ggarrange(plotlist = list(treats_gg, newLegend_gg), ncol = 2, nrow = 1, widths = c(4,1))
    
    plots_lsgg[["Treatment"]] <- outTreats_gg
    
  } # fi treat
  
  if ("Cluster" %in% groupBy_v) {
    plots_lsgg[["Cluster"]] <- DimPlot(subObj, reduction = umap_v, group.by = clust_v, pt.size = 0.2, label = T) + coord_equal() +
      ggtitle(paste0(batch_v, " ", name_v, " Clusters\n", gsub("seurat_clusters_", "", clustName_v))) +
      umapFigureTheme()
  } # fi clust
  
  if (!is.null(outDir_v)) {
    if (is.null(outName_v)) { outName_v <- paste0(batch_v, "_", name_v, "_umap.pdf")}
    pdf(file = file.path(outDir_v, outName_v), onefile = T, width = 10, height = 10)
    invisible(sapply(plots_lsgg, print))
    dev.off()
  } # fi outDir
  
  if (print_v) {
    invisible(sapply(plots_lsgg, print))
  } # fi print
  
  ### Make Tables
  if ('treat' %in% names(summary_lsdt)) {
    print(wrh.rUtils::myKable(data = summary_lsdt$treat, width_v = 20,
                              caption = paste0(name_v, " Cells by Treat")) %>%
            kable_styling(position = "float_left"))
  } # fi treat
  
  if ('pop' %in% names(summary_lsdt)) {
    print(wrh.rUtils::myKable(data = summary_lsdt$pop, width_v = 50,
                              caption = paste0(name_v, " Cells by Pop")) %>%
            kable_styling(position = "left"))
  } # fi treat
  
  if ('combo' %in% names(summary_lsdt)) {
    
    if (is.logical(all.equal(class(summary_lsdt$combo), "list"))) {
      colnames_v <- unique(unlist(sapply(summary_lsdt$combo, colnames)))
      summary_lsdt$combo <- do.call(rbind, sapply(names(summary_lsdt$combo), function(x) {
        y <- summary_lsdt$combo[[x]]
        cols_v <- setdiff(colnames_v, colnames(y))
        for (c_v in cols_v) y[,(c_v) := NA]
        y <- y[which(rowSums(y[,2:ncol(y)], na.rm = T) != 0),]
        y$batch <- x
        return(y)}, simplify = F))
      summary_lsdt$combo <- summary_lsdt$combo[,mget(c("batch", "Pop", setdiff(colnames(summary_lsdt$combo), c("batch", "Pop"))))]
      n_v <- 3
      cols_v <- c("batch", "Pop")
    } else {
      n_v <- 2
      cols_v <- "Pop"
    } # fi
    
    print(wrh.rUtils::kableHighlightColumns(data_dt = summary_lsdt$combo, 
                                cols_v = colnames(summary_lsdt$combo)[n_v:ncol(summary_lsdt$combo)], 
                                condition_v = paste0(">", cutOff_v),
                                caption = paste0(name_v, " Cells by Treat/Pop")))
    
    lowCounts_lsv <- apply(wrh.rUtils::convertDFT(summary_lsdt$combo, col_v = cols_v), 1, function(x) names(which(x <= cutOff_v)))
    lowCounts_lsv <- lowCounts_lsv[which(sapply(lowCounts_lsv, length) > 0)]
    invisible(sapply(names(lowCounts_lsv), function(x) {
      cat(sprintf("Fewer than %s cells for %s\n\t\t%s\n\n", cutOff_v, x, paste(lowCounts_lsv[[x]], collapse = "; ")))
    }))

  } # fi treat
  
  ### Return table
  return(summary_lsdt)
  
} # displaySummary








