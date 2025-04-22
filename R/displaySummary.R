displaySummary <- function(obj, 
                           subObj = NULL, 
                           name_v = NULL, 
                           batch_v, 
                           outDir_v = NULL,
                           popCol_v = "mPop", 
                           treatCol_v = "Treatment", 
                           subCol_v = "sPop", 
                           cPops_v = F,
                           summary_lsv = list(treat = "Treatment", pop = "mPop", combo = c("Treatment", "mPop")), 
                           plots_lsv = list("umap" = c("Population", "Treatment", "Cluster", "byTreat"),
                                            "bar" = c("Population", "Treatment", "TreatPop", "PopTreat", "byTreat"),
                                            "stackBar" = c("TreatPop", "PopTreat"),
                                            "facetUMAP" = c("TreatPop")),
                           displayOrder_v = "decreasing", 
                           cutOff_v = 5, 
                           alpha_v = 0.5, 
                           print_v = T, 
                           save_v = T) {
  #' Display Summary Info
  #' @description
    #' Display summary tables and umaps of provided object.
  #' @param obj seurat object to summarize
  #' @param subObj subset object. need to update getSummaryTables to work on subObj, so i don't need both of these
  #' @param name_v population name. will remove "batch#", if present. See details.
  #' @param batch_v name of batch
  #' @param outDir_v optional output directory to save plots and data
  #' @param popCol_v column name for population annotations
  #' @param treatCol_v column name for treatments
  #' @param subCol_v Column to use in conjunction with name_v to subset data to summarize.
  #' @param cPops_v logical indicating if the sub-populations are the "collapsed" pops or not. Used to get appropriate colors.
  #' @param summary_lsv list of vectors of treatments to summarize. See details.
  #' @param plots_lsv list of vectors of plots to make. See details.
  #' @param alpha_v alpha parameter only used for by-treatment plot.
  #' @param displayOrder_v either 'decreasing' to indicate that identity with most number of cells goes on the bottom, or 'factor' to go by factor levels.
  #' @param cutOff_v numeric value indicating minimum number of cells needed in each group.
  #' @param alpha_v alpha parameter only used for by-treatment plot.
  #' @param print_v logical indicating to print tables and plots to console.
  #' @param save_v logical indicating to save plots and data to outDir_v
  #' @details
    #' This function creates a general summmary report for major/minor populations:
      #' 1. Formatted cell count tables from getSummaryTables()
      #' 1. BarCharts of cell counts
      #' 1. UMAP plots colored by specified metadata
    #'
    #' Tables:  
    #' summary_lsv determines how the data are summarized. It is a named list where
    #' each element's name is user-defined (eventually) and should describe what the output is summarizing.
    #' Refer to the defaults for help. The values of each element correspond to the column or
    #' columns that will be subset from obj@meta.data in order to generate that summary.
    #' 
    #' The default will summarize Treatment alone, major population alone, and treat x pop
    #' together. Right now can only handle 'treat', 'pop' and 'combo'
    #' 
    #' Plots:
    #' plots_lsv determines which classes of plots, and which versions of those classes to make. 
    #' So far, the following are supported (shown in format "listElementName" = c("possible", "values"))
      #' 1. 'umap' = c("Population", "Treatment", "Cluster") 
        #' - 'Population' plot is colored using provided popCol_v
        #' - 'Treatment' plot is colored using provided treatCol_v
        #' - 'Cluster' plot is colored using the seurat_clusters value stored in corresponding wrh.scRNA data object.
      #' 2. 'bar' = c("Population", "Treatment", "TreatPop", "PopTreat")
        #' - Population, Treatment, and Cluster all treated the same
        #' - TreatPop: facet by treatCol_v, color by popCol_v
        #' - PopTreat: facet by popCol_v, color by TreatCol_v
        #' - TreatPop and PopTreat both require 'combo' to be selected for summary_lsv.
        #' - byTreat requires combo. makes one plot per treatment (same results as the treatpop facet, but not in facet form)
      #' 3. 'facetUMAP' = c("TreatPop") 
        #' - Treatment_Population: facet by treatCol_v, color by popCol_v
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
  clustName_v <- paste0(subBatch_v, "_clusterNames_v")
  umapName_v <- paste0(subBatch_v, "_umapNames_v")
  resName_v <- paste0(subBatch_v, "_res_v")
  clust_v <- eval(as.name(clustName_v)); clust_v <- clust_v[which(names(clust_v) == name_v)]
  umap_v <- eval(as.name(umapName_v)); umap_v <- umap_v[which(names(umap_v) == name_v)]
  res_v <- eval(as.name(resName_v)); res_v <- res_v[which(names(res_v) == name_v)]
  
  if (name_v %in% c('batch12', 'batch3')) {
    colorName_v <- "mPopColors_v"
  } else if (name_v == "neoplastic") {
    colorName_v <- paste0(subBatch_v, "_", name_v, "Colors_v")
  } else {
    if (cPops_v) {
      colorName_v <- paste0(subBatch_v, "_", name_v, "Colors_v")
    } else {
      colorName_v <- paste0(subBatch_v, "_", name_v, "BadColors_v")
    } # fi
  } # fi
  
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
  } else if (displayOrder_v == "factor") {
    order_v <- treats_v
  } else {
    stop(sprintf("Either 'factor' or 'decreasing' must be supplied for displayTreat_v. Provided: %s\n", displayOrder_v))
  } # fi displayOrder_v
  
  if (is.null(subObj)) subObj <- obj 
  
  if (length(setdiff(names(colors_v), unique(subObj@meta.data[[popCol_v]]))) > 0) {
    names(colors_v) <- gsub("\\/", "-", gsub(" ", ".", names(colors_v)))
  }
  
  ###
  ### Summary Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Call summary fxn
  summary_lsdt <- getSummaryTables(obj = subObj, subset_v = NULL, subCol_v = subCol_v, popCol_v = popCol_v,
                                   treatCol_v = treatCol_v, summary_lsv = summary_lsv)
  
  if (save_v & !is.null(outDir_v)) {
    tmp_lsdt <- summary_lsdt
    if (subBatch_v == "b12") {
      cols_v <- colnames(tmp_lsdt$combo$Total)
      tmp_lsdt$combo <- sapply(names(tmp_lsdt$combo), function(x) {
        xx <- tmp_lsdt$combo[[x]]
        xx$Batch <- x
        for (c_v in cols_v) { if (!c_v %in% colnames(xx)) xx[[c_v]] <- NA}
        xx <- xx[,mget(c(cols_v, "Batch"))]
        return(xx)
      }, simplify = F, USE.NAMES = T)
      tmp_lsdt$combo <- do.call(rbind, tmp_lsdt$combo)
    } # fi subBatch
    wrh.rUtils::writeCSVorExcel(tmp_lsdt,
                                file_v = file.path(outDir_v, paste0(batch_v, "_", name_v, "_cellCounts.xlsx")))
  } # fi save
  
  if (subBatch_v == "b12") {
    summary_lsdt$combo <- summary_lsdt$combo$Total
  }
  
  ###
  ### Cell Bar Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  plots_lsgg <- list()
  
  if ('bar' %in% names(plots_lsv)) {
    
    barToRun_v <- plots_lsv$bar
    bar_lsgg <- list()
    
    if ('Population' %in% barToRun_v) {
      if (!'pop' %in% names(summary_lsdt)) stop("Missing summary table for population barplot.\n")
      temp <- ggplot(data = summary_lsdt$pop, aes(x = Pop, y = nCells, fill = Pop)) +
        geom_bar(stat = 'identity') + ggtitle("Cells per Population") + umapFigureTheme() + 
        scale_fill_manual(values = colors_v, breaks = names(colors_v))
      
      if (grepl("[Ll]ymphoid|[Mm]yeloid|[Ff]ull", name_v)) {
        temp <- temp + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
          theme(legend.position = "bottom") + 
          guides(fill=guide_legend(ncol=2))
      } # fi grep
      
      if (grepl("[Nn]eoplastic", name_v)) {
        temp <- temp + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
      } # fi grep
      bar_lsgg[['Population']] <- temp
    } # fi population
    
    if ('Treatment' %in% barToRun_v) {
      if (!'treat' %in% names(summary_lsdt)) stop("Missing summary table for treatment barplot.\n")
      bar_lsgg[['Treatment']] <- ggplot(data = summary_lsdt$treat, aes(x = Treat, y = nCells, fill = Treat)) +
        geom_bar(stat = 'identity') + ggtitle("Cells per Treatment") + umapFigureTheme() +
        scale_fill_manual(values = treatColors_v, breaks = names(treatColors_v))
    } # fi treatment
    
    if ('TreatPop' %in% barToRun_v) {
      if (!'combo' %in% names(summary_lsdt)) stop("Missing summary table for TreatPop barplot.\n")
      melt_dt <- melt(summary_lsdt$combo, id.vars = "Pop")
      #melt_dt$variable <- factor(as.character(melt_dt$variable), levels = rev(treats_v)) # want treatment reverse
      temp <- ggplot(data = melt_dt, aes(x = Pop, y = value, fill = Pop)) +
        geom_bar(stat = 'identity') + ggtitle("Cells per Population (by Treatment)") + facet_wrap("~variable", ncol = 2) +
        scale_fill_manual(values = colors_v, breaks = names(colors_v)) + umapFigureTheme() + 
        labs(y = "Cell Count")
      
        temp <- temp + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
          theme(legend.position = "bottom") + 
          guides(fill=guide_legend(ncol=2))
      
      bar_lsgg[['TreatPop']] <- temp
    } # fi treatpop
    
    if ('PopTreat' %in% barToRun_v) {
      if (!'combo' %in% names(summary_lsdt)) stop("Missing summary table for PopTreat barplot.\n")
      melt_dt <- melt(summary_lsdt$combo, id.vars = "Pop")
      bar_lsgg[['PopTreat']] <- ggplot(data = melt_dt, aes(x = variable, y = value, fill = variable)) +
        geom_bar(stat = 'identity') + ggtitle("Cells per Treatment (by Population)") + facet_wrap("~Pop", ncol = 2) +
        scale_fill_manual(values = treatColors_v, breaks = names(treatColors_v)) + umapFigureTheme() +
        labs(x = "Treatment", y = "Cell Count")
    } # fi poptreat
    
    if ('byTreat' %in% barToRun_v) {
      if (!'combo' %in% names(summary_lsdt)) stop("Missing summary table for PopTreat barplot.\n")
      data_dt <- melt(summary_lsdt$combo, id.vars = "Pop")
      byTreatBar_lsgg <- sapply(unique(data_dt$variable), function(x) {
        y <- ggplot(data = data_dt[variable == x,], aes(x = Pop, y = value, fill = Pop)) +
          geom_bar(stat = "identity") + ggtitle(paste0("Cells per Population in ", x)) +
          umapFigureTheme() + scale_fill_manual(values = colors_v, breaks = names(colors_v)) +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        return(y)}, simplify = F, USE.NAMES = T)
    } # fi byTreat
    
    if (print_v) {
      invisible(sapply(bar_lsgg, print))
    } # fi print
    
    if (save_v & !is.null(outDir_v)) {
      pdf(file = file.path(outDir_v, paste0(batch_v, "_", name_v, "_cellCountBar.pdf")), onefile = T, height = 10, width = 10)
      invisible(sapply(bar_lsgg[grep("Population|Treatment", names(bar_lsgg))], print))
      dev.off()
      
      pdf(file = file.path(outDir_v, paste0(batch_v, "_", name_v, "_cellCountFacetBar.pdf")), onefile = T, height = 20, width = 16)
      invisible(sapply(bar_lsgg[grep("TreatPop|PopTreat", names(bar_lsgg))], print))
      dev.off()
      
      if ('byTreat' %in% barToRun_v) {
        pdf(file = file.path(outDir_v, paste0(batch_v, "_", name_v, "_indTreat_cellCountBar.pdf")), onefile = T, height = 10, width = 10)
        invisible(sapply(byTreatBar_lsgg, print))
        dev.off()
      }
    } # fi save
    
    plots_lsgg[['bar']] <- bar_lsgg
    
  } # fi bar
  
  ###
  ### Cell Stacked Bar Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  if ('stackBar' %in% names(plots_lsv)) {
    
    if (!'combo' %in% names(summary_lsv)) stop("Missing summary table for stacked barplots")
    
    stackedToRun_v <- plots_lsv$stackBar
    stackedBar_lsgg <- list()
    
    if ('TreatPop' %in% stackedToRun_v) {
      count_dt <- myT(summary_lsdt$combo, newName_v = "Treatment"); prop_dt <- copy(count_dt)
      prop_dt[, Total := rowSums(.SD), .SDcols = setdiff(colnames(prop_dt), "Treatment")]
      cols_v <- setdiff(colnames(prop_dt), c("Treatment", "Total"))
      prop_dt <- prop_dt[, (cols_v) := lapply(.SD, function(x) x / Total * 100), .SDcols = cols_v]
      prop_dt$Total <- NULL
      
      plotData_lsdt <- list("Counts" = melt(count_dt, id.vars = "Treatment"), "Proportions" = melt(prop_dt, id.vars = "Treatment"))
      
      plotData_lsdt <- sapply(plotData_lsdt, function(x) {
        if (class(x$Treatment) == "character") x$Treatment <- factor(x$Treatment, levels = rev(treats_v))
        return(x)}, simplify = F, USE.NAMES = T)
      
      tempStackedBar_lsgg <- sapply(names(plotData_lsdt), function(x) {
        y <- ggplot(data = plotData_lsdt[[x]], aes(x = Treatment, y = value, fill = variable)) +
          geom_bar(position = "stack", stat = "identity") +
          scale_fill_manual(values = colors_v, breaks = names(colors_v)) +
          ggtitle(paste0("Population ", x, "\nPer Treatment")) + umapFigureTheme() +
          labs(y = gsub("s%", "", x), fill = "Population")
        return(y)}, simplify = F, USE.NAMES = T)
      
      treatPopStacked_gg <- ggpubr::ggarrange(plotlist = tempStackedBar_lsgg, ncol = 2, common.legend = T, legend = "bottom") +
        guides(fill=guide_legend(ncol=2))
      stackedBar_lsgg[['TreatPop']] <- treatPopStacked_gg
  
    } # fi TreatPop
    
    if ('PopTreat' %in% stackedToRun_v) {
      count_dt <- summary_lsdt$combo; prop_dt <- copy(count_dt)
      prop_dt[, Total := rowSums(.SD), .SDcols = setdiff(colnames(prop_dt), "Pop")]
      cols_v <- setdiff(colnames(prop_dt), c("Pop", "Total"))
      prop_dt <- prop_dt[, (cols_v) := lapply(.SD, function(x) x / Total * 100), .SDcols = cols_v]
      prop_dt$Total <- NULL
      
      plotData_lsdt <- list("Counts" = melt(count_dt, id.vars = "Pop"), "Proportions" = melt(prop_dt, id.vars = "Pop"))
      
      plotData_lsdt <- sapply(plotData_lsdt, function(x) {
        if (class(x$Pop) == "character") {
          warning("stackedBar, popTreat, population column is character. haven't tested this yet so check output.")
          x$Pop <- factor(x$Pop, levels = unique(x$Pop)[order(as.numeric(gsub("neo.c", "", x$Pop)))])
        }
        return(x)}, simplify = F, USE.NAMES = T)
      
      tempStackedBar_lsgg <- sapply(names(plotData_lsdt), function(x) {
        y <- ggplot(data = plotData_lsdt[[x]], aes(x = Pop, y = value, fill = variable)) +
          geom_bar(position = "stack", stat = "identity") +
          scale_fill_manual(values = treatColors_v, breaks = names(treatColors_v)) +
          ggtitle(paste0("Treatment ", x, "\nPer Population")) + umapFigureTheme() +
          labs(y = gsub("s%", "", x), fill = "Treatment") + 
          angle_x()
          #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
        return(y)}, simplify = F, USE.NAMES = T)
      
      popTreatStacked_gg <- ggpubr::ggarrange(plotlist = tempStackedBar_lsgg, ncol = 2, common.legend = T, legend = "bottom") +
        guides(fill=guide_legend(ncol=2))
      stackedBar_lsgg[['PopTreat']] <- popTreatStacked_gg
    } # fi TreatPop
    
    if (print_v) {
      invisible(sapply(stackedBar_lsgg, print))
    } # fi print
    
    if (save_v & !is.null(outDir_v)) {
      pdf(file = file.path(outDir_v, paste0(batch_v, "_", name_v, "_cellCountStackedBar.pdf")), onefile = T, width = 16, height = 10)
      invisible(sapply(stackedBar_lsgg, print))
      dev.off()
    } # fi save
    
    plots_lsgg[['stackBar']] <- stackedBar_lsgg
    
  } # fi stackBar
  
  ###
  ### Standard UMAPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###

  ### Maybe try to use standardPlots sometime, but not right now.  
  # if ('umap' | 'facetUMAP' %in% names(plots_lsv)) {
  #   
  #   tmpStdPlots_lsgg <- standardPlots(seurat_obj = subObj, reduction_v = umap_v, clustCol_v = clust_v, res_v = res_v, name_v = name_v,
  #                                     pt.size_v = 0.5, featurePlots_v = T, maxQ_v = "q90")
  #   
  #   umap_lsgg <- list()
  #   umap_lsgg[["Population"]] <- tmpStdPlots_lsgg$clusters
  
  if ('umap' %in% names(plots_lsv)) {
    
    umapToRun_v <- plots_lsv$umap
    umap_lsgg <- list()
    
    if ("Population" %in% umapToRun_v) {
      umap_lsgg[["Population"]] <- DimPlot(subObj, reduction = umap_v, group.by = popCol_v, order = rev(names(colors_v)), pt.size = 0.5, label = F, cols = colors_v) + coord_equal() +
        ggtitle(paste0(batch_v, " ", name_v, " Populations")) + 
        umapFigureTheme() + theme(legend.position = "bottom") + 
        guides(color=guide_legend(override.aes = list(size = 6), ncol=2))
    } # fi pop
    
    if ("Treatment" %in% umapToRun_v) {
      treats_gg <- DimPlot(subObj, reduction = umap_v, group.by = treatCol_v, pt.size = 0.5, order = order_v, label = F, cols = wrh.scRNA::treatColors_v, alpha = alpha_v) + coord_equal() +
        ggtitle(paste0(batch_v, " ", name_v, " Populations")) +
        umapFigureTheme() + theme(legend.position = "none")
      
      legendPlot_gg <- DimPlot(subObj, reduction = umap_v, group.by = treatCol_v, pt.size = 0.5, cols = wrh.scRNA::treatColors_v, order = treats_v) + umapFigureTheme()
      newLegend_gg <- g_legend(legendPlot_gg)
      
      outTreats_gg <- ggpubr::ggarrange(plotlist = list(treats_gg, newLegend_gg), ncol = 2, nrow = 1, widths = c(4,1))
      
      umap_lsgg[["Treatment"]] <- outTreats_gg
      
    } # fi treat
    
    if ("Cluster" %in% umapToRun_v) {
      umap_lsgg[["Cluster"]] <- DimPlot(subObj, reduction = umap_v, group.by = clust_v, pt.size = 0.5, label = T) + coord_equal() +
        ggtitle(paste0(batch_v, " ", name_v, " Clusters\n", gsub("seurat_clusters_", "", clustName_v))) +
        umapFigureTheme()
    } # fi clust
    
    if ('byTreat' %in% umapToRun_v) {
      
      byTreatDim_lsgg <- sapply(treats_v, function(x) {
        cells_v <- rownames(subObj@meta.data[subObj@meta.data[[treatCol_v]] == x,])
        dim_gg <- DimPlot(object = subObj,
                          cells = cells_v,
                          reduction = umap_v,
                          group.by = popCol_v,
                          label = F,
                          cols = colors_v,
                          order = rev(names(colors_v)),
                          pt.size = 0.5) + coord_equal() +
          ggtitle(paste0("Cells per Population in ", x)) + 
          umapFigureTheme() + theme(legend.position = "bottom")
        return(dim_gg)}, simplify = F, USE.NAMES = T)
    } # fi byTreat
    
    if (print_v) {
      invisible(sapply(umap_lsgg, print))
    } # fi print
    
    if (save_v & !is.null(outDir_v)) {
      pdf(file = file.path(outDir_v, paste0(batch_v, "_", name_v, "_umap.pdf")), onefile = T, width = 10, height = 10)
      invisible(sapply(umap_lsgg, print))
      dev.off()
      
      if ('byTreat' %in% umapToRun_v) {
        pdf(file = file.path(outDir_v, paste0(batch_v, "_", name_v, "_indTreat_umapColorPop.pdf")), onefile = T, width = 10, height = 10)
        invisible(sapply(byTreatDim_lsgg, print))
        dev.off()
      }
    } # fi save
    
    plots_lsgg[['umap']] <- umap_lsgg
    
  } # fi umap
  
  ###
  ### Facet UMAP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  if ('facetUMAP' %in% names(plots_lsv)) {
    
    facetToRun_v <- plots_lsv$facetUMAP
    facetUMAP_lsgg <- list()
    
    if ("TreatPop" %in% facetToRun_v) {
      
      subObj@meta.data[[treatCol_v]] <- factor(as.character(subObj@meta.data[[treatCol_v]]), levels = rev(treats_v))
      
      facetUMAP_lsgg[["TreatPop"]] <- DimPlot(subObj, reduction = umap_v, group.by = popCol_v, split.by = treatCol_v, pt.size = 0.4, label = F, cols = colors_v, ncol = 2) +
        ggtitle(paste0(batch_v, " ", name_v, " Populations")) + #umapFigureTheme() +
        theme(legend.position = "right") + 
        guides(color=guide_legend(override.aes = list(size = 6), ncol=1))
      
    } # fi treatPop
    
    if (print_v) {
      invisible(sapply(facetUMAP_lsgg, print))
    } # fi print
    
    if (save_v & !is.null(outDir_v)) {
      pdf(file = file.path(outDir_v, paste0(batch_v, "_", name_v, "_facetUMAP.pdf")), onefile = T, height = 12, width = 8)
      invisible(sapply(facetUMAP_lsgg, print))
      dev.off()
    } # fi save
    
    plots_lsgg[['facetUMAP']] <- facetUMAP_lsgg
    
  } # fi facetUMAP

  ###
  ### Format Summary Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Make Tables
  if ('treat' %in% names(summary_lsdt) & print_v) {
    print(wrh.rUtils::myKable(data = summary_lsdt$treat, width_v = 20,
                              caption = paste0(name_v, " Cells by Treat")) %>%
            kable_styling(position = "float_left"))
  } # fi treat
  
  if ('pop' %in% names(summary_lsdt) & print_v) {
    print(wrh.rUtils::myKable(data = summary_lsdt$pop, width_v = 50,
                              caption = paste0(name_v, " Cells by Pop")) %>%
            kable_styling(position = "left"))
  } # fi treat
  
  if ('combo' %in% names(summary_lsdt) & print_v) {
    
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
  return(list("data" = summary_lsdt, "plots" = plots_lsgg))
  
} # displaySummary








