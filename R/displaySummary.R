displaySummary <- function(obj, subObj = NULL, name_v = NULL, subCol_v = "mPop", batch_v,
                           summary_lsv = list(treat = "Treatment", pop = "mPop", combo = c("Treatment", "mPop")), 
                           outDir_v = NULL, print_v = T, cutOff_v = cutOff_v) {
  #' Display Summary Info
  #' @description
    #' Display summary tables and umaps of provided object.
  #' @param obj seurat object to summarize
  #' @param subObj subset object. need to update getSummaryTables to work on subObj, so i don't need both of these
  #' @param name_v population name. will remove "batch#", if present. See details.
  #' @param subCol_v Column to use in conjunction with name_v to subset data to summarize.
  #' @param batch_v name of batch
  #' @param summary_lsv list of vectors of treatments to summarize. See details.
  #' @param outDir_v optional output directory to save plots and data
  #' @param print_v logical indicating to print tables and plots to console.
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
  #' @return no object returned. Two optional outputs: 
  #' 1. save plots and data to file 
  #' and/or
  #' 2. print plots to console.
  #' @export
  
  ### Handle name stuff
  name_v <- gsub("^.*_", "", name_v)
  subBatch_v <- gsub("atch", "", batch_v)
  if (name_v == 'full') {
    subset_v <- NULL
  } else if (name_v == "lymphoid") {
    subset_v <- "Lymphoid/NK"
  } else {
    subset_v <- simpleCap(name_v)
  }
  
  ### Get package data
  clustName_v <- paste0(subBatch_v, "_clusterNames_v")
  umapName_v <- paste0(subBatch_v, "_umapNames_v")
  if (name_v == 'full') {
    colorName_v <- "mPopColors_v"
  } else {
    colorName_v <- paste0(subBatch_v, "_", name_v, "Colors_v")
  }
  
  clust_v <- eval(as.name(clustName_v)); clust_v <- clust_v[which(names(clust_v) == name_v)]
  umap_v <- eval(as.name(umapName_v)); umap_v <- umap_v[which(names(umap_v) == name_v)]
  
  if (!exists(colorName_v)) {
    colors_v <- NULL
  } else {
    colors_v <- eval(as.name(colorName_v))
  }
  
  ### Call summary fxn
  summary_lsdt <- getSummaryTables(seurat_obj = obj, subset_v = subset_v, subCol_v = subCol_v,
                                   summary_lsv = summary_lsv)
  
  ### Make plots
  if (is.null(subObj)) subObj <- obj
  
  cluster_gg <- DimPlot(subObj, reduction = umap_v, group.by = clust_v, pt.size = 0.1, label = T) + coord_equal() +
    ggtitle(paste0(batch_v, " ", name_v, " Clusters\n", gsub("seurat_clusters_", "", clustName_v)))
  
  pop_gg <- DimPlot(subObj, reduction = umap_v, group.by = subCol_v, pt.size = 0.1, label = T) + coord_equal() +
    ggtitle(paste0(batch_v, " ", name_v, " Populations"))
  
  if (!is.null(outDir_v)) {
    pdf(file = file.path(outDir_v, paste0(batch_v, "_", name_v, "_", subCol_v, "_umap.pdf")), onefile = T)
    print(cluster_gg)
    print(pop_gg)
    dev.off()
  } # fi outDir
  
  if (print_v) {
    print(cluster_gg)
    print(pop_gg)
  } # fi print
  
  ### Make Tables
  if ('treat' %in% names(summary_lsv)) {
    print(wrh.rUtils::myKable(data = summary_lsdt$treat, width_v = 20,
                              caption = paste0(name_v, " Cells by Treat")) %>%
            kable_styling(position = "float_left"))
  } # fi treat
  
  if ('pop' %in% names(summary_lsv)) {
    print(wrh.rUtils::myKable(data = summary_lsdt$pop, width_v = 50,
                              caption = paste0(name_v, " Cells by Pop")) %>%
            kable_styling(position = "left"))
  } # fi treat
  
  if ('combo' %in% names(summary_lsv)) {
    print(wrh.rUtils::kableHighlightColumns(data_dt = summary_lsdt$combo, 
                                cols_v = colnames(summary_lsdt$combo)[2:ncol(summary_lsdt$combo)], 
                                condition_v = paste0(">", cutOff_v),
                                caption = paste0(name_v, " Cells by Treat/Pop")))
    
    lowCounts_lsv <- apply(wrh.rUtils::convertDFT(summary_lsdt$combo, col_v = "Pop"), 1, function(x) names(which(x <= cutOff_v)))
    lowCounts_lsv <- lowCounts_lsv[which(sapply(lowCounts_lsv, length) > 0)]
    invisible(sapply(names(lowCounts_lsv), function(x) {
      cat(sprintf("Fewer than %s cells for %s\n\t\t%s\n\n", cutOff_v, x, paste(lowCounts_lsv[[x]], collapse = "; ")))
    }))

  } # fi treat
  
} # displaySummary








