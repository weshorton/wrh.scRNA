myClusterSweep <- function(seurat_obj, 
                           embedding_v = seurat_obj@reductions$pca@cell.embeddings,
                           reduction_v = "pca",
                           reductionName_v = "umap",
                           ndims_v = 10,
                           res_v = seq(from=0.1, to=1, by=0.1),
                           verbose_v = F,
                           indPlots_v = F,
                           plotDir_v = NULL,
                           name_v = "",
                           print_v = T) {
  #' Cluster Parameter Sweep
  #' @description 
  #' Run a sweep on different resolutions for clustering. Output plots to assess quality of clustering.
  #' Adapted from code provided by Nick Calistri.
  #' @param seruat_obj A seurat object with dimensional reduction embeddings calculated (PCA)
  #' @param embedding_v embeddings for chosen dimensional reduction. Usually: seruat_obj@reductions$pca@cell.embeddings
  #' @param reduction_v which reduction to use. Default is pca (which goes with default embedding). Be sure to change both!
  #' @param reductionName_v name to call reduction on seurat object. Default is "umap", other values likely "umapR" or "umapH" for rliger/harmony
  #' @param ndims_v number of dimensions to use. Default is 10.
  #' @param res_v vector of resolution values to pass to FindClusters()
  #' @param verbose_v logical indicating whether to print seurat function messages or not.
  #' @param indPlots_v logical indicating whether to print and/or save the individual QC plots
  #' @param plotDir_v path to output directory for plots. If NULL, will print to console.
  #' @param name_v Use as a prefix for output plots. Sample1 for example.
  #' @param print_v logical. Should main plots be printed, even if they're also saved? Does not change individual plots
  #' @return Either save plots to file or print to console.
  #' @export
  
  ###
  ### PREP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Remove any clustering columns that are already in the object, to make sure that they don't interfere with QC run
  cols_v <- grep("_res", colnames(seurat_obj@meta.data), value = T)
  if (length(cols_v) > 0) {
    for (col_v in cols_v) seurat_obj@meta.data[[col_v]] <- NULL
  }
  
  ###
  ### CALCULATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  origNames_v <- names(seurat_obj)
  origCols_v <- colnames(seurat_obj@meta.data)
  
  ### Find neighbors using specified number of dimensions
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = reduction_v, dims = 1:ndims_v, verbose = verbose_v)
  
  newNames1_v <- names(seurat_obj)
  newCols1_v <- colnames(seurat_obj@meta.data)
  
  ### Find clusters using the specified resolutions
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = res_v, verbose = verbose_v)
  
  newNames2_v <- names(seurat_obj)
  newCols2_v <- colnames(seurat_obj@meta.data)
  
  ### Calculate UMAP
  if (!reductionName_v %in% names(seurat_obj)) {
    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:ndims_v, reduction = reduction_v, reduction.name = reductionName_v, verbose = F)
  }
  
  ### Compute QC metrics
  seuratClusterQC <- clusterQC(seurat_obj = seurat_obj, embedding_v = embedding_v, ndims_v = ndims_v, reductionName_v = reductionName_v, reduction_v = reduction_v)
  
  ### Melt for plotting
  meltClusterQC <- reshape2::melt(seuratClusterQC$QC, id.vars = "res")
  meltClusterQC$variable <- factor(meltClusterQC$variable, levels = rev(levels(meltClusterQC$variable)))
  
  ### Determine x intercepts
  silX_v <- as.numeric(seuratClusterQC$QC[which.max(seuratClusterQC$QC$meanSilWidth),"res"])
  rmsdX_v <- as.numeric(seuratClusterQC$QC[which.min(seuratClusterQC$QC$meanRMSD),"res"])
  
  ###
  ### OUTPUT INDIVIDUAL PLOTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  if (indPlots_v) {
    if (!is.null(plotDir_v)) {
      
      wrh.rUtils::mkdir(baseDir_v = dirname(plotDir_v), newDir_v = basename(plotDir_v))
      
      ### Dim Plots
      pdf(file = file.path(plotDir_v, paste0(name_v, "_dimPlots.pdf")), onefile = T)
      for (i in 1:length(seuratClusterQC$dimPlot)) print(seuratClusterQC$dimPlot[[i]])
      dev.off()
      
      ###  Plots
      pdf(file = file.path(plotDir_v, paste0(name_v, "_silBoxPlots.pdf")), onefile = T)
      for (i in 1:length(seuratClusterQC$silBox)) print(seuratClusterQC$silBox[[i]])
      dev.off()
      
      ### Dim Plots
      pdf(file = file.path(plotDir_v, paste0(name_v, "_barPlots.pdf")), onefile = T)
      for (i in 1:length(seuratClusterQC$RMSDBar)) print(seuratClusterQC$RMSDBar[[i]])
      dev.off()
      
    } else {
      
      print(seuratClusterQC$dimPlot)
      print(seuratClusterQC$silBox)
      print(seuratClusterQC$RMSDBar)
      
    } # fi is.null(plotDir_v)
  } # fi indPlots_v
  
  ###
  ### OUTPUT SUMMARY PLOTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Plot Sil Width, RMSD, and nClust by resolution
  compareVals_gg <- ggplot(meltClusterQC, aes(x = res, y = value, color = variable)) +
    geom_point() + geom_line() + theme_bw() +
    facet_wrap(~variable, ncol = 3, scales = 'free_y') +
    scale_x_continuous(breaks = seq(from = 0, to = 2.5, by = 0.5)) +
    xlab("Louvain resolution") +
    geom_vline(aes(xintercept=silX_v, color = "Max Sil"), linetype = "dashed") +
    geom_vline(aes(xintercept=rmsdX_v, color = "Min RMSD"), linetype = "dashed") +
    geom_vline(xintercept = silX_v, linetype = 'dashed', color = "red") +
    geom_vline(xintercept = rmsdX_v, linetype = 'dashed', color = "blue") +
    ggtitle(paste0("Sweep Results For ", name_v))
    
  ### Scatterplot of width vs RMSD
  valScatter_gg <- ggplot(seuratClusterQC$QC, aes(x = meanSilWidth, y = meanRMSD, label = res, color = nClust)) +
    geom_point() + ggrepel::geom_text_repel() + theme_bw() +
    ggtitle(paste0("Silhouette Width by RMSD - ", name_v))
  
  ### Scatterplot of nClust vs width/rmsd metric
  bestMetric_gg <- ggplot(seuratClusterQC$QC, aes(x = nClust, y = meanSilWidth/meanRMSD, label = res, color = nClust)) +
    geom_point() + ggrepel::geom_text_repel() + theme_bw() +
    ggtitle(paste0("nClust by width/rmsd metric - ", name_v))
  
  ### Output
  if (!is.null(plotDir_v)) {
    
    ggsave(file.path(plotDir_v, paste0(name_v, "_compSilxRMSD.pdf")), compareVals_gg)
    ggsave(file.path(plotDir_v, paste0(name_v, "_silxRMSDScatter.pdf")), valScatter_gg)
    ggsave(file.path(plotDir_v, paste0(name_v, "_metricScatter.pdf")), bestMetric_gg)
    
    if (print_v) {
      print(compareVals_gg)
      print(valScatter_gg)
      print(bestMetric_gg)
    } # fi print_v
    
  } else {
    print(compareVals_gg)
    print(valScatter_gg)
    print(bestMetric_gg)
  } # fi !is.null(plotDir_v)
  
} # clusterSweep

clusterQC <- function(seurat_obj, embedding_v, ndims_v, reductionName_v, reduction_v) {
  #' Cluster QC
  #' @description 
  #' Given a seurat object and embeddings, calculate qc
  #' @param seurat_obj A seurat object with dimensional reduction embeddings calculated (PCA)
  #' @param embedding_v embeddings for chosen dimensional reduction. Usually: seruat_obj@reductions$pca@cell.embeddings
  #' @param ndims_v number of dimensions to use. Default is 10.
  #' @param reductionName_v name of reduction to use.
  #' @param reduction_v reduction to use
  #' @return List containing the following elements: 
  #'        "QC" - tibble of mean silhouette and RMSD values for each resolution
  #'        "dimPlot" - list of DimPlots, one for each resolution
  #'        "silBox" - list of silhouette boxplots, one for each resolution
  #'        "RMSDBar" - list of RMSD barplots, one for each resolution
  #' @export
  
  ### Grab names of each of the clustering runs
  uniqClusterRunNames_v <- grep("*_res", colnames(seurat_obj@meta.data), value = T)
  
  ### Get a data.frame of all the clustering assignments
  clusterCalls_df <- seurat_obj@meta.data[,uniqClusterRunNames_v, drop = F]
  
  ### Subset embedding
  embedding_v <- embedding_v[,1:ndims_v]
  
  ### Empty list to hold plots
  dim_ls <- box_ls <- bar_ls <- list()
  
  ### Iterate over each and make QC plots
  for (i in 1:length(uniqClusterRunNames_v)) {
    
    ### Get cluster name and values
    currClusterName_v <- uniqClusterRunNames_v[i]
    currClusterCalls_v <- clusterCalls_df[,currClusterName_v]
    
    ### Grab resolution from cluster name
    currRes_v <- as.numeric(stringr::str_remove(currClusterName_v,
                                       pattern = paste0(Seurat::DefaultAssay(seurat_obj), '_snn_res.')))
    
    ### Make reduction name
    #currRedName_v <- paste0(reductionName_v, currRes_v)
    
    ### Calculate UMAP - don't need to recalculate...
    # seurat_obj <- Seurat::RunUMAP(seurat_obj, 
    #                               dims = 1:ndims_v, 
    #                               reduction = reduction_v, 
    #                               reduction.name = reductionName_v, 
    #                               verbose = F)
    
    ### Create DimPlot with clusters labelled.
    dim_ls[[i]] <- Seurat::DimPlot(seurat_obj,
                           group.by = currClusterName_v,
                           reduction = reductionName_v,
                           label = T, pt.size = 0.1) +
      coord_equal() +
      ggtitle(currClusterName_v)
    
    ### Calculate silhouettes for cluster
    currSil <- as.data.frame(bluster::approxSilhouette(x = embedding_v,
                                              clusters = currClusterCalls_v))
    
    ### Boxplot of silhouette widths for each cluster
    box_ls[[i]] <- ggplot(data = currSil, aes(x = cluster, y = width)) +
      geom_boxplot() + theme_classic() +
      ggtitle(paste0("Res: ", currRes_v), subtitle = paste0("Mean sil.width: ", round(mean(currSil$width), digits = 4)))
    
    ### Based on silhouette width, change cluster accordingly
    currBestChoice_v <- ifelse(currSil$width > 0, 
                               as.numeric(as.character(currClusterCalls_v)), 
                               as.numeric(as.character(currSil$other)))
    #currBCOrig_v <- ifelse(currSil$width > 0, currClusterCalls_v, currSil$other)
    
    ### Table of changes
    compareCalls_tab <- table(Assigned=currClusterCalls_v,
                              Closest=currBestChoice_v)
    
    ### Calculate root mean-squared deviation for each cluster (measure of dispersion)
    ### High RMSD -> large internal heterogeneity -> might need further subclustering
    currRMSD <- bluster::clusterRMSD(embedding_v, currClusterCalls_v)
    currRMSD_df <- data.frame(Cluster = names(currRMSD), RMSD = currRMSD)
    
    ### Barplot
    bar_ls[[i]] <- ggplot(data = currRMSD_df, aes(x=Cluster, y = RMSD)) + 
      geom_bar(stat = 'identity') +
      theme_classic() + labs(x = "Cluster", y = "RMSD") +
      ggtitle(paste0("Res: ", currRes_v), subtitle = paste0("Mean RMSD: ", mean(currRMSD)))
    
    ### Output table
    currQC <- dplyr::tibble(res = currRes_v,
                     nClust = length(unique(currSil$cluster)),
                     meanSilWidth = mean(currSil$width),
                     meanRMSD = mean(currRMSD))
    
    if (i == 1) { qcFinal <- currQC } else { qcFinal <- rbind(qcFinal, currQC)}
    
  } # for i
  
  ### Return
  out_ls <- list("QC" = qcFinal,
                 "dimPlot" = dim_ls,
                 "silBox" = box_ls,
                 "RMSDBar" = bar_ls)
  
  return(out_ls)
  
} # clusterQC
