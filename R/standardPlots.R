standardPlots <- function(seurat_obj, reduction_v, clustCol_v, res_v, name_v, pt.size_v = 0.5, featurePlots_v = T,
                          maxQ_v = "q90") {
  #' Standard Seurat Plots
  #' @description
    #' Make some quick plots for unintegrated, harmony, or other objects
  #' @param seurat_obj seurat object to use for plotting
  #' @param reduction_v name of reduction to use for UMAP - usually "umap" or "harmony" 
  #' @param clustCol_v column name of seurat clusters to plot
  #' @param res_v resolution used for clusters
  #' @param name_v name to add to plot titles
  #' @param pt.size_v argument to pt.size of DimPlot and FeaturePlot
  #' @param featurePlots_v logical indicating if set of feature plots should also be made.
  #' @param maxQ_v cut off for feature plots. Use this if the scale is such that you can't see differences.
  #' @return list of ggplot objects
  #' @export
  
  ### To Do - this requires ggpmisc package, can't remember how to do that correctly rn, so just calling it here.
  library(ggpmisc)
  
  ### Empty list to hold plots
  out_lsgg <- list()
  
  ### Default
  out_lsgg[["clusters"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = clustCol_v, label = T, pt.size = pt.size_v) +
    coord_equal() +
    ggtitle(paste0(name_v, " clusters; res - ", res_v))
  
  ### Treatment plots
  if (!"Treatment" %in% colnames(seurat_obj@meta.data)) {
    warning("Treatment column not found. Will not make those plots")
  } else {
    
    ### Color by treatment
    out_lsgg[["cTreat"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "Treatment", label = F, pt.size = pt.size_v) +
      coord_equal() +
      ggtitle(paste0(name_v, " clusters by Treatment"))
    
    ### Split by batch, color treatment
    out_lsgg[["fBatchcTreat"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "Treatment", 
                                          split.by = "batchID", label = F, pt.size = pt.size_v) +
      coord_equal() + theme(legend.position = "bottom") +
      ggtitle(paste0(name_v, " clusters; res - ", res_v))
    
    ### Split by treatment
    out_lsgg[["fTreat"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = clustCol_v, 
                                    split.by = "Treatment", label = F, pt.size = pt.size_v, ncol = 2) +
      coord_equal() + theme(legend.position = "bottom") +
      ggtitle(paste0(name_v, " clusters; res - ", res_v))
  } # fi treat
  
  ### Batch plots
  if (!"batchID" %in% colnames(seurat_obj@meta.data)) {
    warning("batchID column not found. Will not make those plots.")
  } else {
    
    ### Color by batch
    out_lsgg[["cBatch"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "batchID", label = F, pt.size = pt.size_v) +
      coord_equal() +
      ggtitle(paste0(name_v, " clusters by Batch"))
    
    ### Split by batch
    out_lsgg[["fBatch"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = clustCol_v, 
                                    split.by = "batchID", label = F, pt.size = pt.size_v) +
      coord_equal() + theme(legend.position = "bottom") +
      ggtitle(paste0(name_v, " clusters; res - ", res_v))
  } # fi batch
  
  ### TCR/BCR
  if (!"tbcr" %in% colnames(seurat_obj@meta.data)) {
    warning("tbcr column not found. Will not make plot.")
  } else {
    out_lsgg[["tbcr"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "tbcr", cols = c("red", "grey", "blue"), 
                                  label = F, pt.size = pt.size_v) + coord_equal() +
      ggtitle(paste0(name_v, " clusters - TCR/BCR"))
  }
  
  ### Color by original populations
  if (!"ind_mPop" %in% colnames(seurat_obj@meta.data)) {
    warning("ind_mPop column not found. Will not make plot.")
  } else {
    out_lsgg[["origMpop"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "ind_mPop", label = F, pt.size = pt.size_v) +
      coord_equal() + ggtitle(paste0(name_v, "\nindividual-batch Cell Classes"))
  }
  
  ### Color by full data populations
  if (!"full_mPop" %in% colnames(seurat_obj@meta.data)) {
    warning("full_mPop column not found. Will not make plot.")
  } else {
    out_lsgg[["fullMpop"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "full_mPop", label = F, pt.size = pt.size_v) +
      coord_equal() + ggtitle(paste0(name_v, "\nFull Dataset Cell Classes"))
  }
  
  ### Color by major populations
  if (!"mPop" %in% colnames(seurat_obj@meta.data)) {
    warning("mPop column not found. Will not make plot.")
  } else {
    out_lsgg[["mPop"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "mPop", label = F, pt.size = pt.size_v) +
      coord_equal() + ggtitle(paste0(name_v, "\nMajor Population Cell Classes"))
  }
  
  ### Color by sub-populations 
  if (!"sPop" %in% colnames(seurat_obj@meta.data)) {
    warning("sPop column not found. Will not make plot.")
  } else {
    out_lsgg[["sPop"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "sPop", label = F, pt.size = pt.size_v) +
      coord_equal() + ggtitle(paste0(name_v, "\nSub Population Cell Classes"))
  }
  
  ### Color by collapsed sub-populations, if present
  if ("collapsePop" %in% colnames(seurat_obj@meta.data)) {
    out_lsgg[["collapsePop"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "collapsePop", label = F, pt.size = pt.size_v) +
      coord_equal() + ggtitle(paste0(name_v, "\nCollapse Population Cell Classes"))
  }
  
  if ("full_collapsePop" %in% colnames(seurat_obj@meta.data)) {
    out_lsgg[["full_collapsePop"]] <- Seurat::DimPlot(seurat_obj, reduction = reduction_v, group.by = "full_collapsePop", label = F, pt.size = pt.size_v) +
      coord_equal() + ggtitle(paste0(name_v, "\nFull Data Collapse Population Cell Classes"))
  }
  
  if (featurePlots_v) {
  
    ### Cell cycle Feature plot - S phase
    if (!"S.Score" %in% colnames(seurat_obj@meta.data)) {
      warning("S.Score column not found. Will not make plot.")
    } else {
      out_lsgg[["S.Score"]] <- Seurat::FeaturePlot(seurat_obj, reduction = reduction_v, features = "S.Score", max.cutoff = maxQ_v, pt.size = pt.size_v) +
        coord_equal() + ggtitle(paste0(name_v, " S Phase Scores")) +
        scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
    } # fi
    
    ### Cell cycle Feature plot - G2M phase
    if (!"G2M.Score" %in% colnames(seurat_obj@meta.data)) {
      warning("G2M.Score column not found. Will not make plot.")
    } else {
      out_lsgg[["G2M.Score"]] <- Seurat::FeaturePlot(seurat_obj, reduction = reduction_v, features = "G2M.Score", max.cutoff = maxQ_v, pt.size = pt.size_v) +
        coord_equal() + ggtitle(paste0(name_v, " G2M Phase Scores")) +
        scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
    } # fi
    
    ### Feature plot - log UMI
    if (!"logUMI" %in% colnames(seurat_obj@meta.data)) {
      warning("logUMI column not found. Will not make plot.")
    } else {
      out_lsgg[["logUMI"]] <- Seurat::FeaturePlot(seurat_obj, reduction = reduction_v, features = "logUMI", pt.size = pt.size_v) +
        coord_equal() + ggtitle(paste0(name_v, " Log nUMI")) +
        scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
    } # fi
    
    ### Feature plot - log UMI
    if (!"logGene" %in% colnames(seurat_obj@meta.data)) {
      warning("logGene column not found. Will not make plot.")
    } else {
      out_lsgg[["logGene"]] <- Seurat::FeaturePlot(seurat_obj, reduction = reduction_v, features = "logGene", pt.size = pt.size_v) +
        coord_equal() + ggtitle(paste0(name_v, " Log nGene")) +
        scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
    } # fi
    
    ### Feature plot - percent mito
    if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
      warning("percent.mt column not found. Will not make plot.")
    } else {
      out_lsgg[["percent.mt"]] <- Seurat::FeaturePlot(seurat_obj, reduction = reduction_v, features = "percent.mt", pt.size = pt.size_v) +
        coord_equal() + ggtitle(paste0(name_v, " Percent Mt")) +
        scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
    } # fi
    
  } # fi featurePlots_v
  
  ### Add text grob of parameters
  param_dt <- data.table("Parameter" = c("Reduction", "Cluster Column", "Resolution", "Name"),
                         "Value" = c(reduction_v, clustCol_v, res_v, name_v))
  param_gg <- ggplot2::ggplot() + ggplot2::annotate(geom = "table", x = 0, y = 0, label = list(param_dt)) + theme_void()
  out_lsgg[["params"]] <- param_gg
  
  ### Output
  return(out_lsgg)
}