makeDotPlots <- function(seurat_obj, name_v, assay_v = "fullSCT", features_v, cols_v = "RdYlBu", 
                         groupBy_v, splitBy_v = NULL) {
  #' Make Dot Plots
  #' @description standard seurat dotplot wrapper
  #' @param seurat_obj object to use for plot
  #' @param name_v name used for title
  #' @param assay_v assay used for plot
  #' @param features_v what genes to plot
  #' @param cols_v dotplot colors
  #' @param groupBy_v fed into the group.by argument, usually cluster or treatment
  #' @param splitBy_v optional split.
  #' @return ggplot object
  #' @export
  
  ## Double check features
  features_v <- intersect(features_v, rownames(seurat_obj@assays[[assay_v]]@scale.data))
  
  ## Assign ident
  Idents(seurat_obj) <- groupBy_v
  
  if (is.null(splitBy_v)) {
    
    ## Build plot
    plot_gg <- DotPlot(seurat_obj,
                       features = features_v,
                       cols = cols_v,
                       group.by = groupBy_v,
                       dot.scale = 8)
    
    ## Add to plot
    plot_gg <- plot_gg +
      ggtitle(name_v) + angle_x() + coord_flip() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_line(colour = "gray",linetype = "dashed", linewidth=0.35), 
            panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) + 
      theme(legend.position="bottom")
    
    plots_lsgg <- plot_gg
    
  } else {
    
    splits_v <- levels(seurat_obj@meta.data[[splitBy_v]])
    plots_lsgg <- list()
    
    for (i in 1:length(splits_v)) {
      
      currSplit_v <- splits_v[i]
      currCells_v <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[splitBy_v]] == currSplit_v,])
      currObj <- subset(seurat_obj, cells = currCells_v)
      Idents(currObj) <- splitBy_v
      
      ## Build plot
      plot_gg <- DotPlot(currObj,
                         features = features_v,
                         cols = cols_v,
                         group.by = groupBy_v,
                         dot.scale = 8)
      
      ## Add to plot
      plot_gg <- plot_gg +
        ggtitle(name_v) + angle_x() + 
        coord_flip() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(colour = "gray",linetype = "dashed", linewidth=0.35), 
              panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) + 
        theme(legend.position="bottom")
      
      plots_lsgg[[currSplit_v]] <- plot_gg
      
    } # for i
    
  } # fi
  
  return(plots_lsgg)
  
} # makeDotPlots