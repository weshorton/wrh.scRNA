### Silhouette scoring function

silhouetteScore <- function(seurat_obj, resolution_v = 0.6, graph_v = "RNA_snn", ndims_v = 50, reduction_v = "pca", umap_v = "umap") {
  #' Score Silhouettes of Clusters
  #' @description Given a seurat object and desired parameters, calculate clusters and determine silhouette scores for clusters.
  #' @param seurat_obj A seurat object
  #' @param resolution_v Resolution to use
  #' @param graph_v name of graph to use.
  #' @param ndims_v The number of dimensions used to run
  #' @param reduction_v The reduction to use (pca, harmony, or iNMF, usually)
  #' @return Not sure yet
  #' @export
  
  ### Calculate clusters
  seurat_obj <- Seurat::FindClusters(object = seurat_obj,
                             resolution = resolution_v,
                             graph.name = graph_v)
  
  ### DimPlot
  p1_gg <- Seurat::DimPlot(seurat_obj,
                   reduction = umap_v,
                   group.by = "seurat_clusters") +
    coord_equal() +
    ggtitle(paste0("Clusters for res: ", resolution_v, "\nUsing graph: ", graph_v, "\nOn Reduction: ", reduction_v))
  
  ### Grab distance matrix - takes forever! Try future?
  redCols_v <- ncol(seurat_obj@reductions[[reduction_v]]@cell.embeddings)
  
  if (redCols_v >= ndims_v) {
    distance_mat <- dist(seurat_obj@reductions[[reduction_v]]@cell.embeddings[,1:ndims_v])
  } else {
    print("More dimensions requested than available.")
  }
  
  ### Get Clusters and silhouette
  clusters_v <- seurat_obj@meta_data$seurat_clusters
  nClust_v <- length(levels(clusters_v))
  silhouette_v <- cluster::silhouette(as.numeric(clusters_v), dist = distance_mat)
  
  ### Add scores to metadata and get average
  seurat_obj@meta.data$silhouette_score <- silhouette[,3]
  meanSilhouetteScore_v <- mean(seurat_object@meta.data$silhouette_score)
  
  ### Plot
  plot_mat <- seurat_object@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    arrange(seurat_clusters,-silhouette_score) %>%
    mutate(barcode = factor(barcode, levels = barcode))
  
  p2_gg <- ggplot(data = plot_mat) +
    geom_col(aes(x = barcode, y = silhouette_score, fill = seurat_clusters), show.legend = F) +
    geom_hline(yintercept = meanSilhouetteScore_v, color = "red", linetype = "dashed") +
    scale_x_discrete(name = "Cells") + scale_y_continuous(name = "Silhouette Score") +
    #labs(x = "Cells", y = "Silhouette Score") +
    myTheme() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggtitle(paste0("Clusters for res: ", resolution_v, " Using graph: ", graph_v, "\nOn Reduction: ", reduction_v,
                   " Mean Sil: ", meanSilhouetteScore_v, " nClust: ", nClust_v))
  
  
  ### Output
  

    
  
}




