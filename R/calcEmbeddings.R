calcEmbeddings <- function(seurat_obj, clusterName_v, reductionName_v, breaks_v = T, ...) {
  #' Calculate Feature Embeddings
  #' @description Calculate feature embeddings for a given single cell reduction
  #' @param seurat_obj seurat object with data reduction
  #' @param clusterName_v name of meta data column containing cluster identities to calculate embeddings for
  #' @param reductionName_v name of the reduction that embeddings are calculated from
  #' @param breaks_v logical. Should color breaks be adjusted or not
  #' @param ... extra parameters to pass to determineBreaks
  #' @return data.frame of embeddings and prints a heatmap
  #' @export
  
  ### Check that cluster exists
  if (!clusterName_v %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("Cluster name %s not found in seurat object.\n", clusterName_v))
  } # fi
  
  ### Check that reduction exists
  if (!reductionName_v %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction %s not found in seurat object.\n", reductionName_v))
  } # fi
  
  ### Calculate average embeddings per cluster
  clusterEmbedding_tib <- dplyr::tibble(seurat_clusters = seurat_obj@meta.data[[clusterName_v]]) %>%
    cbind(seurat_obj@reductions[[reductionName_v]]@cell.embeddings) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarize(dplyr::across(dplyr::everything(), list(mean)))
  
  ### Adjust color breaks
  if (breaks_v) {
    breaks_lsv <- wrh.rUtils::determineBreaks(data_mat = scale(data.frame(clusterEmbedding_tib[,2:ncol(clusterEmbedding_tib)])),
                                  colors_v = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
                                  ...)
    color_v <- breaks_lsv$colors
    break_v <- breaks_lsv$breaks
  } else {
    color_v <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
    break_v <- NA
  } # fi
  
  ### Heatmap
  print(pheatmap::pheatmap(mat = data.frame(clusterEmbedding_tib[,2:ncol(clusterEmbedding_tib)]),
                     row.names = clusterEmbedding_tib[[clusterName_v]],
                     scale = 'row',
                     color = color_v,
                     breaks = break_v,
                     main = paste0("Avg embeddings for ", clusterName_v, " in ", reductionName_v, " reduction\nrow z-score")))
  
  ### Output
  return(clusterEmbedding_tib)
}
