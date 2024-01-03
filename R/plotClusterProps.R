plotClusterProps = function (seurat_obj, clusters_v = "seurat_clusters", ident_v = NULL, legendName_v, clustColors_v, identColors_v = treatColors_v, sortBySize_v = F) { 
  #' Plot Cluster Proportions
  #' @description Plot proportions of idents falling in clusters. Adapted from
  #' https://www.singlecellcourse.org/scrna-seq-dataset-integration.html#seurat-v3-3-vs-5-10k-pbmc
  #' @param seurat_obj - seurat object with at least one "seurat_clusters" column
  #' @param clusters_v - the name of the "seurat_clusters" column if not default name.
  #' @param ident_v - the name of the metadata column to use as grouping variable
  #' @param legendName_v name to put in legend to label ident_v
  #' @param clustColors_v named vector of colors that correspond to clusters_v column
  #' @param identColors_v named vector of colors that correspond to ident_v
  #' @param sortBySize_v logical, will sort by cluster size rather than factor levels
  #' @return Returns patchwork of two ggplots of cluster data plots
  #' @export
  #' 
  
  ### Tabulate occurrences and melt for plotting
  counts_tab <- table(seurat_obj@meta.data[[clusters_v]], seurat_obj@meta.data[[ident_v]])
  counts_dt <- wrh.rUtils::convertDFT(as.data.frame.matrix(counts_tab), newName_v = "cluster")
  meltCounts_dt <- reshape2::melt(counts_dt, id.vars = "cluster")
  
  ### Sort
  if (length(grep("neo", meltCounts_dt$cluster)) > 0) {
    temp_v <- unique(gsub("neo\\.", "", meltCounts_dt$cluster))
    meltCounts_dt$cluster <- factor(meltCounts_dt$cluster, levels = paste0("neo.", temp_v[order(as.numeric(gsub("[a-z]*", "", temp_v)))]))
  } else {
    if (!is.null(levels(seurat_obj@meta.data[[clusters_v]]))) {
      meltCounts_dt$cluster <- factor(meltCounts_dt$cluster, levels = levels(seurat_obj@meta.data[[clusters_v]]))
    } else {
      meltCounts_dt$cluster <- factor(meltCounts_dt$cluster)
    }
  }
  
  ### Grab size of each cluster
  clustSize_df   <- aggregate(value ~ cluster, data = meltCounts_dt, FUN = sum)
  
  ### Sort clusters by size (should also be cluster name...)
  if (sortBySize_v) {
    sortedLabels_v <- paste(sort(as.integer(clustSize_df$cluster),decreasing = T))
    clustSize_df$cluster <- factor(clustSize_df$cluster,levels = sortedLabels_v)
    meltCounts_dt$cluster <- factor(meltCounts_dt$cluster,levels = sortedLabels_v)
  } else {
    clustSize_df$cluster <- factor(clustSize_df$cluster, levels = levels(meltCounts_dt$cluster))
  }
  
  colnames(meltCounts_dt)[2] <- legendName_v
  
  ### Bar chart of cluster sizes
  p1 <- ggplot(clustSize_df, aes(y= cluster,x = value)) + 
    geom_bar(position="dodge", stat="identity", aes(fill = cluster)) + 
    scale_fill_manual(values = clustColors_v, breaks = names(clustColors_v)) +
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("") +
    ggtitle("Overall Cluster Size")
  
  ### Bar chart of ident_v proportions
  p2 <- ggplot(meltCounts_dt,aes(x=cluster,y=value,fill=!!sym(legendName_v))) + 
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
    scale_fill_manual(values = identColors_v, breaks = names(identColors_v)) +
    ylab("Fraction of cells") + xlab("Cluster") + theme(legend.position="bottom") +
    ggtitle(paste0("Proportion of cells from each ", ident_v, " per cluster"))
  
  p3 <- ggplot(meltCounts_dt, aes(x=!!sym(legendName_v), y = value, fill = cluster)) +
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
    scale_fill_manual(values = clustColors_v, breaks = names(clustColors_v)) +
    ylab("Fraction of cells in each dataset") + xlab("Cluster") + theme(legend.position="bottom") +
    ggtitle(paste0("Proportion of cells from each cluster per ", ident_v))
  
  ### Combine
  out_gg <- ggpubr::ggarrange(plotlist = list(p1, p2, p3), ncol = 3, nrow = 1)
  out_gg <- ggpubr::annotate_figure(out_gg, top = ggpubr::text_grob(label = paste0("Cluster x ", ident_v, " Summary"), size = "20"))
  
  ### Return
  return(out_gg)
  
} # plotClusterProps
