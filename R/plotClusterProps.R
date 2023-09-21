plotClusterProps = function (seurat_obj, clusters_v = "seurat_clusters", ident_v, legendName_v) { 
  #' Plot Cluster Proportions
  #' @description Plot proportions of idents falling in clusters. Adapted from
  #' https://www.singlecellcourse.org/scrna-seq-dataset-integration.html#seurat-v3-3-vs-5-10k-pbmc
  #' @param seurat_obj - seurat object with at least one "seurat_clusters" column
  #' @param clusters_v - the name of the "seurat_clusters" column if not default name.
  #' @param ident_v - the name of the metadata column to use as grouping variable
  #' @param legendName_v name to put in legend to label ident_v
  #' @return Returns patchwork of two ggplots of cluster data plots
  #' @export
  #' 
  
  ### Tabulate occurrences and melt for plotting
  counts_tab <- table(seurat_obj@meta.data[[clusters_v]], seurat_obj@meta.data[[ident_v]])
  counts_dt <- wrh.rUtils::convertDFT(as.data.frame.matrix(counts_tab), newName_v = "cluster")
  meltCounts_dt <- reshape2::melt(counts_dt, id.vars = "cluster")
  meltCounts_dt$cluster <- as.factor(meltCounts_dt$cluster)
  
  ### Grab size of each cluster
  clustSize_df   <- aggregate(value ~ cluster, data = meltCounts_dt, FUN = sum)
  
  ### Sort clusters by size (should also be cluster name...)
  sortedLabels_v <- paste(sort(as.integer(clustSize_df$cluster),decreasing = T))
  clustSize_df$cluster <- factor(clustSize_df$cluster,levels = sortedLabels_v)
  meltCounts_dt$cluster <- factor(meltCounts_dt$cluster,levels = sortedLabels_v)
  
  colnames(meltCounts_dt)[2] <- legendName_v
  
  ### Bar chart of cluster sizes
  p1 <- ggplot(clustSize_df, aes(y= cluster,x = value)) + 
    geom_bar(position="dodge", stat="identity",fill = "grey60") + 
    theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
  
  ### Bar chart of ident_v proportions
  p2 <- ggplot(meltCounts_dt,aes(x=cluster,y=value,fill=!!sym(legendName_v))) + 
    geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
    scale_fill_brewer(palette = "Set2") +
    ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")
  
  ### Combine
  out_gg <- p2 + p1 + patchwork::plot_layout(widths = c(3,1))
  
  ### Return
  return(out_gg)
  
} # plotClusterProps
