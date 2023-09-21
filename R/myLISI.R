myLISI <- function(seurat_obj, reduction_v, metaCols_v, name_v) {
  #' Calc and Plot LISI
  #' @description Wrapper function to run compute_lisi and make plots
  #' @param seurat_obj seurat object to use
  #' @param reduction_v name of reduction to use for embeddings
  #' @param metaCols_v vector of column names to extract from metadata. First must be the grouping variable
  #' (e.g. batchID); second must be the cell classifier variable
  #' @param name_v name to prepend to plots
  #' @details This still needs some work. Right now must have both meta columns (batch and cell class)
  #' @return list of ggplot objects
  #' @export
  
  ### Extract embeddings
  embeddings_df <- seurat_obj@reductions[[reduction_v]]@cell.embeddings
  
  ### Get metadata
  meta_df <- seurat_obj@meta.data[,metaCols_v]
  
  ### Calculate LISI
  lisi_df <- lisi::compute_lisi(X = embeddings_df, meta_data = meta_df, label_colnames = metaCols_v)
  colnames(lisi_df) <- c("iLISI", "cLISI")
  
  ### Merge results with metadata
  lisi_df <- merge(lisi_df, meta_df, by = 0, sort = F)
  
  ### iLISI
  iLISI_gg <- ggplot(data = lisi_df, aes(x = iLISI)) +
    geom_density() + my_theme() +
    ggtitle(paste0(name_v, " 'integration' LISI"))
  
  ### iLISI by population
  iLISI_byPop_gg <- ggplot(data = lisi_df, aes(x = iLISI, color = !!sym(metaCols_v[2]))) +
    geom_density() + my_theme() +
    ggtitle(paste0(name_v, " 'integration' LISI"))
  
  ### iLISI by batch
  iLISI_byBatch_gg <- ggplot(data = lisi_df, aes(x = iLISI, color = !!sym(metaCols_v[1]))) +
    geom_density() + my_theme() +
    ggtitle(paste0(name_v, " 'integration' LISI"))
  
  ### cLISI
  cLISI_gg <- ggplot(data = lisi_df, aes(x = cLISI)) +
    geom_density() + my_theme() +
    ggtitle(paste0(name_v, " 'accuracy' LISI"))
  
  ### Combine and output
  out_lsdfgg <- list("data" = lisi_df, "plots" = list("iLISI" = iLISI_gg, "iLISIv1" = iLISI_byBatch_gg, "iLISIv2" = iLISI_byPop_gg, "cLISI" = cLISI_gg))
  return(out_lsdfgg)
  
}
