labelSigHeatGenes <- function(seurat_obj, genes_lsv, assay_v = "fullSCT") {
  #' Label genes for heatmap
  #' @description Add a * to significant genes for heatmap output
  #' @param seurat_obj seurat object subset for the cells to plot in heatmap
  #' @param genes_lsv list of 2 vectors, 'heat' are the genes to label, 'label' are those to add * to (only done if heat != label)
  #' @param assay_v assay to use. Should usually be fullSCT
  #' @return seurat object with modified labels for genes and also vector of modified genes for doheatmap
  #' @export
  
  ### Check equality of supplied genes (TRUE if equal, character if not)
  allEqual_v <- all.equal(genes_lsv$heat, genes_lsv$label)
  
  if (class(allEqual_v) == "logical") {
    
    out_ls <- list("seurat" = seurat_obj, "genes" = genes_lsv$heat)
    
  } else {
    
    ### Get genes
    allGenes_v <- seurat_obj@assays[[assay_v]]@counts@Dimnames[[1]]
    
    ### Add asterisk for replacement and for output
    labeledGenes_v <- sapply(allGenes_v, function(x) ifelse(x %in% genes_lsv$label, paste0(x, "*"), x))
    outGenes_v <- sapply(genes_lsv$heat, function(x) ifelse(x %in% genes_lsv$label, paste0(x, "*"), x))
    
    ### Add to each slot (not sure how to loop this)
    ### Counts
    if (!all.equal(seurat_obj@assays[[assay_v]]@counts@Dimnames[[1]], allGenes_v)) {
      stop("Count slot doesn't match")
    } else {
      seurat_obj@assays[[assay_v]]@counts@Dimnames[[1]] <- labeledGenes_v
    } # fi
    
    ### data
    if (!all.equal(seurat_obj@assays[[assay_v]]@data@Dimnames[[1]], allGenes_v)) {
      stop("Data slot doesn't match")
    } else {
      seurat_obj@assays[[assay_v]]@data@Dimnames[[1]] <- labeledGenes_v
    } # fi
    
    ### scale.data
    if (!all.equal(rownames(seurat_obj@assays[[assay_v]]@scale.data), allGenes_v)) {
      stop("Scale data slot doesn't match")
    } else {
      rownames(seurat_obj@assays[[assay_v]]@scale.data) <- labeledGenes_v
    } # fi
    
    out_ls <- list("seurat" = seurat_obj, "genes" = unname(outGenes_v))
    
  } #fi all.equal
  
  return(out_ls)
  
} # labelSigHeatGenes