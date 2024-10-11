checkAtlasAnnotations <- function(clusterMarkers_dt, popMap_dt, popMarkerChecks_v, topGenes_v = 15) {
  #' Check Atlas Annotation
  #' @description
  #' Check each seurat cluster and cell type identity in Eric's public annotations
  #' against the genes we have.
  #' @param clusterMarkers_dt data.table containing FindAllMarkers results.
  #' @param popMap_dt data.table mapping seurat clusters to cell types
  #' @param popMarkerChecks_v named vector of genes associated with various populations
  #' @param topGenes_v numeric indicating how many of the top genes to label in the volcano (fed to plotVolcano() labelTop_v argumennt)
  #' @return list of volcano plots of hits and data.table with intersections
  #' @export
  
  clusters_v <- popMap_dt$Cluster
  volcano_lsgg <- intersect_lsv <- list()
  popRes_v <- character()
  
  for (i in 1:length(clusters_v)) {
    
    ### Get data
    currCluster_v <- clusters_v[i]
    currPop_v <- popMap_dt[Cluster == currCluster_v, "Population"]
    currClusterMarkers_dt <- clusterMarkers_dt[cluster == currCluster_v,]
    colnames(currClusterMarkers_dt)[colnames(currClusterMarkers_dt) == "gene"] <- "Gene"
    setorder(currClusterMarkers_dt, -avg_log2FC)
    
    ### Get intersection
    currIntersect_v <- popMarkerChecks_v[popMarkerChecks_v %in% currClusterMarkers_dt$Gene]
    intersect_lsv[[paste0("c", currCluster_v)]] <- c(currPop_v, currIntersect_v)
    
    ### Summarize intersection
    if (length(currIntersect_v) > 0) {
      currIntNames_dt <- as.data.table(table(names(currIntersect_v)))
      currIntNames_dt[,Tot := sapply(V1, function(x) length(grep(x, names(popMarkerChecks_v))))]
      currIntNames_dt[,Pct := round(N/Tot*100,digits=2)]
      setorderv(currIntNames_dt, c("Pct", "Tot"), order = -1)
      currIntRes_v <- paste(setdiff(currIntNames_dt$V1[1:2], NA), collapse = ", ")
      currIntTitle_v <- apply(currIntNames_dt, 1, function(x) paste0("\nFound ", x[2], " of ", x[3], " ", x[1], " Genes"))
    } else {
      currIntTitle_v <- "No Pop Markers Found"
      currIntRes_v <- "none"
    } # fi
    
    popRes_v <- c(popRes_v, currIntRes_v)
    
    ### Make violin
    currVolcano_gg <- plotVolcano(data_dt = currClusterMarkers_dt, geneCol_v = "Gene", ident1_v = paste0("Cluster ", currCluster_v),
                                  labelGenes_v = currIntersect_v, labelTop_v = topGenes_v, labelSize_v = 4, force_v = 40,
                                  title_v = paste0("Check Genes in Cluster ", currCluster_v, "\nCalled ", currPop_v, paste(currIntTitle_v, collapse = "")))
    
    volcano_lsgg[[paste0("c", currCluster_v)]] <- currVolcano_gg
    
  } # for i
  
  intersect_dt <- listToDT(intersect_lsv)
  popRes_dt <- data.table("Cluster" = clusters_v,
                          "Found Pops" = popRes_v)
  popRes_dt <- merge(popMap_dt, popRes_dt, by = "Cluster")
  
  return(list("volcano" = volcano_lsgg, "intersections" = intersect_dt, "popRes" = popRes_dt))
  
} # checkAtlasAnnotation