neoantigenPlots <- function(obj, deg_dt, out_lsls, cells_v, ident1_v = "neo.c8", ident2_v, cellCounts_v) {
  #' Neoantigen Plots
  #' @description Dot and Volcano for specific genes
  #' @param obj corresponding seurat object
  #' @param deg_dt Find markers output fed into wrangleDEG
  #' @param out_lsls output of wrangleDEG
  #' @param cells_v character vector indicating cells included in comparison (e.g. all or 4x)
  #' @param ident1_v name of first identity in comparison. Used for titles and cell counts
  #' @param ident2_v name of second identity. Used for titles and cell counts.
  #' @param cellCounts_v vector of length 2 that indicates number of cells in each identity. Must be in order ident1, ident2
  #' @return list containing one dot plot and one volcano plot for each entry
  #' @export
  
  genes_lsv <- c(out_lsls$genes$up, out_lsls$genes$dn)
  bubble_lsgg <- volcano_lsgg <- list()
  
  for (i in 1:length(genes_lsv)) {
    
    currName_v <- names(genes_lsv)[i]
    currOutName_v <- gsub("\n", " ", currName_v)
    currGenes_v <- genes_lsv[[currName_v]]
    
    if (length(currGenes_v) == 0) {
      
      out_lslsgg <- list("bubble" = list(), "volcano" = list())
      
    } else {
      
      ### Make bubble title
      bubbleTitle_v <- paste0(cells_v, " cells displayed\nGenes are ", length(currGenes_v), " ", currName_v,
                              " significant in\n", cellCounts_v[1], " ", ident1_v, " cells vs. ", 
                              cellCounts_v[2], " ", ident2_v, " cells")
      volcanoTitle_v <- paste0("DEG using ", cells_v, " cells displayed\nGenes are ", length(currGenes_v), " ", currName_v,
                               " significant in\n", cellCounts_v[1], " ", ident1_v, " cells vs. ", 
                               cellCounts_v[2], " ", ident2_v, " cells")
      
      bubble_lsgg[[currOutName_v]] <- DotPlot(object = obj,
                                              assay = "fullSCT",
                                              features = currGenes_v,
                                              cols = "RdYlBu",
                                              dot.scale = 8) +
        coord_flip() +
        #big_label() +
        wrh.rUtils::angle_x() +
        theme(line = element_line(linewidth = 0)) +
        theme(rect = element_rect(linewidth = 0.0)) +
        theme(panel.border = element_rect(linewidth = 0.0)) +
        theme(axis.line = element_line(linewidth = 0.0)) +
        theme(legend.position="bottom", legend.box = "vertical") +
        theme(plot.title = element_text(hjust = 0.5)) +
        ggtitle(bubbleTitle_v)
      #ggtitle(paste0(cells_v, " cells '", comparison_v, "' ", currName_v," (", length(currGenes_v), ")"))
      
      volcano_lsgg[[currOutName_v]] <- plotVolcano(data_dt = deg_dt, ident1_v = ident1_v, labelGenes_v = currGenes_v, labelSize_v = 3,
                                                   title_v = volcanoTitle_v)
    } # for i
    
    out_lslsgg <- list("bubble" = bubble_lsgg, "volcano" = volcano_lsgg)
    
  } # fi
  
  return(out_lslsgg)
  
} # neoantigenPlots