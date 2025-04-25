myTreatFacetUmap <- function(obj, umap_v, vals_v, colors_v, title_v = NULL, size_v = 36, 
                             facetVar_v = "Treatment", groupVar_v = "sPop", plotType_v = "umap") {
  #' My Treatment Facet Umap
  #' @description
    #' custom facet umap for figures
  #' @param obj seurat object
  #' @param umap_v reduction to use
  #' @param vals_v the unique values of facetVar_v (need to be in order)
  #' @param colors_v colors to plot with
  #' @param title_v optional plot title
  #' @param size_v size for plot title
  #' @param facetVar_v variable to facet by (should be 'Treatment')
  #' @param plotType_v either 'umap' or 'scatter'
  #' @export
  
  ### Pattern list for scatterhatch
  patternList_ls <- list()
  patterns_v <- c("-", "|", "/", "\\", "x", "+", "", "-", "|", "/", "\\", "x", "+")
  for (i in 1:length(patterns_v)) {
    patternList_ls[[i]] <- list("pattern" = patterns_v[i], "lineAlpha" = 0.5)
  }
  
  ### Get cells for each value of facetVar
  cells_lsv <- sapply(vals_v, function(x) {
    rownames(obj@meta.data[obj@meta.data[[facetVar_v]] == x,])
    }, simplify = F, USE.NAMES = T)
  
  ### Get UMAP embeddings (for scatter only)
  embeddings_df <- obj@reductions[[umap_v]]@cell.embeddings
  embeddings_df <- merge(embeddings_df, obj@meta.data[,c(groupVar_v, facetVar_v),drop=F], by = 0, sort = F)
  colnames(embeddings_df) <- c("Cell", "X", "Y", groupVar_v, facetVar_v)
  
  ### Get axes limits
  lims_lsv <- list("X" = c(min(embeddings_df[,"X"])*1.01, max(embeddings_df[,"X"])*1.01),
                   "Y" = c(min(embeddings_df[,"Y"])*1.01, max(embeddings_df[,"Y"])*1.01))
  
  ### Get plotting differentiators
  if (length(vals_v) %% 2 != 0) {
    xAxis_v <- rev(rev(vals_v)[1])
    rightCol_v <- vals_v[seq(2, length(vals_v), 2)]
    leftCol_v <- vals_v[seq(1, length(vals_v)-1, 2)]
    keepTitle_v <- leftCol_v[2]
  } else {
    xAxis_v <- rev(rev(vals_v)[1:2])
    rightCol_v <- vals_v[seq(2, length(vals_v), 2)]
    leftCol_v <- vals_v[seq(1, length(vals_v)-1, 2)]
    keepTitle_v <- leftCol_v[2]
  } # fi
  
  ### Go through each value of facetVar and make a umap
  ### Then, depending on its position on the graph, change
  ### its formatting
  outPlot_lsgg <- list()
  for (varName_v in names(cells_lsv)) {
    
    ### Handle plot type
    if (plotType_v == "umap") {
      
      ### Base Plot
      p_gg <- DimPlot(object = obj, cells = cells_lsv[[varName_v]], 
                      reduction = umap_v, label = F, 
                      group.by = groupVar_v, cols = colors_v)
      
    } else if (plotType_v == "scatter") {
      
      ### Embeddings for plot
      currEmbeddings_df <- embeddings_df[embeddings_df[[facetVar_v]] == varName_v,]
      names(patternList_ls) <- levels(currEmbeddings_df[[groupVar_v]])
      
      ### Handle missing group values
      currGroups_v <- unique(currEmbeddings_df[[groupVar_v]])
      currEmbeddings_df[[groupVar_v]] <- factor(as.character(currEmbeddings_df[[groupVar_v]]), levels = intersect(levels(currGroups_v), currGroups_v))
      currPatternList_ls <- patternList_ls[levels(currEmbeddings_df[[groupVar_v]])]
      currColors_v <- colors_v[intersect(names(colors_v), levels(currEmbeddings_df[[groupVar_v]]))]
      
      p_gg <- scatterHatch(data = currEmbeddings_df,
                           x = "X", y = "Y",
                           color_by = groupVar_v,
                           pointAlpha = 1,
                           patternList = currPatternList_ls,
                           colorPalette = currColors_v)
      
    } else {
      
      stop(sprintf("Only 'umap' and 'scatter' are allowed for plotType_v. %s provided.\n", plotType_v))
      
    } # fi plot type
    
    p_gg <- p_gg +
      ggtitle(varName_v) +
      transparentTheme() +
      umapFigureTheme() + 
      coord_equal() +
      scale_y_continuous(limits = lims_lsv$Y) + 
      scale_x_continuous(limits = lims_lsv$X)
    
    ### X-axis
    if (varName_v %in% xAxis_v) {
      p_gg <- p_gg + theme(axis.line.x = element_line(linewidth = 0.25), 
                           axis.ticks.x = element_line(linewidth = 0.25))
    } else {
      p_gg <- p_gg + theme(axis.title.x = element_text(color = "white"),
                           axis.text.x = element_text(color = "white"),
                           axis.line.x = element_line(colour = "white"),
                           axis.ticks.x = element_line(color = "white"))
    } # fi xAxis_v
    
    ### Remove y-axis from right column
    if (varName_v %in% rightCol_v) {
      p_gg <- p_gg + theme(axis.title.y = element_text(color = "white"),
                           axis.text.y = element_text(color = "white"),
                           axis.line.y = element_line(colour = "white"),
                           axis.ticks.y = element_line(color = "white"))
    } # fi right side
    
    ### Remove axis title for all but one of left column; and make y-axis thinner for all
    if (varName_v %in% leftCol_v) {
      p_gg <- p_gg + theme(axis.line.y = element_line(linewidth = .25), 
                           axis.ticks.y = element_line(linewidth = .25))
      if (varName_v != keepTitle_v) {
        p_gg <- p_gg + theme(axis.title.y = element_text(color = "white"))
      } # fi 3xNR
    } # fi left side
    
    outPlot_lsgg[[varName_v]] <- p_gg
  } # for varName_v
  
  ### Combine into 1 plot
  facet_gg <- ggpubr::ggarrange(plotlist = outPlot_lsgg, ncol = 2, nrow = ceiling(length(outPlot_lsgg)/2), common.legend = T)
  
  ### Label
  if (!is.null(title_v)) {
    facet_gg <- ggpubr::annotate_figure(p = facet_gg, top = text_grob(label = title_v, size = size_v))
  } # fi title
  
} # myTreatFacetUmap