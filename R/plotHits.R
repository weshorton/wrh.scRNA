plotHits <- function(obj, deg_lsdt, genes_v, pop_v, name_v) {
  #' Plot Hits
  #' @description Make volcano, violin, and feature plots
  #' @param obj seurat object
  #' @param deg_lsdt list of deg results tables
  #' @param genes_v vector of genes to display
  #' @param pop_v name of the population (e.g. "TAMs" or "b12 TAMs"). Not used for searching so can be anything.
  #' @param name_v name of genes_v for the plot
  #' @details
  #' The names of deg_lsdt should be in format "Treat1_Treat2", the content of each list element
  #' is the specific DEG table associated with Treat1 vs Treat2 comparison.  
  #' obj should be subset to have the same cells that are tested in the Treat1 vs Treat2 comparison.
  #' genes_v is a vector of one or more genes that will be displayed in the plots.
  #' name_v is a name that describes the gene(s) in genes_v
  #' @return list of lists. First list level is the plot type (volcano, violin, feature), while the second level
  #' comprises all of the comparisons provided in deg_lsdt
  #' @export
  
  ###
  ### Volcano plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  volcano_lsgg <- sapply(names(deg_lsdt), function(x) {
    wrh.scRNA::plotVolcano(data_dt = deg_lsdt[[x]], ident1_v = gsub("_.*", "", x),
                           title_v = paste0(name_v, " Genes on\n", pop_v, " - ", gsub("_", " vs. ", x), " DEG"),
                           labelGenes_v = genes_v, labelSize_v = 5)
  }, simplify = F, USE.NAMES = T)
  
  ###
  ### Violin plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  violin_lsgg <- VlnPlot(object = obj, features = genes_v, group.by = "Treatment", combine = F)
  violin_lsgg <- sapply(violin_lsgg, function(x) {
    x + scale_fill_manual(values = treatColors_v, breaks = names(treatColors_v)) +
      big_label() + angle_x()}, simplify = F, USE.NAMES = T)
  violin_gg <- ggpubr::ggarrange(plotlist = violin_lsgg, ncol = 3, nrow = ceiling(length(violin_lsgg)/3),
                                 common.legend = T, legend = "bottom")
  violin_gg <- ggpubr::annotate_figure(violin_gg, 
                                       top = ggpubr::text_grob(paste0("Expression of\n", name_v, " Genes on\n", pop_v), 
                                                               size = 28, face = "bold"))
  
  ###
  ### Feature Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  feature_lsgg <- FeaturePlot(object = obj, features = genes_v, combine = F)
  feature_gg <- ggpubr::ggarrange(plotlist = feature_lsgg, ncol = 2, nrow = ceiling(length(feature_lsgg)/2))
  feature_gg <- ggpubr::annotate_figure(feature_gg, top = paste0("Expression of\n", name_v, " Genes on\n", pop_v))
  
  titles_v <- levels(obj$Treatment)
  splitFeature_lsgg <- sapply(genes_v, function(x) {
    if (!x %in% rownames(obj@assays$RNA) &
        !x %in% rownames(obj@assays$SCT)) return(NULL)
    plots_lsgg <- FeaturePlot(object = obj, features = x, split.by = "Treatment", combine = F)
    names(plots_lsgg) <- titles_v
    plots_lsgg <- sapply(names(plots_lsgg), function(y) {
      z <- plots_lsgg[[y]] + umapFigureTheme() + ggtitle(y)}, simplify = F, USE.NAMES = T)
    plot_gg <- ggpubr::ggarrange(plotlist = plots_lsgg, ncol = 2, nrow = ceiling(length(plots_lsgg)/2))
    plot_gg <- ggpubr::annotate_figure(plot_gg, 
                                       top = ggpubr::text_grob(paste0("Expression of\n", name_v, " Gene ", x, " on\n", pop_v, " by Treatment"), 
                                                               size = 28, face = "bold"))
    return(plot_gg)}, simplify = F, USE.NAMES = T)
  splitFeature_lsgg <- splitFeature_lsgg[which(sapply(splitFeature_lsgg, length) > 0)]
  
  
  ###
  ### Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  out_ls <- list("volcano" = volcano_lsgg,
                 "violin" = violin_gg,
                 "allFeature" = feature_gg,
                 "splitFeature" = splitFeature_lsgg)
  return(out_ls)
  
} # plotHits