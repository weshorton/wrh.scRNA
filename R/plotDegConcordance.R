plotDegConcordance <- function(data1_dt, suffix1_v = NULL, data2_dt, suffix2_v = NULL, geneCol_v = "Gene",
                               pValCol_v = "p_val_adj", lfcCol_v = "avg_log2FC",
                               pThresh_v = 0.05, lfcThresh_v = 0.5, title_v = NULL) {
  #' Plot DEG Concordance
  #' @description
    #' Scatterplot of two different DEG tests, marking genes going in same/diff directions in each
  #' @param data1_dt one DEG table
  #' @param data2_dt another DEG table
  #' @param geneCol_v column name in data1/data2 that indicates the genes
  #' @param pValCol_v column name in data1/data2 that indicates the adjusted p values
  #' @param lfcCol_v column name in data1/data2 that indicates the log fold change
  #' @param pThresh_v value to use for significance thresholding
  #' @param lfcThresh_v value to use for significance thresholding
  #' @return ggplot
  #' @export
  
  # First check column names
  data_lsdt <- list("one" = data1_dt, "two" = data2_dt)
  data_lsdt <- lapply(data_lsdt, function(x) {
    ### Gene columnb
    if (!(geneCol_v %in% colnames(x))) {
      if ("Gene.Name" %in% colnames(x)) {
        colnames(x)[colnames(x) == "Gene.Name"] <- geneCol_v
      } else {
        stop("Check gene name columns!")
      }
    }
    ### P value column
    if (!(pValCol_v %in% colnames(x))) {
      if ("padj" %in% colnames(x)) {
        colnames(x)[colnames(x) == "padj"] <- pValCol_v
      } else {
        stop("Check p value columns!")
      }
    }
    ### LFC column
    if (!(lfcCol_v %in% colnames(x))) {
      if ("log2FoldChange" %in% colnames(x)) {
        colnames(x)[colnames(x) == "log2FoldChange"] <- lfcCol_v
      } else {
        stop("Check lfc columns!")
      }
    }
    return(x)
  })
  
  # Handle suffixes
  if (is.null(suffix1_v)) suffix1_v <- "_1"
  if (is.null(suffix2_v)) suffix2_v <- "_2"
  
  # Now merge them
  merge_dt <- merge(data_lsdt$one[,mget(c(geneCol_v, pValCol_v, lfcCol_v))],
                    data_lsdt$two[,mget(c(geneCol_v, pValCol_v, lfcCol_v))],
                    by = geneCol_v, sort = F, suffixes = c(suffix1_v, suffix2_v))
  
  # New column names
  lfc1_v <- paste0(lfcCol_v, suffix1_v); lfc2_v <- paste0(lfcCol_v, suffix2_v)
  pval1_v <- paste0(pValCol_v, suffix1_v); pval2_v <- paste0(pValCol_v, suffix2_v)
  
  # Add concordance variable
  merge_dt[(get(lfc1_v) > 0 & get(lfc2_v) > 0) |      # both up
             (get(lfc1_v) < 0 & get(lfc2_v) < 0),     # both down
           Concordance := "Concordant"]
  
  merge_dt[(get(lfc1_v) > 0 & get(lfc2_v) < 0) |      # 1 up, 2 down
             (get(lfc1_v) < 0 & get(lfc2_v) > 0),     # 1 down, 2 up
           Concordance := "Discordant"]
  
  # Get data for boxes in the plot
  box_dt <- data.table(xmin = c(-2, 0, -2, 0),
                       xmax = c(0, 2, 0, 2),
                       ymin = c(-2, -2, 0, 0),
                       ymax = c(0, 0, 2, 2),
                       color = c("c", "d", "d", "c"))
  
  # Make plot
  scatter_gg <- ggplot(merge_dt, aes(x = !!sym(lfc1_v), y = !!sym(lfc2_v), 
                                     shape = Concordance, color = Concordance, label = !!sym(geneCol_v))) +
    
    ### Points and labels
    geom_point(size = 3) + 
    scale_color_manual(values = c("blue", "red"), breaks = c("Concordant", "Discordant")) +
    ggrepel::geom_text_repel(max.overlaps = Inf, force = 2, show.legend = F) +
    
    ### Quadrant colors
    geom_rect(data = box_dt, inherit.aes = F, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = color), alpha = 0.1) +
    scale_fill_manual(values = c("blue", "red"), breaks = c("c", "d")) +
    
    ### Quadrant labels
    ggplot2::annotate(geom = "text", label = "Both Down", x = -Inf, y = -Inf, hjust = -0.1, vjust = -1, color = "darkblue") +
    ggplot2::annotate(geom = "text", label = "Both Up", x = Inf, y = Inf, hjust = 1, vjust = 1, color = "darkblue") +
    ggplot2::annotate(geom = "text", label = paste0("Up in ", gsub("_", "", suffix1_v), "\nDown in ", gsub("_", "", suffix1_v)), 
                      x = Inf, y = -Inf, hjust = 1, vjust = -0.2, color = "darkred") +
    ggplot2::annotate(geom = "text", label = paste0("Down in ", gsub("_", "", suffix1_v), "\nUp in ", gsub("_", "", suffix1_v)), 
                      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, color = "darkred") +
    
    
    ### Titles
    labs(title = "DEG Concordance", 
         y = paste0("Avg logFC ", gsub("_", "", suffix2_v)),
         x = paste0("Avg logFC ", gsub("_", "", suffix1_v))) +
    
    ### Misc
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    scale_x_continuous(limits = c(-2, 2)) + scale_y_continuous(limits = c(-2, 2)) +
    my_theme() + guides(fill = "none", shape = guide_legend(override.aes = list(size = 5))) +
    theme(axis.ticks.length = unit(0.3, "cm"))
  
  # Add title
  if (!is.null(title_v)) scatter_gg <- scatter_gg + ggtitle(title_v)
  
  # Return
  return(scatter_gg)
     
}