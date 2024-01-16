plotXYScatter <- function(data_dt, x_v, y_v, title_v = NULL, cor_v = T, scale_v = NULL, facet_v = NULL, color_v = NULL) {
  #' Plot XY Scatter
  #' @description
  #' Plot XY scatter to compare two variables. 
  #' Include correlation, if desired.
  #' @param data_dt data.table to plot. Must have columns corresponding to x_v and y_v that are numeric
  #' @param x_v column name corresponding to x-axis variable
  #' @param y_v column name corresponding to y-axis variable
  #' @param title_v optional plot title
  #' @param cor_v logical. Indicates whether or not to include regression line and correlation
  #' @param scale_v either NULL (no scale), "both", "x", or "y" indicating where to scale.
  #' @param facet_v either NULL (don't facet), or a column from data_dt to use as input to facet_wrap
  #' @param color_v either NULL (don't split by color) or column form data_dt to color by
  #' @details
  #' Create a standard ggplot xy-scatterplot using the provided parameters.
  #' No option to facet right now, best to pair with ggarrange()
  #' Written to compare lfc and pvals of two different FindMarkers() runs, but can be used for any x-y pair.
  #' @return ggplot object
  #' @export 
  
  ### Check classes
  if (class(data_dt[[x_v]]) != "numeric") stop(sprintf("Plotting column %s isn't numeric.\n", x_v))
  if (class(data_dt[[y_v]]) != "numeric") stop(sprintf("Plotting column %s isn't numeric.\n", y_v))
  
  ### Build plot
  plot_gg <- ggplot2::ggplot(data = data_dt, aes(x = !!sym(x_v), y = !!sym(y_v))) +
    geom_point() +
    my_theme()
  
  ### Add title
  if (!is.null(title_v)) {
    plot_gg <- plot_gg + ggtitle(title_v)
  } # fi is.null
  
  ### Scale
  if (!is.null(scale_v)) {
    
    if (scale_v == "x" | scale_v == "both") {
      plot_gg <- plot_gg + scale_x_log10()
    } 
    
    if (scale_v == "y" | scale_v == "both") {
      plot_gg <- plot_gg + scale_y_log10()
    } # fi
    
  } # fi
  
  ### Facet
  if (!is.null(facet_v)) {
    plot_gg <- plot_gg + facet_wrap(as.formula(paste0("~", facet_v)))
  }
  
  ### Color
  if (!is.null(color_v)) {
    plot_gg <- plot_gg + aes(color = !!sym(color_v))
  }
  
  ### Add correlation
  if (cor_v) {
    plot_gg <- plot_gg +
      ggpubr::stat_cor() +
      geom_smooth(method = "lm")
  } # fi cor_v
  
  ### Return
  return(plot_gg)
  
} # plotXYScatter
