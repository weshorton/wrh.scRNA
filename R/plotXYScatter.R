plotXYScatter <- function(data_dt, x_v, y_v, title_v = NULL, cor_v = T) {
  #' Plot XY Scatter
  #' @description
  #' Plot XY scatter to compare two variables. 
  #' Include correlation, if desired.
  #' @param data_dt data.table to plot. Must have columns corresponding to x_v and y_v that are numeric
  #' @param x_v column name corresponding to x-axis variable
  #' @param y_v column name corresponding to y-axis variable
  #' @param title_v optional plot title
  #' @param cor_v logical. Indicates whether or not to include regression line and correlation
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
  
  ### Add correlation
  if (cor_v) {
    plot_gg <- plot_gg +
      ggpubr::stat_cor() +
      geom_smooth(method = "lm")
  } # fi cor_v
  
  ### Return
  return(plot_gg)
  
} # plotXYScatter
