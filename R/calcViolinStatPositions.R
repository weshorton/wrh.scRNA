calcViolinStatPositions <- function(plot_df, measureVar_v, cm_df, levels_v, yOffset_v = 1.56) {
  #' Calculate Violin Stat Positions
  #' @description determine xy coordinates for stat::compare_means results
  #' @param plot_df data.frame used for plotting
  #' @param measureVar_v name of column in plot_df that contains metric to be plotted
  #' @param cm_df data.frame of compare_means output
  #' @param levels_v factor levels. needed to get x-position
  #' @param yOffset_v amount to multiply maximum y-value by to set y-limit
  #' @return named list of length 2:
  #' 1. "cm_df" contains cm_df with the added xPos and yPos columns
  #' 2. "ylim" contains y limits for the plot
  #' @export
  
  # X position
  g1x_v <- as.numeric(factor(cm_df$group1, levels = levels_v))
  g2x_v <- as.numeric(factor(cm_df$group2, levels = levels_v))
  cm_df$xPos <- g1x_v + (g2x_v-g1x_v)/2
  
  # Y position
  yPosBase_v <- round(max(plot_df[[measureVar_v]]), digits = 2)
  yPos_v <- sapply(seq_along(1:nrow(cm_df)), function(x) yPosBase_v + (yPosBase_v*.1)*x)
  cm_df$yPos <- yPos_v

  # Y limits
  ylim_v <- c(min(plot_df[[measureVar_v]]), max(plot_df[[measureVar_v]])*yOffset_v)
  
  # Output
  return(list("cm_df" = cm_df, "yPosBase" = yPosBase_v, "ylim" = ylim_v))
  
} # calcViolinStatPositions

