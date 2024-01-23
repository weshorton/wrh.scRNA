umapFigureTheme <- function() {
  #' UMAP Theme For Figures
  #' @description Larger text all-around, thicker axis.
  #' @export
  
  theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 36),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 28),
          axis.ticks = element_line(linewidth = 2),
          axis.ticks.length = unit(0.25, "cm"),
          axis.line = element_line(linewidth = 2),
          strip.text = element_text(size = 24),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 28))
  
} # umapFigureTheme