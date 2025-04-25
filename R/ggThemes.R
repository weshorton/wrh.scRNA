umapFigureTheme <- function() {
  #' UMAP Theme For Figures
  #' @description Larger text all-around, thicker axis.
  #' @export
  
  theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 36),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 28),
          axis.ticks = element_line(linewidth = 0.25),
          axis.ticks.length = unit(0.25, "cm"),
          axis.line = element_line(linewidth = 0.25),
          strip.text = element_text(size = 24),
          legend.text = element_text(size = 24),
          legend.title = element_text(size = 28))
  
} # umapFigureTheme

waterfallFigureTheme <- function() {
  #' Waterfall Theme For Figures
  #' @description Larger text all-around (even larger than umap), thicker axis.
  #' @export
  
  theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 42),
          axis.text = element_text(size = 36),
          axis.title = element_text(size = 40),
          axis.ticks = element_line(linewidth = 2),
          axis.ticks.length = unit(0.25, "cm"),
          axis.line = element_line(linewidth = 2),
          strip.text = element_text(size = 36),
          legend.text = element_text(size = 36),
          legend.title = element_text(size = 40))
  
} # umapFigureTheme

dotPlotTheme <- function() {
  #' scRNA Dot Plot Theme
  #' @description Customized theme for dotplots. Specific grid and border lines
  #' @export
  
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray",linetype = "dashed", linewidth=0.35), 
        panel.border = element_rect(colour = "black", fill=NA, linewidth=2)) + 
    theme(legend.position="bottom", plot.title = element_text(hjust = 0.5))
  
} # dotPlotTheme


figureDotTheme <- function() {
  #' scRNA Figure Dot Theme
  #' @description
    #' Customized theme for dotplots
  #' @export
  
  #wrh.rUtils::big_label() +
    wrh.rUtils::angle_x() +
    theme(line = element_line(linewidth = 0)) +
    theme(rect = element_rect(linewidth = 0.0)) +
    theme(panel.border = element_rect(linewidth = 0.0)) +
    theme(axis.line = element_line(linewidth = 0.0)) +
    theme(legend.position="bottom", legend.box = "vertical")
    
} # figureDotTheme

transparentTheme <- function() {
  #' Transparent Everything
  #' @description try to make as few artifacts as possible
  #' @export
  
  my_theme() + theme(
    panel.background = element_rect(fill = "transparent", 
                                    colour = NA_character_), # necessary to avoid drawing panel outline
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    plot.background = element_rect(fill = "transparent",
                                   colour = NA_character_), # necessary to avoid drawing plot outline
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent"))
}