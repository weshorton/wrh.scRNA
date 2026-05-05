#' Plot expression of a receptor's ligands by other cell types as a chord plot
#'
#' Creates a chord plot of expression of ligands that can activate a specified
#' receptor where chord widths correspond to mean ligand expression by the cluster.
#'
#' @param signaling_dt Data.table of expression
#' @param dest name of destination gene
#' @param expression_threshold Minimum mean expression value of a ligand by a cell type for a chord to be rendered between the cell type and the receptor
#' @param cellCol_v column name for cell types - usually either 'Pop' or 'Treatment'
#' @param cell_idents Vector of cell types from cluster assignments in the domino object to be included in the plot.
#' @param cell_colors Named vector of color names or hex codes where names correspond to the plotted cell types and the color values
#' @param title if provided, final title will be paste0(receptor/ligand, " Signaling")
#' @param label_arcs logical indicating if arc populations should be labeled
#' @param multi_plot numeric vector that, if present indicates which plot in a list of plots this one is (for legend output)
#' @param destExpr_v expression controling arc width for destination gene
#' @param max_v optional global maximum
#' @return renders a circos plot to the active graphics device
#' @export circos_ligand_receptor_general
#' 
newCircos <- function(signaling_dt,
                      dest,
                      arc_genes,
                      expression_threshold = 0.01, 
                      cellCol_v = "Pop",
                      cell_idents = NULL,
                      cell_colors = NULL, 
                      title = "Signaling", 
                      label_arcs = T, 
                      multi_plot = NULL,
                      destExpr_v, # testing
                      max_v = NULL) {
  
  ### Helper fxn
  ggplot_col_gen <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    return(grDevices::hcl(h = hues, l = 65, c = 100)[1:n])
  }
  
  ### Convert
  signaling_df <- as.data.frame(signaling_dt)
  
  ### Subset
  signaling_df <- signaling_df[signaling_df[[cellCol_v]] %in% cell_idents,]
  signaling_df <- signaling_df[grep(paste(arc_genes, collapse = "|"), signaling_df$origin),]
  
  ### Scale - depends on if global maximum is provided
  if (is.null(max_v)) {
    signaling_df$scaled.mean.expression <- signaling_df$mean.expression / max(signaling_df$mean.expression)
    scaledDestExp_v <- destExpr_v / max(signaling_df$mean.expression)
  } else {
    signaling_df$scaled.mean.expression <- signaling_df$mean.expression / max_v
    scaledDestExp_v <- destExpr_v / max_v
  }
  
  ### Handle cell idents
  if (is.null(cell_idents)) {
    cell_idents <- unique(signaling_dt[[cellCol_v]])
  }
  
  ### Handle cell colors
  if (is.null(cell_colors)) {
    cell_colors <- ggplot_col_gen(n = length(cell_idents))
    names(cell_colors) <- cell_idents
  } else {
    cell_colors <- cell_colors[cell_idents]
    if (length(cell_colors) != length(cell_idents)) stop("mismtach")
  }
  
  ### Ligand arcs
  arc_df <- signaling_df[, c("origin", "destination")]
  arc_df["ligand.arc"] <- 1
  
  # receptor arc will always sum to 4 no matter how many ligands and cell idents are plotted
  # arc_df["receptor.arc"] <- 4 / (nrow(signaling_df))
  arc_df["receptor.arc"] <- 2 / (nrow(signaling_df)) # try 2 to make this portion narrower
  
  # name grouping based on [cell_ident]
  nm <- c(dest, arc_df$origin)
  sub_nm <- sapply(nm[-1], function(x) {
    y <- strsplit(x, split = "-")[[1]]
    y <- y[1:length(y)-1]
    y <- paste0(y, collapse = "-")}, USE.NAMES = F)
  group <- structure(c(nm[1], sub_nm), names = nm)
  
  # order group as a factor with the receptor coming first
  group <- factor(group, levels = c(
    #dest, sort(unique(sub_nm)) # alphabetical order of the other cell idents
    dest, cell_idents # use provided order
  ))
  
  # colors for ligand chords
  arc_colors <- ggplot_col_gen(length(arc_genes))
  names(arc_colors) <- arc_genes
  
  grid_col <- c("#FFFFFF") # hide the arc corresponding to the receptor by coloring white
  for (i in 1:length(arc_colors)) {
    grid_col <- c(grid_col, rep(arc_colors[i], length(cell_idents)))
  }
  names(grid_col) <- c(dest, signaling_df$origin)
  
  ### Determine which dest expr to use
  yesScale_v <- ((max(signaling_df$mean.expression) > 1) | destExpr_v > 1)
  if (yesScale_v) {
    destExpr_v <- scaledDestExp_v
  } # fi
  
  # Start plot
  circlize::circos.clear()
  circlize::circos.par(start.degree = 0, circle.margin = 0.5)
  circlize::chordDiagram(arc_df,
                         group = group, grid.col = grid_col, link.visible = FALSE, annotationTrack = c("grid"),
                         preAllocateTracks = list(track.height = circlize::mm_h(4), track.margin = c(circlize::mm_h(2), 0)), big.gap = 2
  )
  
  ### Make title
  if (!is.null(title)) title(paste0(dest, " ", title))
  
  ### Draw each link
  for (send in signaling_df$origin) {
    
    ### Skip if not above threshold
    if (signaling_df[signaling_df$origin == send, "mean.expression"] < expression_threshold) next
    
    ### Determine if scaled or unscaled should be used
    if (yesScale_v) {
      expr <- signaling_df[signaling_df$origin == send, "scaled.mean.expression"]
      max_width <- signif(max(signaling_df$mean.expression), 2)
    } else {
      expr <- signaling_df[signaling_df$origin == send, "mean.expression"]
      max_width <- 1
    } # fi
    
    ### Handle if global max is provided
    max_width <- ifelse(is.null(max_v), max_width, signif(max_v, 2))
    
    ### Draw link
    circlize::circos.link(sector.index1 = send,
                          point1 = c(0.5 - (expr / 2), 0.5 + (expr / 2)), # do this to make sure arc is centered.
                          sector.index2 = dest,
                          point2 = c(1 - (destExpr_v/2), 1 + (destExpr_v/2)), # just have this be 2 if we don't want Mif width
                          col = paste0(grid_col[[send]], "88"))
    
  } # for send
  
  
  sector_names <- circlize::get.all.sector.index()
  sub_sector_names <- c(sector_names[1], sapply(sector_names[-1], function(x) {
    y <- strsplit(x, split = "-")[[1]]
    y <- y[1:length(y)-1]
    y <- paste0(y, collapse = "-")}, USE.NAMES = F))
  cell_sectors <- cell_idents[cell_idents %in% sub_sector_names]
  for (cell in cell_sectors) {
    row_pick <- sector_names[grepl(paste0("^", cell), sector_names)]
    if (length(row_pick)) {
      if (label_arcs) { cell_text <- cell } else { cell_text <- NULL }
      circlize::highlight.sector(sector_names[grepl(paste0("^", cell, "-"), sector_names)],
                                 track.index = 1,
                                 col = cell_colors[cell], text = cell_text, cex = 1, facing = "inside", text.col = "black",
                                 niceFacing = FALSE, text.vjust = -1.5
      )
    }
  }
  
  circlize::highlight.sector(sector_names[grepl(paste0("^", dest, "$"), sector_names)],
                             track.index = 1,
                             col = "#FFFFFF", text = dest, cex = 1.5, facing = "clockwise", text.col = "black", niceFacing = TRUE,
                             pos = 4
  )
  
  # create legends
  lgd_cells <- ComplexHeatmap::Legend(
    at = as.character(cell_idents), type = "grid", legend_gp = grid::gpar(fill = cell_colors),
    title_position = "topleft", title = "cell identity"
  )
  lgd_ligands <- ComplexHeatmap::Legend(
    at = arc_genes, type = "grid", legend_gp = grid::gpar(fill = arc_colors), title_position = "topleft",
    title = "ligand"
  )
  #chord_width <- 10 / (4 + length(cell_idents) * length(arc_genes)) # not entirely sure where 10 and 4 come from
  chord_width <- 10 / (2 + length(cell_idents) * length(arc_genes)) # 2 instead?
  lgd_chord <- ComplexHeatmap::Legend(
    at = c(expression_threshold, max_width), col_fun = circlize::colorRamp2(c(
      expression_threshold,
      max_width
    ), c("#DDDDDD", "#DDDDDD")), legend_height = grid::unit(chord_width, "in"), title_position = "topleft",
    title = "expression"
  )
  if (is.null(multi_plot)) {
    lgd_list_vertical <- ComplexHeatmap::packLegend(lgd_cells, lgd_ligands, lgd_chord)
    ComplexHeatmap::draw(lgd_list_vertical, x = grid::unit(0.02, "npc"), y = grid::unit(0.98, "npc"), just = c("left", "top"))
  } else {
    lgd_list_vertical <- ComplexHeatmap::packLegend(lgd_cells, lgd_ligands, lgd_chord)
    #lgd_list_vertical <- lgd_chord
    if (multi_plot %% 2 == 0) { 
      ComplexHeatmap::draw(lgd_list_vertical, x = grid::unit(0.9, "npc"), y = grid::unit(0.98-(0.13*(multi_plot-2)), "npc"), just = c("left", "top"))
    } else { 
      ComplexHeatmap::draw(lgd_list_vertical, x = grid::unit(0.02, "npc"), y = grid::unit(0.98-(0.13*(multi_plot-1)), "npc"), just = c("left", "top"))
    }
    
  }
  
} # newCircos