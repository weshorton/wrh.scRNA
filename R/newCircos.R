#' Plot expression of a receptor's ligands by other cell types as a chord plot
#'
#' Creates a chord plot of expression of ligands that can activate a specified
#' receptor where chord widths correspond to mean ligand expression by the cluster.
#'
#' @param signaling_dt Data.table of signaling expression. Must contain columns
#'   for cell type (cellCol_v), 'origin', 'destination', and 'mean.expression'.
#' @param dest name of destination gene
#' @param arc_genes name(s) of ligands to include as arcs
#' @param expression_threshold Minimum mean expression value of a ligand by a cell type for a chord to be rendered between the cell type and the receptor
#' @param cellCol_v column name for cell types - usually either 'Pop' or 'Treatment'
#' @param cell_idents Vector of cell types from cluster assignments in the domino object to be included in the plot.
#' @param cell_colors Named vector of color names or hex codes where names correspond to the plotted cell types and the color values
#' @param title if provided, final title will be paste0(receptor/ligand, " Signaling")
#' @param label_arcs logical indicating if arc populations should be labeled
#' @param multi_plot numeric vector that, If provided, indicates position of this
#'   plot in a multi-panel layout for legend placement. Even = right side,
#'   odd = left side.
#' @param destExpr_v expression controlling arc width for destination gene
#' @param max_v optional global maximum. If NULL, uses local maximum
#' @param plotHeight_v height in inches of plot. Used to make the scale legend appropriately sized
#' @details
#' Sender arc widths are always set to 1.0 (unit width) and colored by ligand.
#' The destination arc width is set to 2 / nrow(signaling_df) so it remains
#' narrow relative to the sender arcs. Expression values are scaled to [0, 1]
#' relative to max_v (or local max) before being used as link widths.
#'
#' The lgd_chord legend bar height is proportional to the sender arc width
#' relative to the total arc space: 1 / (nrow(signaling_df) + 1).
#' @return renders a circos plot to the active graphics device
#' @export newCircos
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
                      destExpr_v,
                      max_v = NULL,
                      plotHeight_v = 7) {
  
  ### Helper fxn
  ggplotColGen <- function(n) {
    hues_v <- seq(15, 375, length = n + 1)
    return(grDevices::hcl(h = hues_v, l = 65, c = 100)[1:n])
  } # ggplotColGen
  
  ### ============================================================
  ### Data preparation
  ### ============================================================
  
  ### Convert
  signaling_df <- as.data.frame(signaling_dt)
  
  ### Subset
  signaling_df <- signaling_df[signaling_df[[cellCol_v]] %in% cell_idents,]
  signaling_df <- signaling_df[grep(paste(arc_genes, collapse = "|"), signaling_df$origin),]
  
  ### Handle cell idents
  if (is.null(cell_idents)) {
    cell_idents <- unique(signaling_dt[[cellCol_v]])
  }
  
  ### ============================================================
  ### Expression scaling
  ### ============================================================
  
  ### Scale expression values to [0, 1] relative to max
  ### Use global max_v if provided, otherwise use local maximum
  
  ### Get appropriate max
  localMax_v <- max(signaling_df$mean.expression)
  scaleMax_v <- ifelse(is.null(max_v), localMax_v, max_v)
  
  ### Apply scale
  signaling_df$scaled.mean.expression <- signaling_df$mean.expression / scaleMax_v
  scaledDestExpr_v <- destExpr_v / scaleMax_v
  
  ### Use scaled values if any expression exceeds 1 (i.e., not already in [0,1])
  yesScale_v <- (localMax_v > 1) | (destExpr_v > 1)
  if (yesScale_v) {
    destExpr_v <- scaledDestExpr_v
  } # fi
  
  ### ============================================================
  ### Colors
  ### ============================================================
  
  ### Handle cell colors
  if (is.null(cell_colors)) {
    cell_colors <- ggplot_col_gen(n = length(cell_idents))
    names(cell_colors) <- cell_idents
  } else {
    cell_colors <- cell_colors[cell_idents]
    if (length(cell_colors) != length(cell_idents)) stop("mismtach between names of cell_colors and cell_idents")
  } # fi
  
  ### Arc (ligand) colors — one color per arc gene
  arc_colors_v <- ggplotColGen(length(arc_genes))
  names(arc_colors_v) <- arc_genes
  
  ### Grid colors: destination sector is white (hidden), senders colored by ligand
  grid_col_v <- c("#FFFFFF", rep(arc_colors_v, each = length(cell_idents)))
  names(grid_col_v) <- c(dest, signaling_df$origin)
  
  ### ============================================================
  ### Sector layout
  ### ============================================================
  
  ### Build arc data.frame: origin -> destination with arc widths
  arc_df <- signaling_df[, c("origin", "destination")]
  arc_df$ligand.arc   <- 1                              ### Sender arc: fixed unit width
  arc_df$receptor.arc <- 2 / nrow(signaling_df)        ### Dest arc: narrow, fixed total (was originally 4, changed to 2 to be narrower)
  
  ### Build grouping: maps each sector name to its parent cell type
  ### The destination sector maps to itself; senders map to their cell type
  ### (Do the sapply/strsplit instead of gsub in case there are other "-" in the names)
  sectorNames_v <- c(dest, arc_df$origin)
  senderParents_v <- sapply(sectorNames_v[-1], function(x) {
    parts_v <- strsplit(x, split = "-")[[1]]
    paste0(parts_v[1:(length(parts_v) - 1)], collapse = "-")
  }, USE.NAMES = FALSE)
  
  group_v <- c(sectorNames_v[1], senderParents_v)
  names(group_v) <- sectorNames_v
  
  ### Order group factor: destination first, then cell types in provided order
  group_v <- factor(group_v, levels = c(dest, cell_idents))
  
  ### Determine max width for legend (if we don't scale, the max is just 1)
  max_width_v <- ifelse(yesScale_v, 
                        signif(max(signaling_df$mean.expression), 2),
                        1)
  
  ### have to handle if global max is provided
  max_width_v <- ifelse(is.null(max_v),
                        max_width_v,
                        signif(max_v, 2))
  
  ### ============================================================
  ### Initialize circos plot
  ### ============================================================
  
  circlize::circos.clear()
  circlize::circos.par(start.degree = 0, circle.margin = 0.5)
  
  circlize::chordDiagram(arc_df,
                         group = group_v, 
                         grid.col = grid_col_v, 
                         link.visible = FALSE,  # drawn manually below
                         annotationTrack = c("grid"),
                         preAllocateTracks = list(track.height = circlize::mm_h(4), 
                                                  track.margin = c(circlize::mm_h(2), 0)), 
                         big.gap = 2)
  
  ### Make title
  if (!is.null(title)) graphics::title(paste0(dest, " ", title))
  
  ### ============================================================
  ### Draw links
  ### ============================================================
  
  for (send in signaling_df$origin) {
    
    ### Get expression
    sendExpr_v <- signaling_df[signaling_df$origin == send, "mean.expression"]
    
    ### Skip if expression is below threshold
    if (sendExpr_v < expression_threshold) next
    
    ### Use scaled or raw expression for link width
    if (yesScale_v) {
      linkWidth_v <- signaling_df[signaling_df$origin == send, "scaled.mean.expression"]
    } else {
      linkWidth_v <- sendExpr_v
    } # fi
    
    ### Draw link - sender arcs are centered at 0.5, dest arc is centered at 1.0
    ### For future - point2 can be a single value if we don't want an arc for the dest
    circlize::circos.link(sector.index1 = send,
                          point1 = c(0.5 - (linkWidth_v / 2), 0.5 + (linkWidth_v / 2)), # do this to make sure arc is centered.
                          sector.index2 = dest,
                          point2 = c(1 - (destExpr_v/2), 1 + (destExpr_v/2)), 
                          directional = -1,
                          arr.type = "big.arrow",
                          col = paste0(grid_col_v[[send]], "88")) # 88 gives 53% opacity
    
  } # for send
  
  ### ============================================================
  ### Sector highlighting and labels
  ### ============================================================
  
  sector_names <- circlize::get.all.sector.index()
  
  ### Identify which cell types are actually present as sectors
  ### Same substitution for senderParents above
  sectorParents_v <- c(sectorNames_v[1],                        # destination
                       sapply(sectorNames_v[-1], function(x) {
                         parts_v <- strsplit(x, split = "-")[[1]]
                         paste0(parts_v[1:(length(parts_v) - 1)], collapse = "-")
                         }, USE.NAMES = FALSE))
  
  presentCells_v <- cell_idents[cell_idents %in% unique(sectorParents_v)]
  
  ### Highlight each cell type's sectors with its color
  for (cell_v in presentCells_v) {
    
    ### Build grep - must be at beginning and followed by -
    grep_v <- paste0("^", cell_v, "-")
    
    ### Get sectors that correspond with this cell
    cellSectors_v <- grep(grep_v, sectorNames_v, value = T) 
    
    ### Skip if none
    if (length(cellSectors_v) == 0) next
    
    ### Make label if option is selected
    cellLabel_v <- if (label_arcs) cell_v else NULL
    
    ### Do the highlight
    circlize::highlight.sector(sector.index = cellSectors_v,
                               track.index = 1,
                               col = cell_colors[cell_v],
                               text = cellLabel_v, 
                               cex = 1, 
                               facing = "inside", 
                               text.col = "black", 
                               niceFacing = F, 
                               text.vjust = -1.5)
  } # for cell_v
  
  ### Highlight destination sector (white background, labeled)
  ### CUSTOM - ONLY FOR MIF PROJECT AS OF NOW
  circlize::highlight.sector(sector.index = sector_names[grepl(paste0("^", dest, "$"), sector_names)],
                             track.index = 1,
                             col = "#FFFFFF",
                             text = "Myeloid",
                             cex = 1.5,
                             facing = "outside",
                             text.col = "black",
                             niceFacing = F,
                             text.vjust = 3
                             )
  
  ### Highlight destination sector (white background, labeled)
  circlize::highlight.sector(sector.index = sector_names[grepl(paste0("^", dest, "$"), sector_names)],
                             track.index = 1,
                             col = "#FFFFFF",
                             text = dest,
                             cex = 1.5,
                             facing = "clockwise",
                             text.col = "black",
                             niceFacing = T,
                             #text.vjust = -1.5
                             )
  
  ### ============================================================
  ### Legends
  ### ============================================================
  
  ### Cell identity legend
  lgd_cells <- ComplexHeatmap::Legend(at = as.character(cell_idents),
                                      type = "grid",
                                      legend_gp = grid::gpar(fill = cell_colors),
                                      title_position = "topleft",
                                      title = "cell identity")
  
  ### Ligand legend
  lgd_ligands <- ComplexHeatmap::Legend(at = arc_genes,
                                        type = "grid",
                                        legend_gp = grid::gpar(fill = arc_colors_v),
                                        title_position = "topleft",
                                        title = "ligand")
  
  ###
  ### Expression Scale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### This one is a bit complicated.
  ### Height represents the max sender arc width as a fraction of the total arc space
  ### Total arc space is the sender arcs and the destination arc
  ### total sender arcs is simply nrow(signaling_df)
  ### Destination arc is hardcoded to be 2
  ### Have to then convert that to inches based on the plot height
  arcSpace_v <- nrow(signaling_df) + 2
  arcFraction_v <- 1 / arcSpace_v
  #chordWidth_v <- plotHeight_v / arcSpace_v
  chordWidth_v <- 10 / (4 + nrow(signaling_df))
  
  ### Build legend using chordWidth as height of legend
  lgd_chord <- ComplexHeatmap::Legend(at = c(expression_threshold, max_width_v),
                                      col_fun = circlize::colorRamp2(breaks = c(expression_threshold, max_width_v), 
                                                                     colors = c("#DDDDDD", "#DDDDDD")),
                                      legend_height = grid::unit(chordWidth_v, "in"),
                                      title_position = "topleft",
                                      title = "expression")
  
  # ### Above doesn't seem to be working. Try using grid instead
  # ### Total arc space in canvas units (full canvas = 2.0)
  # arcSpace_v <- nrow(signaling_df) + 2
  # #totalArcSpace_n <- nrow(signaling_df) + 2
  # 
  # ### Height of one sender arc in canvas coordinates
  # ### Canvas spans 2.0 units total (-1 to 1)
  # barHeight_v <- (1 / arcSpace_v) * 2
  # #barHeight_n <- (1.0 / totalArcSpace_n) * 2.0
  # 
  # ### Position the bar in the lower-left of the canvas (avoid the other legends)
  # xLeft_v   <- -1.3
  # yBottom_v <- -1.2
  # 
  # graphics::rect(
  #   xleft   = xLeft_v,
  #   ybottom = yBottom_v,
  #   xright  = xLeft_v + 0.15,       ### Bar width — adjust for aesthetics
  #   ytop    = yBottom_v + barHeight_v,
  #   col     = "#DDDDDD",
  #   border  = "black",
  #   xpd     = NA                     ### Allow drawing outside plot region
  # )
  # 
  # ### Label above the bar
  # graphics::text(
  #   x      = xLeft_v + 0.075,
  #   y      = yBottom_v + barHeight_v + 0.05,
  #   labels = paste0("expr = ", signif(max_width_v, 2)),
  #   cex    = 0.7,
  #   adj    = c(0.5, 0),
  #   xpd    = NA
  # )
  # 
  # ### Section label below the bar
  # graphics::text(
  #   x      = xLeft_v + 0.075,
  #   y      = yBottom_v - 0.05,
  #   labels = "arc width",
  #   cex    = 0.6,
  #   adj    = c(0.5, 1),
  #   col    = "grey30",
  #   xpd    = NA
  # )
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ### Pack legends
  lgd_list_vertical <- ComplexHeatmap::packLegend(lgd_cells, lgd_ligands, lgd_chord)
  
  ### ============================================================
  ### Multi-Plot Legends
  ### ============================================================
  
  if (is.null(multi_plot)) {
    
    ### Single plot: top-left corner
    ComplexHeatmap::draw(lgd_list_vertical, 
                         x = grid::unit(0.02, "npc"), 
                         y = grid::unit(0.98, "npc"), 
                         just = c("left", "top"))
  } else {
    
    ### Multi-plot: alternate left/right based on plot index parity
    ### Even indices go on the right, odd indices go on the left
    if (multi_plot %% 2 == 0) {
      xPos_v <- grid::unit(0.9, "npc")
      yPos_v <- grid::unit(0.98-(0.13*(multi_plot-2)), "npc")
    } else { 
      xPos_v <- grid::unit(0.02, "npc")
      yPos_v <- grid::unit(0.98-(0.13*(multi_plot-1)), "npc")
    } # fi
    
    ### Draw the legend using above locations
    ComplexHeatmap::draw(lgd_list_vertical, 
                         x = xPos_v, 
                         y = yPos_v,
                         just = c("left", "top"))
    
  } # fi
  
  ### ============================================================
  ### Output
  ### ============================================================
  
  ### Return processed data invisibly for downstream use
  # invisible(list(
  #   signaling_df  = signaling_df,
  #   arc_colors_v  = arc_colors_v,
  #   cell_colors_v = cell_colors
  # ))
  
} # newCircos