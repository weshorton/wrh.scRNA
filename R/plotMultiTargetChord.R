#' Plot Asymmetric Multi-Target Chord Diagram
#'
#' @description
#' Generates a chord diagram for a single sender population linking to all
#' receiver populations. 
#'
#' @param circData_dt Cleaned data.table.
#' @param senderPop_v Population to act as the sender
#' @param receiverPop_v The populations to act as receiver. NULL (default) uses all.
#' @param senderDataCol_v column that has values for sender
#' @param receiverDataCol_v column that has values for receiver
#' @param popCol_v character vector indicating which column of circData_dt has sector info
#' @param sectorColors_v optional vector to color sectors
#' @param arcColorCol_v column indicating what to color connectors by (usually treatment)
#' @param arcColors_v Named vector of colors.
#' @param legendType_v one of 'both', 'sector', 'arc', or 'none' indicating what to include in legend
#' @param cex_v size of labels
#' @param axisTicks_v logical indicating to draw ticks or not
#' @param tickScaler_v what value should the tick mark scale bar show?
#' @param specialLabel_v not sure how to add two sector labels, so have to use this and run a second time for combo figure.
#' @param startDegree_v degree on circle that arcs begin from. Default is 270, but have to change depending on relative sizes of each arc
#' @return print plot
#' @details Run with specialLabel_v = F for standard output. Change to true to get the sender/receiver labels that go 
#' outside the plot.
#' @export
plotMultiTargetChord <- function(circData_dt, 
                                 senderPop_v, 
                                 receiverPop_v = NULL,
                                 senderDataCol_v, receiverDataCol_v,
                                 popCol_v = "Pop",
                                 sectorColors_v = NULL,
                                 arcColorCol_v = "Treatment",
                                 arcColors_v = NULL,
                                 legendType_v = "both",
                                 axisTicks_v = F,
                                 tickScaler_v = 5,
                                 cex_v = 1,
                                 specialLabel_v = F,
                                 startDegree_v = 270) {
  
  ### Wrangle input data
  links_dt <- prepareMultiTargetData(circData_dt, senderPop_v, receiverPop_v, 
                                     senderDataCol_v, receiverDataCol_v, popCol_v, arcColorCol_v)
  
  ###
  ### Sector Sizes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Sectors need x-axis limits that go from 0 to the total width of the sender/reciever
  ### The input table has these divided by the links (treatments), so need to sum the individual parts
  ### For polar coordinate plots, the x-axis is the arc/sector length
  
  ### Define Sectors - have to have suffixes to avoid duplicate error
  senderSector_v <- paste0(senderPop_v, "_Sender")
  receiverSectors_v <- paste0(unique(links_dt$receiverPop), "_Receiver")
  
  ### Sender
  senderMax_v <- sum(links_dt$senderWidthPart) # sum(unique(links_dt[,mget(c(senderDataCol_v, arcColorCol_v))][,senderDataCol_v]))
  
  ### Each receiver sector is the sum of all the arcs going to it
  receiverLimits_dt <- links_dt[, .(xmax = sum(receiverWidthPart)), by = receiverPop]
  receiverLimits_dt[, sectorName := receiverPop]
  
  ### Combine all sectors, which gives total polar coordinate limits
  allSectors_v <- c(senderSector_v, receiverSectors_v)
  allXmax_v <- c(senderMax_v, receiverLimits_dt$xmax)
  
  ###
  ### Colors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Link Colors
  arcs_v <- unique(links_dt[[arcColorCol_v]])
  if (is.null(arcColors_v)) {
    arcColors_v <- stats::setNames(object = grDevices::hcl.colors(length(arcs_v), "Zissou 1"), 
                                   nm = arcs_v)
  } else {
    arcColors_v <- arcColors_v[names(arcColors_v) %in% arcs_v]
    if (length(arcColors_v) != length(arcs_v)) stop("Provided arc colors don't match with arcs. Fix or don't provide.\n")
    arcColors_v <- arcColors_v[match(arcs_v, names(arcColors_v))]
  }
  
  ### Sector Colors
  if (is.null(sectorColors_v)) {
    sectors_v <- unique(links_dt[[popCol_v]])
    sectorColors_v <- stats::setNames(object = grDevices::hcl.colors(length(sectors_v), "viridis"), 
                                      nm = sectors_v)
  } ### End if
  
  ###
  ### Label location ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### If the sector is small, it's possible that the label will be wider than the sector itself
  ### so it will need to be offset. The minimum width depends on the length of the label and cex
  
  ### Get length of labels
  cleanSectors_v <- gsub("_Sender|_Receiver", "", allSectors_v)
  sectorNchar_v <- nchar(cleanSectors_v)
  
  ### Ratio of sector width to label size
  sectorLabelRatio_v <- allXmax_v / (sectorNchar_v * cex_v)
  
  sectorLabelLocation_v <- numeric(length = length(cleanSectors_v))
  for (i in 1:length(sectorLabelRatio_v)) {
    currRatio_v <- sectorLabelRatio_v[i]
    otherRatio_v <- ifelse(i == 1, 
                           sectorLabelRatio_v[length(sectorLabelRatio_v)],
                           sectorLabelRatio_v[(i-1)])
    
    if (currRatio_v > 0.25) {
      sectorLabelLocation_v[i] <- 0.5
    } else {
      if (otherRatio_v <= 0.25) {
        sectorLabelLocation_v[i] <- 1.25
      } else {
        sectorLabelLocation_v[i] <- 1.75
      } # fi otherRatio
    } # fi currRatio
  } # for i
  
  ###
  ### Legend ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  legendArcColors_v <- arcColors_v[intersect(names(arcColors_v), unique(links_dt[[arcColorCol_v]]))]
  legendSectorColors_v <- sectorColors_v[intersect(names(sectorColors_v), unique(gsub("_(Sender|Receiver)$", "", allSectors_v)))]
  
  if (legendType_v == "both") {
    legendColors_v <- c(legendArcColors_v, legendSectorColors_v)
  } else if (legendType_v == "arc") {
    legendColors_v <- legendArcColors_v
  } else if (legendType_v == "sector") {
    legendColors_v <- legendSectorColors_v
  } else if (legendType_v != "none") {
    warning(sprintf("Unknown legendType_v provided - %s. Not drawing legend.\n", legendType_v))
    legendType_v <- "none"
    legendColors_v <- c(legendArcColors_v, legendSectorColors_v)
  }
  
  lgd <- ComplexHeatmap::Legend(labels = names(legendColors_v),
                                type = "points",
                                background = "white",
                                size = unit(5, 'mm'),
                                row_gap = unit(0.5, 'mm'),
                                legend_gp = grid::gpar(col = legendColors_v))
  
  ###
  ### Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Make sure plot area is cleared
  circlize::circos.clear()
  
  ### Have to add this to avoid errors when sector widths are very small
  circos.par(cell.padding = c(0,0))
  
  ### Large gap between sender side (left) and receiver side (right)
  gaps_v <- c(10, rep(2, length(receiverSectors_v) - 1), 10)
  circlize::circos.par(gap.after = gaps_v, 
                       start.degree = startDegree_v,
                       canvas.xlim = c(-1.1, 1.1),
                       canvas.ylim = c(-1.1, 1.1))
  
  ### Initialize sectors (anything added here is searchable with get.cell.meta.data)
  circlize::circos.initialize(sectors = allSectors_v, 
                              xlim = cbind(rep(0, length(allXmax_v)), allXmax_v))
  
  ### Draw Outer Track (Sectors)
  circlize::circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    
    ### Get sector-specific info from initialization
    sect <- circlize::get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim <- circlize::get.cell.meta.data("xlim")
    
    ### Other plotting info
    cleanName_v <- gsub("_(Sender|Receiver)$", "", sect)
    sectorColor_v <- sectorColors_v[cleanName_v]
    sectorY_v <- sectorLabelLocation_v[i]
    # sectorY_v <- 1.25
    
    ### Build sector 'rectangle'
    circlize::circos.rect(xlim[1], 0, xlim[2], 1, col = sectorColor_v, border = "black")
    
    ### Add ticks
    if (axisTicks_v) {
      circlize::circos.axis(h = "top",
                            major.tick.length = 0.2,
                            labels.cex = 0.4,
                            direction = "inside")
    }
    
    # Labels
    circlize::circos.text(
      x = mean(xlim),  # why mean? AI weirdness?
      #x = xlim,
      y = sectorY_v,
      labels = cleanName_v,
      cex = cex_v, 
      facing = "bending.inside", 
      niceFacing = TRUE)
  }, track.height = 0.15, bg.border = NA)
  
  ### Calculate Running Totals for Link Placement
  links_dt[, senderStart := c(0, head(cumsum(senderWidthPart), -1))]
  links_dt[, senderEnd := senderStart + senderWidthPart]
  
  links_dt[, receiverStart := c(0, head(cumsum(receiverWidthPart), -1)), by = receiverPop]
  links_dt[, receiverEnd := receiverStart + receiverWidthPart]
  
  ### Draw Links
  for (i in seq_len(nrow(links_dt))) {
    row <- links_dt[i]
    
    # Skip drawing if both ends are essentially zero to avoid circlize errors
    if (row$senderWidthPart < 1e-6 && row$receiverWidthPart < 1e-6) next
    
    circlize::circos.link(
      sector.index1 = senderSector_v, 
      point1 = c(row$senderStart, row$senderEnd),
      sector.index2 = paste0(row$receiverPop, "_Receiver"),
      point2 = c(row$receiverStart, row$receiverEnd),
      col = grDevices::adjustcolor(arcColors_v[row$Treatment], alpha.f = 0.5),
      directional = 1,
      arr.type = "big.arrow",
      border = NA
    )
  } ### End for
  
  ### Add labels
  if (specialLabel_v) {
    circlize::highlight.sector(sector.index = senderSector_v,
                               track.index = 1,
                               col = "#FFFFFF",
                               text = paste0("Sender: ", senderPop_v, "; Ligand: ", senderDataCol_v),
                               cex = 1, 
                               facing = "bending.inside", 
                               text.col = "black", 
                               niceFacing = T, 
                               text.vjust = -3)
    
    circlize::highlight.sector(#sector.index = receiverSectors_v[ceiling(length(receiverSectors_v)/2)],
                               sector.index = receiverSectors_v[3],
                               track.index = 1,
                               col = "#FFFFFF",
                               text = paste0("          Receiving Clusters; Receptor: ", receiverDataCol_v),
                               cex = 1, 
                               facing = "bending.inside", 
                               text.col = "black", 
                               niceFacing = T, 
                               text.vjust = -6)
  } # fi
  
  graphics::title(main = paste("Links originating from:",senderPop_v, "\n", senderDataCol_v, " to ", receiverDataCol_v), line = -1.5)
  
  if (legendType_v != "none") ComplexHeatmap::draw(lgd, x = unit(0.15, "npc"), y = unit(0.2, "npc"), just = "right")
  
  ### Try adding scale
  addCircosScaleBar(scaleValue = tickScaler_v,
                    totalMax = max(allXmax_v),
                    label = as.character(tickScaler_v),
                    #xPos = -1.4,
                    xPos = -1,
                    yPos = 1)
  
  circlize::circos.clear()
} ### End plotMultiTargetChord