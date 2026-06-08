#' Add a scale bar legend to a circos plot
#'
#' @description
#' Draws a horizontal scale bar in the margin of an existing circos plot to
#' indicate arc width units. Must be called after circos.track() and before
#' circos.clear().
#'
#' @param scaleValue Numeric scalar. The width value the bar represents.
#' @param totalMax Numeric scalar. The maximum sector xlim value in the plot,
#'   used to calculate the bar's proportional length.
#' @param label Character scalar. Label to display above the bar.
#' @param xPos Numeric scalar. X position of the bar's left edge in canvas coords.
#' @param yPos Numeric scalar. Y position of the bar's bottom edge in canvas coords.
#' @param barHeight Numeric scalar. Height of the bar rectangle.
#' @param col Character scalar. Fill color of the bar.
#'
#' @details
#' Canvas coordinates in circlize default to c(-1, 1) on both axes unless
#' overridden with canvas.xlim/canvas.ylim. Place the bar in the extra space
#' created by expanding the canvas.
#'
#' @return Invisibly NULL. Called for side effects.
#' @export
addCircosScaleBar <- function(scaleValue,
                              totalMax,
                              label = NULL,
                              xPos = -1.3,
                              yPos = -1.25,
                              barHeight = 0.03,
                              col = "grey40") {
  
  ### Calculate bar width as a proportion of the canvas
  ### Canvas spans 2 units (-1 to 1), so scale accordingly
  barWidth <- (scaleValue / totalMax) * 0.5
  
  ### Draw the bar rectangle
  graphics::rect(
    xleft   = xPos,
    ybottom = yPos,
    xright  = xPos + barWidth,
    ytop    = yPos + barHeight,
    col     = col,
    border  = col,
    xpd     = NA
  )
  
  ### Draw end ticks
  graphics::segments(
    x0  = xPos,
    y0  = yPos - 0.02,
    x1  = xPos,
    y1  = yPos + barHeight + 0.02,
    col = col, xpd = NA
  )
  graphics::segments(
    x0  = xPos + barWidth,
    y0  = yPos - 0.02,
    x1  = xPos + barWidth,
    y1  = yPos + barHeight + 0.02,
    col = col, xpd = NA
  )
  
  ### Label above the bar
  if (is.null(label)) {
    label <- as.character(scaleValue)
  } ### End if
  
  graphics::text(
    x      = xPos + barWidth / 2,
    y      = yPos + barHeight + 0.05,
    labels = label,
    cex    = 0.7,
    adj    = c(0.5, 0),
    xpd    = NA
  )
  
  ### Section header below the bar
  graphics::text(
    x      = xPos + barWidth / 2,
    y      = yPos - 0.05,
    labels = "Expression",
    cex    = 0.6,
    adj    = c(0.5, 1),
    col    = "grey30",
    xpd    = NA
  )
  
  invisible(NULL)
} ### End addCircosScaleBar