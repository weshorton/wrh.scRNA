vlnPlotCompareMeans <- function(data_df, name_v, indVar_v = "Treatment", measureVar_v, groupBy_v = NULL, method_v = "wilcox.test", sigOnly_v = F, labelTop_v = NULL,
                                padjMethod_v = "holm", displayP_v = "padj", popCol_v = "sPop", comp_dt = NULL, colors_v = NULL, nrow_v = 1, labelSize_v = 5) {
  #' Violin Plot Compare Means
  #' @description
  #' Make a violin plot of some metric and compare groups
  #' @param data_df data.frame containing some combination of groups and measure variables. 
  #' @param name_v name to describe data in data_df. Will be in plot title.
  #' @param indVar_v character vector of independent variable, see details. Treatment is default
  #' @param measureVar_v column name of column in data_dt that contains the metric to be plotted
  #' @param groupBy_v optional value to facet plots by
  #' @param method_v passed to stat::compare_means 'method' argument
  #' @param sigOnly_v logical indicating if signifianct (padj < 0.05) comparisons only should be plotted
  #' @param labelTop_v numeric vector used to limit (or extend) number of comparisons labeled. See details.
  #' @param padjMethod_v passed to stat::compare_means 'p.adjust.method' argument
  #' @param displayP_v vector. If "padj", will display adjusted p-value. if numeric, will round p-value to that many digits.
  #' @param popCol_v either 'sPop' or 'collapsePop'
  #' @param comp_dt data.table of comparisons to subset calculated stats for. If null, will do all. See details.
  #' @param colors_v named vector of hex codes. Names must be values of indVar_v
  #' @param nrow_v number of rows to coerce output into.
  #' @param labelSize_v size of compare means results
  #' @details
  #' indVar_v is used to build a formula with measureVar_v: as.formula(paste0(measureVar_v, " ~ ", indVar_v))  
  #' 
  #' comp_dt has two columns: group1 and group2. The values must be valid values in data_dt[[indVar_v]]  
  #' 
  #' if groupBy_v is null, just one plot, if not, make one plot per unique value in group by and arrange with ggarrange   
  #' 
  
  #' @return violin plot ggplot
  #' @export
  #' 
  
  ### Old details stuff that I don't think I want. New way is to never truncate with labelTop_v
  # If sigOnly_v == T and there are no significant comparisons, if labelTop_v is set, the top N comparisons will be labeled instead
  # If sigOnly_v == T and there are multiple significant comparisons, if labelTop_v is set and is LESS than the number of significant comparisons, they will be truncated
  # If sigOnly_v == F and labelTop_v is set, the full table will be truncated for the top N
  
  ### Build a formula and run compare means
  formula_v <- as.formula(paste0(measureVar_v, " ~ ", indVar_v))
  cm_df <- as.data.table(compare_means(formula = formula_v, data = data_df, group.by = groupBy_v, method = method_v, p.adjust.method = padjMethod_v))
  setorder(cm_df, p.adj)
  
  ### Merge with custom
  if (!is.null(comp_dt)) cm_df <- merge(cm_df, comp_dt, by = c("group1", "group2"), sort = F)
  
  ### Subset for sig, if desired
  if (sigOnly_v) {
    
    ### Subset for significant genes
    sig_df<- cm_df[cm_df$p.adj < 0.05,]
    
    ### If labelTop is set, pad sig_df as needed
    if (!is.null(labelTop_v)) {
      if (nrow(sig_df) < labelTop_v) {
        sig_df <- rbind(sig_df, cm_df[1:(labelTop_v-nrow(sig_df))])
      }
    }
    
    cm_df <- sig_df
  } # fi sig Only
  
  ### Subset for label top
  if (!is.null(labelTop_v) & !sigOnly_v) {
    cm_df <- cm_df[1:labelTop_v]
  }
    
  
  ### Handle colors
  if (is.null(colors_v)) {
    colors_v <- getColors(n_v = length(unique(data_df[[indVar_v]])))
    names(colors_v) <- unique(data_df[[indVar_v]])
  }
  
  ### Get plotting info
  plotInfo_ls <- calcViolinStatPositions(plot_df = data_df, measureVar_v = measureVar_v, cm_df = cm_df,
                                         levels_v = intersect(names(colors_v), unique(data_df[[indVar_v]])))
  cm_df <- plotInfo_ls$cm_df
  yPosBase_v <- plotInfo_ls$yPosBase
  ylim_v <- plotInfo_ls$ylim
  
  if (ylim_v[2] < max(cm_df$yPos)) ylim_v[2] <- max(cm_df$yPos)
  
  ### Make base plot
  if (is.null(groupBy_v)) {
    
    plot_gg <- ggplot2::ggplot(data = data_df, aes(x = !!sym(indVar_v), y = !!sym(measureVar_v))) +
      geom_boxplot(aes(fill = !!sym(indVar_v)), width = 0.9) +
      scale_y_continuous(limits = ylim_v) +
      scale_fill_manual(values = colors_v, breaks = names(colors_v)) +
      wrh.rUtils::my_theme() + wrh.rUtils::angle_x() +
      theme(legend.position = "bottom", axis.text = element_text(size=20), axis.title.y = element_text(size = 20), plot.title = element_text(size=16),
            axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      ylab(measureVar_v) + ggtitle(name_v) +
      stat_compare_means(label = c("p.format"), size = labelSize_v, label.x = 2)
    
    ### Add stats
    if (displayP_v == "padj") {
      plot_gg <- plot_gg +
        ggsignif::geom_signif(data = cm_df, inherit.aes = F, mapping = aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = yPos), manual = T, textsize = labelSize_v)
    } else if (displayP_v == "none") {
      plot_gg <- plot_gg
    } else if (is.logical(all.equal(class(displayP_v), "numeric"))) {
      plot_gg <- plot_gg +
        ggpubr::geom_signif(data = cm_df, inherit.aes = F, mapping = aes(xmin = group1, xmax = group2, 
                                                                         annotations = round(p, digits = displayP_v), y_position = yPos), manual = T, textsize = labelSize_v)
    } # fi displayP_v
    
  } else {
    
    ### One plot per level
    plots_lsgg <- sapply(levels(data_df[[groupBy_v]]), function(x) {
      ### Subset
      subCM_df <- cm_df[cm_df[[groupBy_v]] == x,]
      subData_df <- data_df[data_df[[groupBy_v]] == x,]
      ### Y position
      if (nrow(subCM_df) > 0) subCM_df$yPos <- sapply(seq_along(1:nrow(subCM_df)), function(x) yPosBase_v + (yPosBase_v*.1)*x)
      ### Plot
      plot_gg <- ggplot2::ggplot(data = subData_df, aes(x = !!sym(indVar_v), y = !!sym(measureVar_v))) +
        geom_boxplot(aes(fill = !!sym(indVar_v)), width = 1) +
        scale_y_continuous(limits = ylim_v) +
        scale_fill_manual(values = colors_v, breaks = names(colors_v)) +
        wrh.rUtils::my_theme() + wrh.rUtils::angle_x() +
        theme(legend.position = "bottom", axis.text = element_text(size=20), axis.title.y = element_text(size = 20), plot.title = element_text(size=16), axis.title.x = element_blank(),
              axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        ylab(measureVar_v) + ggtitle(x) +
        stat_compare_means(label = c("p.format"), size = 3, label.x = 2)
      ### Add stats
      if (displayP_v == "padj") {
        plot_gg <- plot_gg +
          ggpubr::geom_signif(data = subCM_df, aes(xmin = group1, xmax = group2, annotations = p.adj, y_position = yPos), manual = T, textsize = 2)
      } else if (is.logical(all.equal(class(displayP_v), "numeric"))) {
        plot_gg <- plot_gg +
          ggpubr::geom_signif(data = subCM_df, aes(xmin = group1, xmax = group2, annotations = round(p, digits = displayP_v), y_position = yPos), manual = T, textsize = 2)
      } # fi displayP_v
      ### Remove y-axis if not 
      if (x != levels(data_df[[popCol_v]])[1]) {
        plot_gg <- plot_gg + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.line.y = element_blank())
      }
      return(plot_gg)}, simplify = F, USE.NAMES = T)
    
    ### Combine
    plot_gg <- ggpubr::ggarrange(plotlist = plots_lsgg, nrow = nrow_v, ncol = ceiling(length(plots_lsgg)/nrow_v), widths = c(2, rep(1, (length(plots_lsgg)-1))), common.legend = T)
    plot_gg <- annotate_figure(p = plot_gg, top = text_grob(label = paste0(name_v), size = 32))
      
  } # fi is.null(groupBy_v)
  
  ### Return
  return(plot_gg)
  
} # vlnPlotCompareMeans
