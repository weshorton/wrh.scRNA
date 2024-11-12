compareVsEachWithATAC <- function(atacGenes_lsv, baseDEG_ls, pop_v, treat1_v = "4x", treat2_v = "3xR", lfc_v = 0.5, pval_v = 0.05,
                                  immunoAntigens_v, allAntigens_v) {
  #' Compare Vs Each DEGs with ATAC
  #' @param atacGenes_lsv list of atac DA results (one up and one down)
  #' @param baseDEG_ls list of DEG results (e.g. vsEach$batch3_neoplastic)
  #' @param pop_v name of population to grab from baseDEG_ls (e.g. c8 for vsEach$batch3_neoplastic$c8)
  #' @param treat1_v name of first treatment from comparison
  #' @param treat2_v name of second treatment from comparison
  #' @param lfc_v log fold change cut-off
  #' @param pval_v p value cut-off
  #' @param immunoAntigens_v vector of immunogenic antigens
  #' @param allAntigens_v vector of all antigens
  #' @details atacGenes_lsv must be list in format list("Up" = someVector, "Down" = someVector)
  #' @return tbd
  #' @export
  
  ### Get up and down genes from DEGs
  deg_dt <- baseDEG_ls[[pop_v]][[treat1_v]][[treat2_v]]
  up_v <- deg_dt[avg_log2FC > 0.5 & p_val_adj < 0.05, Gene]
  down_v <- deg_dt[avg_log2FC < -0.5 & p_val_adj < 0.05, Gene]
  
  ### Compare with atac
  upShare_v <- intersect(atacGenes_lsv$Up, up_v)
  downShare_v <- intersect(atacGenes_lsv$Down, down_v)
  
  ### Check with immunogenic antigens
  upShareImmuno_v <- intersect(upShare_v, immunoAntigens_v)
  downShareImmuno_v <- intersect(downShare_v, immunoAntigens_v)
  
  ### Check with other antigens
  otherAntigens_v <- setdiff(allAntigens_v, immunoAntigens_v)
  upShareOther_v <- intersect(upShare_v, otherAntigens_v)
  downShareOther_v <- intersect(downShare_v, otherAntigens_v)
  
  ### List output
  share_lslsv <- list("Up" = list("All" = upShare_v, "ImmunoAntigen" = upShareImmuno_v, "OtherAntigen" = upShareOther_v),
                      "Down" = list("All" = downShare_v, "ImmunoAntigen" = downShareImmuno_v, "OtherAntigen" = downShareOther_v))
  
  ### Table Output
  share_dt <- data.table("Direction" = c(rep("Up", 3), rep("Down", 3)),
                         "Compare" = rep(c("All", "ImmunoAntigen", "OtherAntigen"), 2),
                         "N" = c(length(upShare_v), length(upShareImmuno_v), length(upShareOther_v),
                                 length(downShare_v), length(downShareImmuno_v), length(downShareOther_v)),
                         "Pop" = rep(pop_v, 6))
  
  ### Return
  return(list("summary" = share_dt,
              "genes" = share_lslsv))
  
} # compareVsEachWithATAC