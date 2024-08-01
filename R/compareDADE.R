compareDADE <- function(atacDAG_ls, baseDEG_ls, pop_v, treat1_v = "4x", treat2_v = "3xR", lfc_v = 0.5, pval_v = 0.05,
                                  immunoAntigens_v = NULL, allAntigens_v = NULL) {
  #' Compare Vs Each DEGs with ATAC
  #' @param atacDAG_ls list of DAG results (element name format is paste0(treat1_v, "_vs_", treat2_v))
  #' @param atacGenes_lsv list of atac DA results (one up and one down)
  #' @param baseDEG_ls list of DEG results (e.g. vsEach$batch3_neoplastic or global)
  #' @param pop_v name of population to grab from baseDEG_ls (e.g. c8 for vsEach$batch3_neoplastic$c8, or batch3_neoplastic for global$batch3_neoplastic)
  #' @param treat1_v name of first treatment from comparison
  #' @param treat2_v name of second treatment from comparison
  #' @param lfc_v log fold change cut-off
  #' @param pval_v p value cut-off
  #' @param immunoAntigens_v vector of immunogenic antigens
  #' @param allAntigens_v vector of all antigens
  #' @details atacGenes_lsv must be list in format list("Up" = someVector, "Down" = someVector)
  #' @return tbd
  #' @export
  
  ###
  ### DA/DE Gene overlap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  ### Get up and down genes from atac
  dag_dt <- atacDAG_ls[[paste0(treat1_v, "_vs_", treat2_v)]]
  atacGenes_lsv <- list("Up" = dag_dt[avg_log2FC > 0.5 & p_val_adj < 0.05, Gene],
                        "Down" = dag_dt[avg_log2FC < -0.5 & p_val_adj < 0.05, Gene])
  
  ### Get up and down genes from DEGs
  deg_dt <- baseDEG_ls[[pop_v]][[treat1_v]][[treat2_v]]
  up_v <- deg_dt[avg_log2FC > 0.5 & p_val_adj < 0.05, Gene]
  down_v <- deg_dt[avg_log2FC < -0.5 & p_val_adj < 0.05, Gene]
  
  ### Compare with atac
  upShare_v <- intersect(atacGenes_lsv$Up, up_v)
  downShare_v <- intersect(atacGenes_lsv$Down, down_v)
  
  ### List output
  share_lslsv <- list("Up" = list("All" = upShare_v),
                      "Down" = list("All" = downShare_v))
  
  ###
  ### Antigen overlap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  if (!is.null(immunoAntigens_v)) {
    ### Check with immunogenic antigens
    upShareImmuno_v <- intersect(upShare_v, immunoAntigens_v)
    downShareImmuno_v <- intersect(downShare_v, immunoAntigens_v)
    
    share_lslsv$Up$ImmunoAntigen <- upShareImmuno_v
    share_lslsv$Down$ImmunoAntigen <- downShareImmuno_v
    
  } # fi
  
  if (!is.null(allAntigens_v)) {
    
    if (!is.null(immunoAntigens_v)) {
      otherAntigens_v <- setdiff(allAntigens_v, immunoAntigens_v)
      name_v <- "OtherAntigen"
    } else {
      otherAntigens_v <- allAntigens_v
      name_v <- "AllAntigen"
    } # fi
    
    upShareOther_v <- intersect(upShare_v, otherAntigens_v)
    downShareOther_v <- intersect(downShare_v, otherAntigens_v)
    
    share_lslsv$Up[[name_v]] <- upShareOther_v
    share_lslsv$Down[[name_v]] <- downShareOther_v
    
  } # fi
  
  ###
  ### Outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ###
  
  share_dt <- data.table("Direction" = c(rep("Up", length(share_lslsv$Up)), rep("Down", length(share_lslsv$Down))),
                         "Compare" = rep(names(share_lslsv$Up), 2),
                         "N" = c(sapply(share_lslsv$Up, length), sapply(share_lslsv$Down, length)),
                         "Pop" = rep(pop_v, sum(sapply(share_lslsv, length))))
  
  ### Results output
  rnaRes_dt <- deg_dt[Gene %in% c(upShare_v, downShare_v),]; setorder(rnaRes_dt, "p_val_adj")
  atacRes_dt <- dag_dt[Gene %in% c(upShare_v, downShare_v),]; setorder(atacRes_dt, "p_val_adj")
  
  ### Merge results
  cols_v <- c("Gene", "avg_log2FC", "p_val_adj")
  merge_dt <- merge(rnaRes_dt[,mget(cols_v)], atacRes_dt[,mget(cols_v)],
                    by = "Gene", sort = F, suffixes = c("_RNA", "_ATAC"))
  merge_dt <- merge_dt[,mget(c("Gene", "avg_log2FC_RNA", "avg_log2FC_ATAC", "p_val_adj_RNA", "p_val_adj_ATAC"))]
  
  res_ls <- list("RNA" = rnaRes_dt, "ATAC" = atacRes_dt, "merge" = merge_dt)
  
  ### Return
  return(list("summary" = share_dt,
              "genes" = share_lslsv,
              "res" = res_ls))
  
} # compareVsEachWithATAC