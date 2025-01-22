wrangleDEG <- function(deg_dt, subList_lsv = list("neoantigen" = allNeoantigens_v, 
                                                  "immunogenic" = allImmunoNeoantigens_v, 
                                                  "top65Immuno" = sortAllTop65ImmunoNeoantigens_v, 
                                                  "atac4x" = atacNeoantigenList_v)) {
  #' Wrangle DEG
  #' @description Wrangle DEGs with neoantigen labels and subset
  #' @param deg_dt data.table converted from FindMarkers data.frme output
  #' @param subList_lsv list of genes to subset data by
  #' @details This is used in the neoplasticDEG_checkAntigen reports. the referenced antigens in subList_lsv are loaded there.
  #' @return list of results subset variously
  #' @export
  
  ### Add columns
  for (i in 1:length(subList_lsv)) {
    
    ### Get info
    currColumnName_v <- names(subList_lsv)[i]
    currAntigens_v <- subList_lsv[[currColumnName_v]]
    
    ### Add column
    deg_dt[[currColumnName_v]] <- F
    deg_dt[Gene %in% currAntigens_v, currColumnName_v] <- T
    
  } # for i
  
  ### Split deg into up and down (and significant)
  up_dt <- deg_dt[avg_log2FC > 0 & p_val_adj < 0.05,]; up_lsdt <- list()
  dn_dt <- deg_dt[avg_log2FC < 0 & p_val_adj < 0.05,]; dn_lsdt <- list()
  
  up_lsdt <- upGenes_lsv <- dn_lsdt <- dnGenes_lsv <- list()
  for (i in 1:length(subList_lsv)) {
    
    ### Get info
    currColumnName_v <- names(subList_lsv)[i]
    currAntigens_v <- subList_lsv[[currColumnName_v]]
    
    ### Subset
    up_lsdt[[currColumnName_v]] <- up_dt[get(currColumnName_v) == T,]
    dn_lsdt[[currColumnName_v]] <- dn_dt[get(currColumnName_v) == T,]
    
    ### Extract Genes
    upGenes_lsv[[paste0("Sig Up\n", currColumnName_v)]] <- up_lsdt[[currColumnName_v]]$Gene
    dnGenes_lsv[[paste0("Sig Down\n", currColumnName_v)]] <- dn_lsdt[[currColumnName_v]]$Gene
    
  } # for i
  
  ### Output
  out_lsls <- list("dt" = list("up" = up_lsdt, "dn" = dn_lsdt),
                   "genes" = list("up" = upGenes_lsv, "dn" = dnGenes_lsv))
  
  return(out_lsls)
  
} # wrangle DEG