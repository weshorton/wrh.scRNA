addVDJ <- function(seurat_obj, type_v, dir_v, barcodeString_v = "1", flag_v = T) {
  #' Add 10x VDJ Data to Seurat Obj
  #' @description Add 10x VDJ T or B data to an initialized GEX Seurat object. Modified from function found:
  #' https://ucdavis-bioinformatics-training.github.io/2020-Advanced_Single_Cell_RNA_Seq/data_analysis/VDJ_Analysis_fixed
  #' @param seurat_obj gene expression seurat object
  #' @param type_v either "t" for TCR or "b" for BCR
  #' @param dir_v directory where "filtered_contig_annotations.csv" and "clonotypes.csv" are found.
  #' @param barcodeString_v see details. CHARACTER string indicating what cell barcodes to grab. 
  #' Default is 1, other values used when loading from multiple samples)
  #' @param flag_v default is True, which will add a column "hasT" or "hasB" (depends on type_v) that will have TRUE if that cell has
  #' TCR/BCR info and FALSE if it does not.
  #' @details Cell barcodes are of the format [ATCG]+\-1 in an individual seurat object. If cellranger aggregate is run on the TCR/BCR output,
  #' then it solves the problem of barcode duplication by changing the "1" at the end with another number. In the case I used to work on this (batch3),
  #' the replaced values are the sample numbers - I'm not sure if this is because of arguments provided in aggregate, 
  #' or if it just goes through numerically and does it. Regardless, if you're unsure take a peak at the barcode suffixes in contigs_dt
  #' @return Exact same seurat object with new metadata columns with clonotype information: 
  #' [type_v]_clonotype_id and [type_v]_cdr3s_aa
  #' @export
  
  ### Load files
  whichDir_v <- ifelse(type_v == "t", "vdj_t", "vdj_b")
  contigs_dt <- fread(file.path(dir_v, whichDir_v, "filtered_contig_annotations.csv"))
  clonotypes_dt <- fread(file.path(dir_v, whichDir_v, "clonotypes.csv"))
  
  ### Subset to get one row for each barcode and just retain clonotype id
  ### So now we have all the cells that had a clonotype and what the id of that clonotype is
  ### Can have the same clonotype id in multiple rows if multiple cells had it
  contigs_dt <- contigs_dt[!duplicated(contigs_dt$barcode), mget(c("barcode", "raw_clonotype_id"))]
  colnames(contigs_dt)[colnames(contigs_dt) == "raw_clonotype_id"] <- "clonotype_id"
  
  ### Merge amino acid sequences from clonotypes table
  contigs_dt <- merge(contigs_dt, clonotypes_dt[,mget(c("clonotype_id", "cdr3s_aa", "cdr3s_nt"))],
                      by = "clonotype_id", sort = F)
  
  ### Reformat with barcodes as rownames for integration with Seurat
  contigs_df <- wrh.rUtils::convertDFT(contigs_dt, col_v = "barcode", rmCol_v = T)
  
  ### Add type designator to column names
  colnames(contigs_df) <- paste(type_v, colnames(contigs_df), sep = "_")
  
  ### Subset for barcode string
  grep_v <- paste0("\\-", barcodeString_v)
  barcodes_v <- grep(grep_v, rownames(contigs_df), value = T)
  subContigs_df <- contigs_df[rownames(contigs_df) %in% barcodes_v,]
  
  ### Replace barcode string with 1
  rownames(subContigs_df) <- stringr::str_replace(rownames(subContigs_df), pattern = barcodeString_v, replacement = "1")
  
  ### Add flags
  if (flag_v) {
    col_v <- paste0("has", toupper(type_v))
    subContigs_df[[col_v]] <- !is.na(subContigs_df[[1]])
  }
  
  ### Add to object
  seurat_obj <- Seurat::AddMetaData(object = seurat_obj, metadata = subContigs_df)
  
  ### Change NA to false if adding flags
  if (flag_v) {
    seurat_obj@meta.data[is.na(seurat_obj@meta.data[[col_v]]), col_v] <- FALSE
  }
  
  ### Return
  return(seurat_obj)
} # addVDJ
