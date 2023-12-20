#############################################
### Combo 4x scRNAseq Object Descriptions ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################

#' Major Population Colors
#' 
#' Colors used on UMAPs for major populations
#' 
#' @format ## `mPopColors_v`
#' Named vector of length 6
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Population Name}
#' }
"mPopColors_v"

###
### BATCH 1+2 POPULATIONS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

#' Batch1+2 Lymphoid Colors
#' 
#' Colors used on UMAPs for b1+b2 lymphoid populations
#' 
#' @format ## `b12_lymphoidColors_v`
#' Named vector of length 11
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Sub Population Name}
#' }
"b12_lymphoidColors_v"

#' Batch1+2 Lymphoid Colors 2
#' 
#' Colors used on UMAPs for b1+b2 lymphoid populations 
#' when collapsing like sub-pops.
#' 
#' @format ## `b12_lymphoid2Colors_v`
#' Named vector of length 6
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Collapsed Sub Population Name}
#' }
"b12_lymphoid2Colors_v"

#' Batch1+2 Lymphoid Pop Map
#' 
#' Data.table that matches original sPops
#' with collapsed populations
#' 
#' @format ## `b12_lymphoidPopMap_dt`
#' data.table with 11 rows and 2 columns
#' \describe{
#'    \item{sPop}{original sPop labels}
#'    \item{collapsePop}{collapsed population labels}
#' }
"b12_lymphoidPopMap_dt"

#' Batch1+2 Myeloid Colors
#' 
#' Colors used on UMAPs for b1+b2 myeloid populations
#' 
#' @format ## `b12_myeloidColors_v`
#' Named vector of length 7
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Sub Population Name}
#' }
"b12_myeloidColors_v"

#' Batch1+2 Myeloid Colors 2
#' 
#' Colors used on UMAPs for b1+b2 myeloid populations 
#' when collapsing like sub-pops.
#' 
#' @format ## `b12_myeloid2Colors_v`
#' Named vector of length 6
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Collapsed Sub Population Name}
#' }
"b12_myeloid2Colors_v"

#' Batch1+2 Myeloid Pop Map
#' 
#' Data.table that matches original sPops
#' with collapsed populations
#' 
#' @format ## `b12_lymphoidPopMap_dt`
#' data.table with 7 rows and 2 columns
#' \describe{
#'    \item{sPop}{original sPop labels}
#'    \item{collapsePop}{collapsed population labels}
#' }
"b12_myeloidPopMap_dt"

###
### BATCH 3 POPULATIONS  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

#' Batch3 Lymphoid Full Colors
#' 
#' Colors used on UMAPs for b3 lymphoid populations (with Neoplastic)
#' 
#' @format ## `b3_lymphoidFullColors_v`
#' Named vector of length 12
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Sub Population Name}
#' }
"b3_lymphoidFullColors_v"

#' Batch3 Lymphoid Colors
#' 
#' Colors used on UMAPs for b3 lymphoid populations
#' 
#' @format ## `b3_lymphoidColors_v`
#' Named vector of length 11
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Sub Population Name}
#' }
"b3_lymphoidColors_v"

#' Batch3 Lymphoid Colors 2
#' 
#' Colors used on UMAPs for b3 lymphoid populations 
#' when collapsing like sub-pops.
#' 
#' @format ## `b3_lymphoid2Colors_v`
#' Named vector of length 5
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Collapsed Sub Population Name}
#' }
"b3_lymphoid2Colors_v"

#' Batch3 Lymphoid Pop Map
#' 
#' Data.table that matches original sPops
#' with collapsed populations
#' 
#' @format ## `b3_lymphoidPopMap_dt`
#' data.table with 11 rows and 2 columns
#' \describe{
#'    \item{sPop}{original sPop labels}
#'    \item{collapsePop}{collapsed population labels}
#' }
"b3_lymphoidPopMap_dt"

#' Batch3 Myeloid Colors
#' 
#' Colors used on UMAPs for b3 myeloid populations
#' 
#' @format ## `b3_myeloidColors_v`
#' Named vector of length 11
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Sub Population Name}
#' }
"b3_myeloidColors_v"

#' Batch3 Myeloid Colors 2
#' 
#' Colors used on UMAPs for b3 myeloid populations 
#' when collapsing like sub-pops.
#' 
#' @format ## `b3_myeloid2Colors_v`
#' Named vector of length 7
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Collapsed Sub Population Name}
#' }
"b3_myeloid2Colors_v"

#' Batch3 Myeloid Pop Map
#' 
#' Data.table that matches original sPops
#' with collapsed populations
#' 
#' @format ## `b3_lymphoidPopMap_dt`
#' data.table with 11 rows and 2 columns
#' \describe{
#'    \item{sPop}{original sPop labels}
#'    \item{collapsePop}{collapsed population labels}
#' }
"b3_myeloidPopMap_dt"

###
### TREATMENT COLORS AND ASSIGNMENTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

#' Treatment Colors
#' 
#' Colors used on plots for treatments
#' 
#' @format ## `treatColors_v`
#' Named vector of length 9
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Treatment Name}
#' }
"treatColors_v"

#' Batch1+2 Treatment Labels
#' 
#' Treatment names
#' 
#' @format ## `b12_treats_v`
#' vector of length 7
#' \describe{
#'    \item{values}{"4x", "3xR", "3xNR", "2x", "Ent", "PTX", "Unt"}
#' }
"b12_treats_v"

#' Batch3 Treatment Labels
#' 
#' Treatment names
#' 
#' @format ## `b3_treats_v`
#' vector of length 8
#' \describe{
#'    \item{values}{"4x", "3xR", "3xNR", "2x", "Ent", "PTX", "NT", "d80"}
#' }
"b3_treats_v"

###
### BATCH 1+2 OBJECT ATTRIBUTES  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

#' Batch1+2 Cluster names
#' 
#' Column names found in batch12 seurat object
#' metadata that contain the clusters identified
#' from the Seurat analysis pipeline. The names
#' are formated "seurat_clusters_d[dim]_res[res]
#' where dim = dimensions used in reduction and
#' res = resolution used for clustering
#' 
#' @format ## `b12_clusterNames_v`
#' Named vector of length 4
#' \describe{
#'    \item{values}{meta.data column name}
#'    \item{names}{seurat object name}
#' }
"b12_clusterNames_v"

#' Batch1+2 UMAP names
#' 
#' UMAP names found in names() of batch12 seurat object
#' 
#' @format ## `b12_umapNames_v`
#' Named vector of length 4
#' \describe{
#'    \item{values}{umap name from seurat object}
#'    \item{names}{seurat object name}
#' }
"b12_umapNames_v"

#' Batch1+2 Resolution values
#' 
#' Resolution values used for clustering in batch12 seurat object
#' 
#' @format ## `b12_res_v`
#' Named vector of length 4
#' \describe{
#'    \item{values}{resolution}
#'    \item{names}{seurat object name}
#' }
"b12_res_v"

###
### BATCH 3 OBJECT ATTRIBUTES  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

#' Batch3 Cluster names
#' 
#' Column names found in batch3 seurat object
#' metadata that contain the clusters identified
#' from the Seurat analysis pipeline. The names
#' are formated "seurat_clusters_d[dim]_res[res]
#' where dim = dimensions used in reduction and
#' res = resolution used for clustering
#' 
#' @format ## `b3_clusterNames_v`
#' Named vector of length 4
#' \describe{
#'    \item{values}{meta.data column name}
#'    \item{names}{seurat object name}
#' }
"b3_clusterNames_v"

#' Batch3 UMAP names
#' 
#' UMAP names found in names() of batch3 seurat object
#' 
#' @format ## `b3_umapNames_v`
#' Named vector of length 4
#' \describe{
#'    \item{values}{umap name from seurat object}
#'    \item{names}{seurat object name}
#' }
"b3_umapNames_v"

#' Batch3 Resolution values
#' 
#' Resolution values used for clustering in batch3 seurat object
#' 
#' @format ## `b3_res_v`
#' Named vector of length 4
#' \describe{
#'    \item{values}{resolution}
#'    \item{names}{seurat object name}
#' }
"b3_res_v"

###
### MISC STUFF  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

#' WES Sample - Treatment Mapping
#' 
#' data.table that maps sample IDs from Whole Exome
#' Sequencing to treatment labels
#' 
#' @format ## `wesSample_scRNATreat_dt`
#' data.table with 22 rows and 2 columns
#' \describe{
#'    \item{Tumor_Sample_Barcode}{Sample ID in WES data}
#'    \item{Treatment}{Treatment code}
#' }
"wesSample_scRNATreat_dt"

#' scRNA/scATAC Treatment Code Mapping
#' 
#' data.table that maps shorthand treatment
#' codes to long-form "_"-separated codes. Note
#' that PTX_ENT_aCSF1R_unt is in sc
#' 
#' @format ## `fullAbbrTreatmentMap_dt`
#' data.table with 10 rows and 3 columns
#' \describe{
#'    \item{Treatment}{scRNAseq Treatment Code}
#'    \item{longTreatment}{long-form treatment code used in scATAC}
#'    \item{Name}{Full Rx Name}
#' }
"fullAbbrTreatmentMap_dt"

#' Feature Plot Color Scale
#' 
#' Color scale to use in Feature plots
#' blue - green - yellow - orange - red
#' low - med - high
#' Provided by Sushma Nagaraj
#' 
#' @format ## `fpCols_v`
#' vector of hex codes of length 9
#' \describe{
#'    \item{value}{Hex Color}
#'  }
"fpCols_v"

#' GSEA Plot Color Scale
#' 
#' Color scale to use in GSEA plots
#' blue - down
#' red - up
#' 
#' @format ## `gseaColors_v`
#' named vector of color names of length 2
#' \describe{
#'    \item{values}{red, blue}
#'    \item{names}{UP, DOWN}
#'  }
"gseaColors_v"

#' Volcano Color Scale - One group
#' 
#' Color scale to use in standard
#' volcano plots
#' red - higher expression
#' blue - lower expression
#' grey - no change in expression
#' 
#' @format ## `oneGroupVolcanoColors_v`
#' named vector of colors of length 3
#' \describe{
#'    \item{values}{grey, blue, red}
#'    \item{names}(NO, DOWN, UP)
#'  }
"oneGroupVolcanoColors_v"

#' Volcano Color Scale - Two groups
#' 
#' Color scale to use in comparison (two-group)
#' volcano plots
#' red - higher expression in g1
#' blue - lower expression in g1
#' grey - no change in expression in g1
#' darkred - higher expression in g2
#' darkblue - lower expression in g2
#' darkgrey - no change in expression in g2
#' 
#' @format ## `twoGroupVolcanoColors_v`
#' named vector of colors of length 6
#' \describe{
#'    \item{values}{grey, blue, red, darkgrey, darkblue, darkred}
#'    \item{names}{NO_1, DOWN_1, UP_1, NO_2, DOWN_2, UP_2}
#'  }
"twoGroupVolcanoColors_v"

#' Original Data Lymphoid Colors
#' 
#' Colors used on UMAPs for b1+b2 lymphoid populations
#' from original Ana/Rossin analysis
#' 
#' @format ## `orig_lymphoidColors_v`
#' Named vector of length 13
#' \describe{
#'    \item{values}{Hex Color}
#'    \item{names}{Sub Population Name}
#' }
"orig_lymphoidColors_v"
