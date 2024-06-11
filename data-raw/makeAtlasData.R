###
### COMBO 4X TREATMENT AND POPULATION COLORS
###

setwd("/Users/hortowe/my_tool_repos/wrh.scRNA")
library(colorspace)
library(data.table)

###
### B2 Sub-Populations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### B2 Lymphoid
atlas_b2_lymphoidColors_v <- c("NK cells" = "#FFFFB3", "CD4" = "#BEBADA", "CD4 Treg" = "#80B1d3", "Activated T cells" = "#BC80BD", 
                               "CD8" = "#8DD3C7", "B Cells" = "#FB8072", "CD8 Trm" = "#FDB462")

### B2 Myeloid
atlas_b2_myeloidColors_v <- c(`Inflammatory monocytes` = "#E7298A", TAMs = "#66A61E", `TAMs Lyve1+` = "#80B1D3", `TAMs TNF+` = "#D95F02", 
                              `Non-classical monocytes` = "#1B9E77", `Immunosuppressive myeloid` = "#7570B3", DCs = "#E6AB02")

### B2 Neoplastic
atlas_b2_neoplasticColors_v <- c(neo.c0 = "#1B9E77", neo.c1 = "#D95F02", neo.c2 = "#7570B3", neo.c3 = "#E7298A", neo.c4 = "#66A61E")

### B2 Other
atlas_b2_otherColors_v <- c(CAF = "#E99A2C", Endothelial = "#E7B5F5", Myoepithelial = "#91008D")

usethis::use_data(atlas_b2_lymphoidColors_v, overwrite = T)
usethis::use_data(atlas_b2_myeloidColors_v, overwrite = T)
usethis::use_data(atlas_b2_neoplasticColors_v, overwrite = T)
usethis::use_data(atlas_b2_otherColors_v, overwrite = T)

###
### B3 Sub-Populations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### B3 Lymphoid
atlas_b3_lymphoidColors_v <- c(`Cytotoxic CD8 Trm` = "#FDB462", CD4 = "#BEBADA", `B cells` = "#FB8072", `Activated ICOS+ PD1+ T cells` = "#BC80BD", 
                               `DP CD8-CD4 ICOS-CTLA4` = "#B3DE69", `CD4 Treg` = "#80B1d3")

### B3 Myeloid
atlas_b3_myeloidColors_v <- c(`Non-classical monocytes` = "#1B9E77", `Monocytes Trem2+` = "#D95F02", TAMs = "#66A61E", 
                              `TAMs Lyve1+` = "#80B1D3", `Inflammatory monocytes` = "#E7298A", `Immunosuppressive myeloid` = "#7570B3", DCs = "#E6AB02")

### B3 Neoplastic
atlas_b3_neoplasticColors_v <- c(neo.c0 = "#1B9E77", neo.c1 = "#D95F02", neo.c2 = "#7570B3", neo.c3 = "#E7298A", neo.c4 = "#66A61E", 
                                 neo.c5 = "#E6AB02", neo.c6 = "#A6761D", neo.c7 = "#666666", neo.c8 = "#8DD3C7")

### B3 Other
atlas_b3_otherColors_v <- c(CAF = "#E99A2C", Endothelial = "#E7B5F5", Myoepithelial = "#91008D")

usethis::use_data(atlas_b3_lymphoidColors_v, overwrite = T)
usethis::use_data(atlas_b3_myeloidColors_v, overwrite = T)
usethis::use_data(atlas_b3_neoplasticColors_v, overwrite = T)
usethis::use_data(atlas_b3_otherColors_v, overwrite = T)

###
### Gene Sets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

atlasGeneLists_lsv <- list(Neoplastic = c("Anxa6", "Ccl2", "Abcb1a", "Irs1", "Mapk1", "Mapk3"), 
                           MMTvPyMT = c("Sncg", "Htra3", "S100a4", "Ltc4s", "Pparg", "Nnmt", "Ifngr1", "Sparc", "Tgfbi", "Rcn3", "Cd302", "Cd300a", "Atp1a3", "Hpse", "C5ar1"), 
                           Myeloid = "Ccr2")

usethis::use_data(atlasGeneLists_lsv, overwrite = T)