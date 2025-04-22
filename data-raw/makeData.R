###
### COMBO 4X TREATMENT AND POPULATION COLORS
###

setwd("/Users/hortowe/my_tool_repos/wrh.scRNA")
library(colorspace)
library(data.table)
library(AnnotationHub)

### Colors are mainly from Tol Discrete and ditto package

###
### scRNAseq Major populations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

mPopColors_v <- c("#D1BBD7", "#882E72", "#5289C7", "#90C987", "#F7F056", "#F1932D")
names(mPopColors_v) <- c("CAF", "Endothelial", "Myoepithelial", "Lymphoid/NK",
                         "Myeloid", "Neoplastic")

usethis::use_data(mPopColors_v, overwrite = T)

###
### B12 Sub-Populations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### B12 Lymphoid
b12_lymphoidBadColors_v <- c(RColorBrewer::brewer.pal(11, "Set3"))
names(b12_lymphoidBadColors_v) <- c("CD8 T Resident Memory Cells", "NK Cells", "CD4 Tregs", "Memory B Cells", "Effector T Cells",
                                 "Cytotoxic CD8 T Resident Memory Cells", "Proliferating Cytotoxic T Cells", "Cytotoxic CD8 T Cells",
                                 "PD1/CXCR6/ICOS-activated Memory T Cells", "Th2 Cells", "Plasma Cells")

### B12 Lymphoid Collapsed
b12_lymphoidColors_v <- c("#AE76A3", "#1965B0", "#7BAFDE", "#CAE0AB", "#F6C141", "#E8601C")
names(b12_lymphoidColors_v) <- c("NK cells", "B cells", "T cells", "CD8", "CD4", "APC")

### B12 Lymphoid Pop Map
b12_lymphoidPopMap_dt <- data.table("sPop" = c("CD4 Tregs", "Th2 Cells", "CD8 T Resident Memory Cells", "Cytotoxic CD8 T Resident Memory Cells", "Cytotoxic CD8 T Cells",
                                               "Memory B Cells", "NK Cells", "Plasma Cells", "Effector T Cells", "Proliferating Cytotoxic T Cells",
                                               "PD1/CXCR6/ICOS-activated Memory T Cells"),
                                    "collapsePop" = c("CD4", "CD4", "CD8", "CD8", "CD8", "B cells", "NK cells", "APC", "T cells", "T cells", "T cells"))
setkey(b12_lymphoidPopMap_dt, "collapsePop")


### B12 Myeloid
b12_myeloidBadColors_v <- RColorBrewer::brewer.pal(7, "Dark2")
names(b12_myeloidBadColors_v) <- c("Non-classical monocytes", "Recruited monocytes",
                                "Immunosuppressive myeloid", "Inflammatory monocytes",
                                "TAMs", "cDC1s", "Migratory DCs")

### B12 Myeloid Collapsed
b12_myeloidColors_v <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
names(b12_myeloidColors_v) <- c("Inflammatory monocytes", "Immunosuppressive myeloid", "Non-classical monocytes", 
                                 "DCs", "Recruited monocytes", "TAMs")

### B12 Myeloid Pop Map
b12_myeloidPopMap_dt <- data.table("sPop" = c("cDC1s", "Migratory DCs", "Immunosuppressive myeloid", "Inflammatory monocytes", "Non-classical monocytes", "Recruited monocytes", "TAMs"),
                                   "collapsePop" = c("DCs", "DCs", "Immunosuppressive myeloid", "Inflammatory monocytes", "Non-classical monocytes", "Recruited monocytes", "TAMs"))
setkey(b12_myeloidPopMap_dt, "collapsePop")

### B12 Neoplastic
b12_neoplasticColors_v <- c("#AE76A3", "#1965B0", "#7BAFDE", "#CAE0AB", "#F6C141", "#E8601C")
names(b12_neoplasticColors_v) <- paste0("neo.c", 0:4)

### Save all B12
usethis::use_data(b12_lymphoidBadColors_v, overwrite = T)
usethis::use_data(b12_lymphoidColors_v, overwrite = T)
usethis::use_data(b12_lymphoidPopMap_dt, overwrite = T)

usethis::use_data(b12_myeloidBadColors_v, overwrite = T)
usethis::use_data(b12_myeloidColors_v, overwrite = T)
usethis::use_data(b12_myeloidPopMap_dt, overwrite = T)

usethis::use_data(b12_neoplasticColors_v, overwrite = T)

###
### B3 Sub-Populations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### B3 Lymphoid
b3_lymphoidBadColors_v <- c(RColorBrewer::brewer.pal(11, "Set3"), RColorBrewer::brewer.pal(7, "Greys")[7])
names(b3_lymphoidBadColors_v) <- c("CD8 T cells", "Cytotoxic CD8 Trm", "PD1/CXCR6/ICOS-activated Trm (A)",
                                "CD4 Tregs", "CD4 Memory T cells", "B cells", "Dysfunctional T cells",
                                "CD8 Memory T cells", "CD4 Trm", "Proliferating Cytotoxic CD8 T cells",
                                "PD1/CXCR6/ICOS-activated Trm (B)", "Neoplastic")

### B3 Lymphoid Collapsed
b3_lymphoidColors_v <- c("#AE76A3", "#1965B0", "#7BAFDE", "#CAE0AB", "#F6C141")
names(b3_lymphoidColors_v) <- c("CD8", "PD1/CXCR6/ICOS-activated Trm", "B cells", "Dysfunctional T cells", "CD4")

### B3 Lymphoid Pop Map
b3_lymphoidPopMap_dt <- data.table("sPop" = c("B cells", "CD4 Tregs", "CD4 Memory T cells", "CD4 Trm", "CD8 T cells", "Cytotoxic CD8 Trm", "CD8 Memory T cells",
                                              "Proliferating Cytotoxic CD8 T cells", "Dysfunctional T cells", "PD1/CXCR6/ICOS-activated Trm (A)", "PD1/CXCR6/ICOS-activated Trm (B)"),
                                   "collapsePop" = c("B cells", rep("CD4", 3), rep("CD8", 4), "Dysfunctional T cells", rep("PD1/CXCR6/ICOS-activated Trm", 2)))
setkey(b3_lymphoidPopMap_dt, "collapsePop")

### B3 Myeloid
set3_v <- RColorBrewer::brewer.pal(12, "Set3")
set1_v <- RColorBrewer::brewer.pal(9, "Set1")
b3_myeloidBadColors_v <- c(set3_v[1], set1_v[1], set3_v[c(2:3,5:8,4)], set1_v[c(2,3)])
names(b3_myeloidBadColors_v) <- c("Non-classical Monocytes", "Immunosuppressive Myeloid",
                               "Recruited Monocytes", "Inflammatory Monocytes", "Lyve1 Fcgr3 TAMs",
                               "M2-like TAMs", "M1-like TAMs", "Proliferating TAMs", "Resident Macrophages",
                               "cDC1", "DCs")

### B3 Myeloid collapsed
b3_myeloidColors_v <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AD7700")
names(b3_myeloidColors_v) <- c("Resident macrophages", "Recruited monocytes", "Inflammatory monocytes", 
                                "TAMs", "Immunosuppressive myeloid", "DCs", "Non-classical monocytes")

### B3 Myeloid Pop Map
b3_myeloidPopMap_dt <- data.table("sPop" = c("Immunosuppressive Myeloid", "Inflammatory Monocytes", "Non-classical Monocytes", "Recruited Monocytes", "Resident Macrophages", 
                                             "Lyve1 Fcgr3 TAMs", "M2-like TAMs", "M1-like TAMs", "Proliferating TAMs", "cDC1", "DCs"),
                                  "collapsePop" = c("Immunosuppressive myeloid", "Inflammatory monocytes", "Non-classical monocytes", "Recruited monocytes", "Resident macrophages",
                                                    rep("TAMs", 4), rep("DCs", 2)))

### B3 Neoplastic
b3_neoplasticColors_v <- c("#D1BBD7", "#AE76A3", "#882E72", "#1965B0", "#5289C7", "#7BAFDE", "#90C987", "#CAE0AB", "#F7F056", "#F6C141", "#F1932D", "#E8601C", "#DC050C")
names(b3_neoplasticColors_v) <- paste0("neo.c", 0:12)

### Save all B3
usethis::use_data(b3_lymphoidBadColors_v, overwrite = T)
usethis::use_data(b3_lymphoidColors_v, overwrite = T)
usethis::use_data(b3_lymphoidPopMap_dt, overwrite = T)

usethis::use_data(b3_myeloidBadColors_v, overwrite = T)
usethis::use_data(b3_myeloidColors_v, overwrite = T)
usethis::use_data(b3_myeloidPopMap_dt, overwrite = T)

usethis::use_data(b3_neoplasticColors_v, overwrite = T)

###
### TREATMENT COLORS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

treatColors_v <- c("d80" = "#868686", "Unt" = "#000000", "NT" = "#000000", "PTX" = "#F8EC1F", "Ent" = "#488CCB", "2x" = "#EE3224", "3xNR" = "#E09AC4", "3xR" = "#732B90", "4x" = "#F7931E")
b12_treats_v <- c("4x", "3x", "3xR", "3xNR", "2x", "Ent", "PTX", "Unt")
b3_treats_v <- c("4x", "3x", "3xR", "3xNR", "2x", "Ent", "PTX", "NT", "d80")

usethis::use_data(treatColors_v, overwrite = T)
usethis::use_data(b12_treats_v, overwrite = T)
usethis::use_data(b3_treats_v, overwrite = T)

###
### SEURAT OBJECT NAMES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### B12
b12_clusterNames_v <- c("batch12" = "seurat_clusters_d41_res0.6", "lymphoid" = "seurat_clusters_d42_res0.4", 
                        "myeloid" = "seurat_clusters_d19_res0.1", "neoplastic" = "seurat_clusters_d15_res0.1")
b12_umapNames_v <- sapply(b12_clusterNames_v, function(x) gsub("d", "umap", grep("d", unlist(strsplit(x, split = "_")), value = T)))
b12_res_v <- sapply(b12_clusterNames_v, function(x) as.numeric(gsub("r\\.|res", "", grep("r\\.|res", unlist(strsplit(x, split = "_")), value = T))))


### B3
b3_clusterNames_v <- c("batch3" = "seurat_clusters_d50_r.5", "lymphoid" = "seurat_clusters_d16_res0.5", 
                       "myeloid" = "seurat_clusters_d42_res0.3", "neoplastic" = "seurat_clusters_d42_res0.2")
b3_umapNames_v <- sapply(b3_clusterNames_v, function(x) gsub("d", "umap", grep("d", unlist(strsplit(x, split = "_")), value = T)))
b3_res_v <- sapply(b3_clusterNames_v, function(x) as.numeric(gsub("r\\.|res", "", grep("r\\.|res", unlist(strsplit(x, split = "_")), value = T))))

### Output
usethis::use_data(b12_clusterNames_v, overwrite = T)
usethis::use_data(b12_umapNames_v, overwrite = T)
usethis::use_data(b12_res_v, overwrite = T)

usethis::use_data(b3_clusterNames_v, overwrite = T)
usethis::use_data(b3_umapNames_v, overwrite = T)
usethis::use_data(b3_res_v, overwrite = T)

###
### Misc Stuff ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Mapping of WES samples names to scRNAseq treatment names
wesSample_scRNATreat_dt <- data.table("Tumor_Sample_Barcode" = c("B5", "B6", "C1", "C3", "D4", "D5", "D6", "H4", "H5", "H6",
                                                                 "E4", "E5", "E6", "F1", "F4", "F5", "G4", "G5", "G6", "I4", "I5", "I6"),
                                      "Treatment" = c(rep("d80", 2), rep("NT", 2), rep("PTX", 3), rep("Ent", 3), rep("2x", 3), 
                                                      rep("3xNR", 3), rep("3xR", 3), rep("4x", 3)))

usethis::use_data(wesSample_scRNATreat_dt, overwrite = T)

### Mapping of explicit treatment names
fullAbbrTreatmentMap_dt <- data.table("Treatment" = c("d80", "NT", "Unt", "PTX", "Ent", "2x", "3xNR", "3xR", "4x", "3xENT"),
                                      "longTreatment" = c(NA, "unt_unt_unt_unt", "unt_unt_unt_unt", "PTX_unt_unt_unt", 
                                                          "unt_ENT_unt_unt", "PTX_unt_aCSF1R_unt", "PTX_unt_aCSF1R_aPD1N", 
                                                          "PTX_unt_aCSF1R_aPD1R", "PTX_ENT_aCSF1R_aPD1R", "PTX_ENT_aCSF1R_unt"),
                                      "Name" = c("Untreated d80", "Untreated", "Untreated", "Paclitaxel", "Entinostat", 
                                                 "PTX/aCSF1R", "PTX/aCSF1R/aPD-1 Non-Regressors", "PTX/aCSF1R/aPD-1 Regressors", 
                                                 "PTX/aCSF1R/aPD-1/Entinostat", "PTX/aCSF1R/Entinostat"))

fullAbbrTreatmentMap_dt <- merge(fullAbbrTreatmentMap_dt, data.table("Treatment" = names(treatColors_v), "HexColor" = treatColors_v),
                                 by = "Treatment", all = T, sort = F)

usethis::use_data(fullAbbrTreatmentMap_dt, overwrite = T)

### Plot Colors
fpCols_v = c("#2C7BB6", "#00A6CA","#00CCBC","#90EB9D","#FFFF8C","#F9D057","#F29E2E","#E76818","#D7191C") # feature
gseaColors_v <- c("UP" = "red", "DOWN" = "blue") # gsea
oneGroupVolcanoColors_v <- c("NO" = "grey", "DOWN" = "blue", "UP" = "red")
twoGroupVolcanoColors_v <- c("NO_1" = "grey", "DOWN_1" = "blue", "UP_1" = "red", "NO_2" = "darkgrey", "DOWN_2" = "darkblue", "UP_2" = "darkred")


usethis::use_data(fpCols_v, overwrite = T)
usethis::use_data(gseaColors_v, overwrite = T)
usethis::use_data(oneGroupVolcanoColors_v, overwrite = T)
usethis::use_data(twoGroupVolcanoColors_v, overwrite = T)

###
### Gene Lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Main Markers
tumors <- c("Epcam","Krt8","Krt18","Krt14","Cldn4","Cldn3","Aldh1a3","Kit","Rrm2")
caf <- c("Mfap5","S100a4","Fap","Itgb1","Acta2","Mme","Clu","Pdgfra","Pdgfrb","Sparc","Vim","Col3a1","Gsn","Dcn")
myoepi <- c("Krt17","Fst","Sphk1","Actg2","Serpinb5","Krt5","Myh11","Cnn1","Pxn","Lamb3","Lamc2","Mdfi","Grwd1","Smim3","Mylk","Oxtr","Krt4")
endo <- c("Pecam1","Fabp4","Aqp7","Cdh5","S1pr1","Prox1","Lyve1","Pdpn","Esam")
mono <- c("Cx3cr1","Ccr2","Trem1", "Trem2", "Lyz2")
dc <- c("Xcr1","Itgax","Itgae","Ftl1","Btla","H2-Ab1","Ccr7")
tam <- c("Adgre1","Ccr5","Il6","Il1b","Cd86", "Mrc1","Il10","Cd163","Arg1")
mdsc <- c("Itgam", "Plac8","Cd14","Cd84","Tnf","Ly6c1", "S100a8","S100a9", "Ly6g")
tcells <- c("Cd3e","Cd3d","Cd3g","Cd8a","Cd4","Foxp3","Il2ra","Ctla4","Pdcd1","Icos")
bcells <- c("Cd79a", "Cd79b","Ms4a1","Cd19")
nk <- c("Nkg7", "Klrb1c","Gzmb","Ncr1")

### Combine
markers_lsv <- list("tumor" = tumors, "Caf" = caf, "MyoEpi" = myoepi, "Endo" = endo, "Mono" = mono,
                    "DC" = dc, "TAMs" = tam, "MDSC" = mdsc, "Tcells" = tcells, "Bcells" = bcells, "NK" = nk)

### Number of columns for each feature plot set
markerCols_lsv <- list("tumor" = 3, "Caf" = 3, "MyoEpi" = 3, "Endo" = 3, "Mono" = 3,
                 "DC" = 3, "TAMs" = 3, "MDSC" = 3, "Tcells" = 3, "Bcells" = 2, "NK" = 2)

### Get vector of names
names_v <- unlist(lapply(names(markers_lsv), function(x) rep(x, length(markers_lsv[[x]]))))

### Combine all the markers and add names
dotMarkers_v <- unique(unlist(markers_lsv))
names(dotMarkers_v) <- names_v

usethis::use_data(markers_lsv, overwrite = T)
usethis::use_data(markerCols_lsv, overwrite = T)
usethis::use_data(dotMarkers_v, overwrite = T)

###
### Markers for figures
###

### Major Pop Dotplot
majorPopGenes_v <- c("Epcam", "Krt18", "Krt8", "Itgam","Itgax", "Ly6g", "Adgre1", "Arg1",
                     "Cd3e", "Cd8a", "Cd4", "Ms4a1", "Cd79a", "Esam", "S1pr1", "Pecam1", 
                     "Mylk", "Myh11", "Col3a1", "Dcn", "Gsn", "Fap", "Mfap5")

### Myeloid Minor Pop Dotplot
myeloidMinorPopGenes_lsv <- list("TAMs" = c("Itgam", "Adgre1", "Trem2", "Il10", "Cx3cr1"),
                                 "Recruited" = c("Ftl1", "Itgax", "Cx3cr1"),
                                 "Immunosuppressive" = c("S100a8", "S100a9", "Il1b"),
                                 "Inflammatory" = c("Plac8", "Ccr5", "Ccr2", "Itgam"),
                                 "ncTam" = c("Adgre1", "Itgam", "Mrc1", "Cd163"),
                                 "DC"= c("Cxcr1", "Batf3", "Btla", "Itgae"))

myeloidMinorPopGenes_v <- unique(unlist(myeloidMinorPopGenes_lsv))

### Lymphoid Minor Pop Dotplot
lymphoidMinorPopGenes_lsv <- list("NK" = c("Ncr1", "Klrb1c", "Nkg7"),
                                  "T cells" = c("Cd3g", "Cd3e", "Cd4", "Cd8a"),
                                  "CD8" = "Cd3g", "Cd3e", "Cd8a", "Gzmb",
                                  "CD4" = "Cd3g", "Cd3e", "Cd4",
                                  "B cells" = c("Cd19", "Ms4a1", "Cd79a", "Cd79b", "Cd86", "H2-Ab1"),
                                  "Plasma cells" = c("Cd19", "Ms4a1", "Cd79a", "Cd79b"))
lymphoidMinorPopGenes_v <- unique(unlist(lymphoidMinorPopGenes_lsv))

usethis::use_data(majorPopGenes_v, overwrite = T)
usethis::use_data(myeloidMinorPopGenes_v, overwrite = T)
usethis::use_data(lymphoidMinorPopGenes_v, overwrite = T)

###
### Cell Cycle Genes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###

### Load genes and view
cell_cycle_genes <- read.csv(text = RCurl::getURL("https://raw.githubusercontent.com/hbc/tinyatlas/refs/heads/master/cell_cycle/Mus_musculus.csv"))

### Get ensembl data
### Grab gene name - gene ID mapping from ensembl
ah <- AnnotationHub()                                                         # Connect to hub
ahDb <- query(ah, pattern = c("Mus musculus", "EnsDb"), ignore.case = TRUE)   # Grab mouse db
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)                         # Get latest annotations
edb <- ah[[id, force=T]]                                                      # Download db
annotations <- genes(edb, return.type = "data.frame")                         # Get genes
annotations <- annotations %>%                                                # Select annotations
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

### Combine with cell-cycle marker IDs and split
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
s_genes <- cell_cycle_markers %>% dplyr::filter(phase == "S") %>% dplyr::pull("gene_name")
g2m_genes <- cell_cycle_markers %>% dplyr::filter(phase == "G2/M") %>% dplyr::pull("gene_name")

usethis::use_data(s_genes, overwrite = T)
usethis::use_data(g2m_genes, overwrite = T)
