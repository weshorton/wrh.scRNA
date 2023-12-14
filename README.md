# wrh.scRNA

Collection of helper functions that are used in scRNAseq analyses and lazy-loaded objects containing color assignments, treatment mappings, and other project-specific metadata.

A suffix of \_v means that object is a lazy-loaded vector, a suffix of \_dt refers to a lazy-loaded data.table, and no suffix denotes a function 

# Usage

## Install

```
devtools::install_github("weshorton/wrh.scRNA")
```

## Use Metadata

```
library(wrh.scRNA)

# View major population colors
mPopColors_v

# View mapping between treatment naming conventions
fullAbbrTreatmentMap_dt
```

## View all metadata*

```
> grep("_v|_dt", ls("package:wrh.scRNA", all.names=TRUE), value = T)

 [1] "b12_clusterNames_v"      "b12_lymphoid2Colors_v"  
 [3] "b12_lymphoidColors_v"    "b12_lymphoidPopMap_dt"  
 [5] "b12_myeloid2Colors_v"    "b12_myeloidColors_v"    
 [7] "b12_myeloidPopMap_dt"    "b12_res_v"              
 [9] "b12_treats_v"            "b12_umapNames_v"        
[11] "b3_clusterNames_v"       "b3_lymphoid2Colors_v"   
[13] "b3_lymphoidColors_v"     "b3_lymphoidFullColors_v"
[15] "b3_lymphoidPopMap_dt"    "b3_myeloid2Colors_v"    
[17] "b3_myeloidColors_v"      "b3_myeloidPopMap_dt"    
[19] "b3_res_v"                "b3_treats_v"            
[21] "b3_umapNames_v"          "fpCols_v"               
[23] "fullAbbrTreatmentMap_dt" "gseaColors_v"           
[25] "mPopColors_v"            "oneGroupVolcanoColors_v"
[27] "treatColors_v"           "twoGroupVolcanoColors_v"
[29] "wesSample_scRNATreat_dt"
```

### View all functions*

```
> grep("_v|_dt", ls("package:wrh.scRNA", all.names=TRUE), value = T, invert = T)

 [1] "addVDJ"                  "calcEmbeddings"         
 [3] "calcViolinStatPositions" "checkGenesInData"       
 [5] "checkNames"              "clusterQC"              
 [7] "compareMarkerResults"    "comparePCs"             
 [9] "displaySummary"          "getHeatGenes"           
[11] "getSummaryTables"        "labelSigHeatGenes"      
[13] "leadingEdgeVolcano"      "loadData"               
[15] "loadMsigdbr"             "makePairs"              
[17] "makeRunTable"            "myClusterSweep"         
[19] "myLISI"                  "orderTreatments"        
[21] "plotClusterProps"        "plotGSEA"               
[23] "plotVolcano"             "plotXYScatter"          
[25] "readGlobalMarkers"       "readVsEachMarkers"      
[27] "runFGSEA"                "silhouetteScore"        
[29] "standardPlots"           "vlnPlotCompareMeans"    
[31] "volcanoWrangleMarkers"   "wrangleCollapsePops"
```

*example results are from 2023/12/14
