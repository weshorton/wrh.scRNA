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
