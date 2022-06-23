---
jupyter:
  jupytext:
    formats: ipynb,Rmd
    text_representation:
      extension: .Rmd
      format_name: rmarkdown
      format_version: '1.2'
      jupytext_version: 1.11.4
  kernelspec:
    display_name: R
    language: R
    name: ir353
---




```r
setwd(here::here())
```


```r
library("metacell")
library("tgstat")
library("tgconfig")
library("Matrix")
library("dplyr")
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```


```r
scdb_init("scrna_db/", force_reinit = T)
```

```
## initializing scdb to scrna_db/
```

#### Remove extraembryonic ectoderm and parietal endoderm cells
#### Merge cells with the wildtype atlas


```r
source("scripts/pipeline/remove_exe_ectoderm_and_parietal_endo_cls.R")
source("scripts/pipeline/merge_umi_mat_with_wt10_umi_mat.R")

for (mat_nm in c("dko12_chim", "dko13_chim", "dko23_chim")) {
    print(mat_nm)
    mat <- scdb_mat(mat_nm)
    mat <- remove_exe_ectoderm_and_parietal_endoderm(mat_query = mat)
    mat_new <- merge_umi_mat_with_wt10(scmat = mat)
    scdb_add_mat(id = paste0(mat_nm, "_wt10"), mat = mat_new)
}
```

```
## [1] "dko12_chim"
## [1] "dko13_chim"
## [1] "dko23_chim"
```

#### Generate single-cell balanced kNN-graph
This step needs the file  
_data/tet_tko.bad_genes.txt_


```r
source("scripts/pipeline/gen_cgraph.R")
```

```
## initializing scdb to scrna_db/
```

```r
gen_cgraph("dko12_chim_wt10", force_recompute_cgraph = F)
```

```
## [1] TRUE
```

```r
gen_cgraph("dko13_chim_wt10", force_recompute_cgraph = F)
```

```
## [1] TRUE
```

```r
gen_cgraph("dko23_chim_wt10", force_recompute_cgraph = F)
```

```
## [1] TRUE
```

#### Transfer cell type annotation to TKO cells


```r
source("scripts/pipeline/transfer_cell_type_annotation.R")
transfer_color_chimera_tetraploid("dko12_chim_wt10")
transfer_color_chimera_tetraploid("dko13_chim_wt10")
transfer_color_chimera_tetraploid("dko23_chim_wt10")
```

#### Time chimera embryos based on their host and KO cells


```r
source("scripts/pipeline/transfer_time_annotation.R")
```

```
## Loading required package: slam
```

```
## Loading required package: sparsesvd
```

```r
source("scripts/pipeline/dko_chim_timing.R")
dko_chim_timing("dko12_chim_wt10")
dko_chim_timing("dko13_chim_wt10")
dko_chim_timing("dko23_chim_wt10")
```


