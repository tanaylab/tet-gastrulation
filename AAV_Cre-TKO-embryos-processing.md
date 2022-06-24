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

# Processing of whole-embryo Tet TKO embryos 
(AAV/Cre-delivered to Tet triple floxed zygotes)




```r
setwd(here::here())
```


```r
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
library("metacell")
library("tgstat")
```


```r
scdb_init("scrna_db", force_reinit = T)
```

```
## initializing scdb to scrna_db
```


```r
mat <- scdb_mat("tko_germline")
```


```r
# for compatibility with later scripts, change column name 'cell_genotype' to 'cell_type'
md <- mat@cell_metadata

if ("cell_genotype" %in% colnames(md)) {
    md <- rename(md, cell_type = cell_genotype)
}

mat@cell_metadata <- md
scdb_add_mat(id = "tko_germline", mat = mat)
```

#### Remove extraembryonic ectoderm and parietal endoderm cells
#### Merge cells with the wildtype atlas


```r
source("scripts/pipeline/remove_exe_ectoderm_and_parietal_endo_cls.R")
source("scripts/pipeline/merge_umi_mat_with_wt10_umi_mat.R")

mat <- remove_exe_ectoderm_and_parietal_endoderm(mat_query = mat)
mat_new <- merge_umi_mat_with_wt10(scmat = mat)
scdb_add_mat(id = "tko_germline_wt10", mat = mat_new)
```

#### Generate single-cell balanced kNN-graph
This step needs the file  
_data/tet_tko.bad_genes.txt_


```r
# 4. generate cgraph
source("scripts/pipeline/gen_cgraph.R")
```

```
## initializing scdb to scrna_db/
```

```r
gen_cgraph(mat_nm = "tko_germline_wt10", force_recompute_cgraph = F)
```

```
## [1] TRUE
```

#### Transfer cell type annotation to TKO cells


```r
source("scripts/pipeline/transfer_cell_type_annotation.R")
transfer_color_chimera_tetraploid("tko_germline_wt10")
```

#### Time single embryos


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
source("scripts/pipeline/tko_germline_timing.R")
tko_germline_timing("tko_germline_wt10")
```


