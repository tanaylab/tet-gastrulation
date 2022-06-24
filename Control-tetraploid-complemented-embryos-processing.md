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

# Processing of tetraploid-complemented control embryos 




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


```r
mat1 <- scdb_mat("control_tetra")
```

#### Remove extraembryonic ectoderm and parietal endoderm cells


```r
source("scripts/pipeline/remove_exe_ectoderm_and_parietal_endo_cls.R")
mat1 <- remove_exe_ectoderm_and_parietal_endoderm(mat_query = mat1)
```

#### Merge cells with the wildtype atlas
We first merge this matrix with the tetraploid-complemented control embryos from Mittnenzweig, Mayshar et al. Cell (2021)


```r
source("scripts/pipeline/preprocessing_foxc_experiment_control_tetraploid_mat.R")
preprocessing_tetraploid_control_mat_from_embflow_paper()
```


```r
source("scripts/pipeline/merge_umi_mat_with_wt10_umi_mat.R")
```


```r
source("scripts/pipeline/merge_umi_mat_with_wt10_umi_mat.R")
mat2 <- scdb_mat("control_tetra_from_embflow_paper_f")

mat_new <- merge_two_scmat(mat1 = mat1, mat2 = mat2)

scdb_add_mat("control_tetra_all", mat = mat_new)
```

Next merge matrix with WT cells


```r
mat_new <- scdb_mat("control_tetra_all")
mat_new <- merge_umi_mat_with_wt10(scmat = mat_new)
scdb_add_mat(id = "control_tetra_all_wt10", mat = mat_new)
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
gen_cgraph(mat_nm = "control_tetra_all", force_recompute_cgraph = F)
```

```
## [1] TRUE
```

```r
gen_cgraph(mat_nm = "control_tetra_all_wt10", force_recompute_cgraph = F)
```

```
## [1] TRUE
```

#### Transfer cell type annotation to chimera cells


```r
source("scripts/pipeline/transfer_cell_type_annotation.R")
transfer_color_chimera_tetraploid(mat_nm = "control_tetra_all_wt10")
```

#### Time tetraploid-complemented embryos based on their injected control cells


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
source("scripts/pipeline/tetraploid_timing.R")
tetra_timing("control_tetra_all_wt10")
```

#### Export infered time and cell type into table


```r
source("scripts/pipeline/create_metadata.R")
create_control_tetraploid_metadata()
```

