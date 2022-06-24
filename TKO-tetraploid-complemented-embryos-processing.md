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

# Processing of tetraploid-complemented TKO embryos




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
mat <- scdb_mat("tko_tetra")
```

#### Remove extraembryonic ectoderm and parietal endoderm cells


```r
source("scripts/pipeline/remove_exe_ectoderm_and_parietal_endo_cls.R")
mat <- remove_exe_ectoderm_and_parietal_endoderm(mat_query = mat)
```

#### Merge cells with the wildtype atlas


```r
source("scripts/pipeline/merge_umi_mat_with_wt10_umi_mat.R")
mat_new <- merge_umi_mat_with_wt10(scmat = mat)
scdb_add_mat(id = "tko_tetra_wt10", mat = mat_new)
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
gen_cgraph(mat_nm = "tko_tetra_wt10", force_recompute_cgraph = F)
```

```
## [1] TRUE
```

#### Transfer cell type annotation to TKO cells


```r
source("scripts/pipeline/transfer_cell_type_annotation.R")
transfer_color_chimera_tetraploid(mat_nm = "tko_tetra_wt10")
```

#### Time tetraploid-complemented embryos based on their injected TKO cells


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
tetra_timing("tko_tetra_wt10")
```

#### Export infered time and cell type into table


```r
source("scripts/pipeline/create_metadata.R")
create_tko_tetraploid_metadata()
```

#### Generate a separate UMI matrix and metacell object for the TKO cells


```r
source("scripts/pipeline/generat_ko_mat.R")
generate_ko_mat_tetraploid()
```


```r
source("scripts/pipeline/generic_mc.R")
tgconfig::override_params(config_file = "config/tet_chim.yaml", package = "metacell")

bad_genes <- read.table("data/tet_tko.bad_genes.txt", sep = "\t", stringsAsFactors = F)
bad_genes <- bad_genes[, 1]

mat_nm <- "tko_tetra_ko"

generate_mc(mat_nm, add_bad_genes = bad_genes)
```


```r
# Run the line below if you want to transfer cell type annotation from the WT atlas to the tko_tetra_ko metacell object
#
# wt_atlas <- mcell_gen_atlas(
#     mat_id = "sing_emb_wt10",
#     mc_id = "sing_emb_wt10_recolored",
#     gset_id = "sing_emb_wt10",
#     mc2d_id = "sing_emb_wt10_recolored"
# )

# cmp <- mcell_proj_on_atlas(
#     mat_id = "tko_tetra_ko",
#     mc_id = "tko_tetra_ko",
#     atlas = wt_atlas,
#     fig_cmp_dir = paste0("figs/atlas_projection_", "tko_tetra_ko"), recolor_mc_id = "tko_tetra_ko", ten2mars = F
# )
```


```r
tgconfig::set_param(param = "mcell_mc2d_T_edge", value = 0.03, package = "metacell")
mcell_mc2d_force_knn(mc2d_id = "tko_tetra_ko", mc_id = "tko_tetra_ko", graph_id = "tko_tetra_ko", graph_parametric = T, feats_gset = "tko_tetra_ko")
```

```
## got 67 feat genes for mc graph construction
```

```r
mcell_mc2d_plot("tko_tetra_ko")
```

```
## png 
##   2
```




