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


```r
mat <- scdb_mat("tko_chim")
```

#### Remove cells belonging to E8.5 chimera embryos. Only use chimera embryos from E7.5


```r
e85_embryos <- paste0("TKO29_", c(13:19))
ig_cells <- colnames(mat@mat)[mat@cell_metadata[colnames(mat@mat), "embryo"] %in% e85_embryos]
ig_cells <- c(mat@ignore_cells, ig_cells)

mat <- scm_ignore_cells(scmat = mat, ig_cells = ig_cells)
```

#### Remove extraembryonic ectoderm and parietal endoderm cells


```r
source("scripts/pipeline/remove_exe_ectoderm_and_parietal_endo_cls.R")
mat <- remove_exe_ectoderm_and_parietal_endoderm(mat_query = mat)
scdb_add_mat(id = "tko_chim", mat = mat)
```

#### Merge cells with the wildtype atlas


```r
source("scripts/pipeline/merge_umi_mat_with_wt10_umi_mat.R")
mat_new <- merge_umi_mat_with_wt10(scmat = mat)
scdb_add_mat(id = "tko_chim_wt10", mat = mat_new)
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
gen_cgraph(mat_nm = "tko_chim_wt10", force_recompute_cgraph = F)
```

```
## [1] TRUE
```

#### Transfer cell type annotation to chimera cells


```r
source("scripts/pipeline/transfer_cell_type_annotation.R")
transfer_color_chimera_tetraploid(mat_nm = "tko_chim_wt10")
```

#### Time chimera embryos based on their host and KO cells


```r
source("scripts/pipeline/tko_chim_timing.R")
source("scripts/pipeline/transfer_time_annotation.R")
```

```
## Loading required package: slam
```

```
## Loading required package: sparsesvd
```

```r
tko_chim_timing()
# tko_chimera_plot_cor_time_dist_per_emb(mat_nm = "tko_chim_wt10")
```

#### Export infered time and cell type into table


```r
source("scripts/pipeline/create_metadata.R")
create_tko_chimera_metadata()
```

#### Generate a separate UMI matrix and metacell object for the TKO cells


```r
source("scripts/pipeline/generat_ko_mat.R")
generate_ko_mat_chimera()
```


```r
source("scripts/pipeline/generic_mc.R")
tgconfig::override_params(config_file = "config/tet_chim.yaml", package = "metacell")

bad_genes <- read.table("data/tet_tko.bad_genes.txt", sep = "\t", stringsAsFactors = F)
bad_genes <- bad_genes[, 1]

mat_nm <- "tko_chim_ko"

generate_mc(name = mat_nm, add_bad_genes = bad_genes)
```


```r
wt_atlas <- mcell_gen_atlas(
    mat_id = "sing_emb_wt10",
    mc_id = "sing_emb_wt10_recolored",
    gset_id = "sing_emb_wt10",
    mc2d_id = "sing_emb_wt10_recolored"
)
```


```r
cmp <- mcell_proj_on_atlas(
    mat_id = "tko_chim_ko",
    mc_id = "tko_chim_ko",
    atlas = wt_atlas,
    fig_cmp_dir = paste0("figs/atlas_projection_", "tko_chim_ko"), recolor_mc_id = "tko_chim_ko", ten2mars = F
)
```


```r
tgconfig::set_param(param = "mcell_mc2d_T_edge", value = 0.03, package = "metacell")
mcell_mc2d_force_knn(mc2d_id = "tko_chim_ko", mc_id = "tko_chim_ko", graph_id = "tko_chim_ko", graph_parametric = T, feats_gset = "tko_chim_ko")
```

```
## got 74 feat genes for mc graph construction
```

```r
mcell_mc2d_plot("tko_chim_ko")
```

```
## png 
##   2
```


