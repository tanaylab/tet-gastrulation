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

# FACS gating of cell population




```r
setwd(here::here())
```


```r
library("metacell")
library("tgstat")
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
source("scripts/pipeline/proj_on_wt9_atlas.R")
source("scripts/pipeline/add_clone_information_to_scmat_metadata.R")
```


```r
scfigs_init("figs")
scdb_init("scrna_db/", force_reinit = TRUE)
```

```
## initializing scdb to scrna_db/
```

#### Tet TKO chimera embryos


```r
source("scripts/gating/tko_chim_gating.R")
```

#### Tet TKO tetraploid embryos


```r
source("scripts/gating/tko_tetra_gating.R")
```

####  Control tetraploid embryos


```r
source("scripts/gating/control_tetra_gating.R")
```

####  Tet 1/2 DKO chimera embryos


```r
source("scripts/gating/dko12_chim_gating.R")
dko12_gating()
```

####  Tet 1/3 DKO chimera embryos


```r
source("scripts/gating/dko13_chim_gating.R")
dko13_gating()
```

####  Tet 2/3 DKO chimera embryos


```r
source("scripts/gating/dko23_chim_gating.R")
dko23_gating()
```


