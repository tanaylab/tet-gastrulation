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
scdb_init("scrna_db/", force_reinit = TRUE)
```

```
## initializing scdb to scrna_db/
```

### Import Tet TKO chimera plates into metacell
This notebook requires that the UMI tables for each plate are available in the folder **data/umi.tables** . If they are not available, please run the cell below to download the necessary data


```r
# If you need to download the umi tables, please run the two lines below
if (0) {
    source("../scripts/download_data.R")
    download_umi_tables()
}
```


```r
source("scripts/pipeline/import_plates.R")
```


```r
mat_nm <- "tko_chim"

metadata <- read.csv(file = "data/plate_metadata.csv", stringsAsFactors = F)
genotype <- "TKO"
treatment <- c("Chimera assay")
metadata <- metadata[metadata$treatment %in% treatment, ]
metadata <- metadata[grep(pattern = genotype, x = metadata$genotype), ]

# remove additionally the plate 190730_P03. No embryo and facs data available. Only relevant for TKO Chimera
metadata <- metadata[metadata$plate != "190730_P03", ]

import_plates(mat_nm = mat_nm, metadata = metadata)
```

```
## will read 190303_P02
```

```
## will read 190318_P01
```

```
## will read 190318_P02
```

```
## will read 190318_P03
```

```
## will read 190318_P04
```

```
## will read 190318_P05
```

```
## will read 190318_P06
```

```
## will read 190318_P07
```

```
## will read 190318_P09
```

```
## will read 190318_P11
```

```
## will read 190318_P12
```

```
## will read 190322_P01
```

```
## will read 190322_P02
```

```
## will read 190322_P04
```

```
## will read 190325_P01
```

```
## will read 190325_P02
```

```
## will read 190325_P03
```

```
## will read 190325_P04
```

```
## will read 190325_P05
```

```
## will read 190325_P07
```

```
## will read 190325_P08
```

```
## will read 190730_P01
```

```
## will read 190730_P02
```

```
## will read 190730_P04
```

```
## will read 190730_P05
```

```
## will read 190730_P06
```

```
## will read 190730_P07
```

```
## will read 190730_P08
```

```
## will read 190730_P09
```

```
## will read 190730_P10
```

```
## will read 190730_P11
```

```
## will read 190806_P04
```

```
## will read 190806_P05
```

```
## will read 190806_P06
```

```
## will read 20190806_P01
```

```
## will read 20190806_P02
```

```
## will read 20190806_P03
```

```
## will read 20190806_P07
```

```
## will read 20190806_P08
```

```
## will read 20190806_P09
```

```
## will read 1106TKO29_01
```

```
## will read 1106TKO29_02
```

```
## will read 1106TKO29_03
```

```
## will read 1106TKO29_04
```

```
## will read 1106TKO29_05
```

```
## will read 1106TKO29_06
```

```
## will read 1106TKO29_07
```

```
## will read 1106TKO29_08
```

```
## will read 1106TKO29_09
```

```
## will read 1106TKO29_10
```

```
## will read 1106TKO29_11
```

```
## will read 1106TKO29_12
```

```
## will read 1106TKO29_13
```

```
## will read 1106TKO29_14
```

```
## will read 1106TKO29_15
```

```
## will read 191108_P01
```

```
## will read 191108_P02
```

```
## will read 191108_P03
```

```
## will read 191108_P04
```

```
## will read 191108_P05
```

```
## will read 191108_P06
```

```
## will read 191108_P07
```

```
## will read 191108_P08
```

```
## will read 20200106_TKO26_P01_NOVA
```

```
## will read 20200106_TKO26_P02_NOVA
```

```
## will read 20200106_TKO26_P03_NOVA
```

```
## will read 20200106_TKO26_P04_NOVA
```

```
## will read 20200106_TKO26_P05_NOVA
```

```
## will read TKO26_13_P01
```

```
## will read TKO26_13_P02
```

```
## will read TKO26_13_P03
```

```
## will read TKO26_13_P04
```

```
## will read TKO26_13_P05
```

```
## will read TKO26_13_P06
```

```
## will read TKO26_13_P07
```

```
## will read TKO26_13_P08
```

```
## will read TKO26_13_P09
```

```
## will read TKO26_13_P10
```

```
## will read TKO26_13_P11
```

```
## will read TKO26_13_P12
```

```
## will read TKO26_13_P15
```

```
## will read TetTKO23_P01
```

```
## will read TetTKO23_P02
```

```
## will read TetTKO23_P03
```

```
## will read TetTKO23_P04
```

```
## will read TetTKO23_P05
```

```
## will read TetTKO23_P06
```

```
## will read TetTKO23_P07
```

```
## will read TetTKO23_P08
```

```
## will read TKO23_P01
```

```
## will read TKO23_P02
```

```
## will read TKO23_P04
```

```
## will read TKO23_P05
```

```
## will read TKO23_P06
```

```
## [1] TRUE
```

### Import Tet TKO tetraploid-complemented embryos


```r
mat_nm <- "tko_tetra"

metadata <- read.csv(file = "data/plate_metadata.csv", stringsAsFactors = F)
genotype <- "TKO"
treatment <- c("Tetraploid complementation assay")
metadata <- metadata[metadata$treatment %in% treatment, ]
metadata <- metadata[grep(pattern = genotype, x = metadata$genotype), ]

import_plates(mat_nm = mat_nm, metadata = metadata)
```

```
## will read TKO26_tetra_191111_P01_NOVAseq
```

```
## will read TKO26_tetra_191111_P02
```

```
## will read TKO26_tetra_191111_P03
```

```
## will read TKO26_tetra_191111_P04
```

```
## will read 191125_P01_NOVA
```

```
## will read 191125_P02
```

```
## will read 191125_P03
```

```
## will read 191125_P04
```

```
## will read 4N_P01
```

```
## will read 4N_P06
```

```
## will read 4N26_01
```

```
## will read 4N26_02
```

```
## will read 4N26_03
```

```
## will read 4N26_04
```

```
## will read 4N26_05
```

```
## will read 4N26_06
```

```
## will read 4N26_07
```

```
## will read 4N26_08
```

```
## will read 4N_P07
```

```
## will read 190819_P01
```

```
## will read 190819_P02
```

```
## will read 190819_P03
```

```
## will read 190819_P04
```

```
## will read 190916_P07
```

```
## will read 190916_P08
```

```
## will read 190916_P09
```

```
## will read 4N_P04
```

```
## will read 4N_P05
```

```
## will read 4N_P08
```

```
## will read 4N_P10
```

```
## will read 4N_P11
```

```
## will read 4N_P12
```

```
## will read 4N_P13
```

```
## will read 4N_P14
```

```
## will read 4N29_01
```

```
## will read 4N29_02
```

```
## will read 4N29_03
```

```
## will read 4N29_04
```

```
## will read 4N29_05
```

```
## will read 4N29_06
```

```
## will read 4N29_07
```

```
## will read 4N29_08
```

```
## will read 4N29_09
```

```
## will read 4N29_10
```

```
## will read 4N29_11
```

```
## will read 4N29_12
```

```
## will read 4N29_13
```

```
## will read 4N29_14
```

```
## [1] TRUE
```


```r
# There are four mixed plates containing both tetraploid-complemented control and Tet TKO embryos
# We additionally ignore the cells from control embryos
mat <- scdb_mat("tko_tetra")

ctrl_cells <- mat@cell_metadata$cell[grep(pattern = "Ctrl", x = mat@cell_metadata$embryo)]

mat <- scm_ignore_cells(scmat = mat, ig_cells = union(mat@ignore_cells, ctrl_cells))

scdb_add_mat(id = "tko_tetra", mat = mat)
```

### Import Control tetraploid-complemented embryos


```r
mat_nm <- "control_tetra"

metadata <- read.csv(file = "data/plate_metadata.csv", stringsAsFactors = F)
genotype <- "Control"
treatment <- c("Tetraploid complementation assay")
metadata <- metadata[metadata$treatment %in% treatment, ]
metadata <- metadata[grep(pattern = genotype, x = metadata$genotype), ]

import_plates(mat_nm = mat_nm, metadata = metadata)
```

```
## will read 191125_P04
```

```
## will read 190819_P01
```

```
## will read 190916_P07
```

```
## will read 190916_P09
```

```
## will read 190916_P01
```

```
## will read 190916_P02
```

```
## will read 190916_P03
```

```
## will read 190916_P04
```

```
## will read 190916_P05
```

```
## will read 190916_P06
```

```
## will read 0217_P01
```

```
## will read 0217_P02
```

```
## will read 0217_P03
```

```
## will read 0217_P04
```

```
## will read 0217_P05
```

```
## will read 191125_P05
```

```
## will read 191125_P06_NOVA
```

```
## will read 191125_P07
```

```
## [1] TRUE
```


```r
# There are four mixed plates containing both tetraploid-complemented control and Tet TKO embryos
# We additionally ignore the cells from TKO embryos
mat <- scdb_mat("control_tetra")

tko_cells <- mat@cell_metadata$cell[grep(pattern = "TKO", x = mat@cell_metadata$embryo)]
length(tko_cells)
```

```
## [1] 809
```

```r
mat <- scm_ignore_cells(scmat = mat, ig_cells = union(mat@ignore_cells, tko_cells))

scdb_add_mat(id = "control_tetra", mat = mat)
```

### Import Tet 1/2 DKO chimera embryos


```r
mat_nm <- "dko12_chim"

metadata <- read.csv(file = "data/plate_metadata.csv", stringsAsFactors = F)
genotype <- "Tet1/2 DKO"
treatment <- c("Chimera assay")
metadata <- metadata[metadata$treatment %in% treatment, ]
metadata <- metadata[grep(pattern = genotype, x = metadata$genotype), ]

import_plates(mat_nm = mat_nm, metadata = metadata)
```

```
## will read Tet12_0228_P01
```

```
## will read Tet12_0228_P02
```

```
## will read Tet12_0228_P03
```

```
## will read Tet12_0228_P04
```

```
## will read Tet12_0228_P05
```

```
## will read Tet12_0228_P06
```

```
## will read Tet12_0228_P07
```

```
## will read Tet12_0228_P08
```

```
## will read Tet12_0228_P09
```

```
## will read Tet12_0228_P10
```

```
## will read Tet12_0303_P01
```

```
## will read Tet12_0303_P02
```

```
## will read Tet12_0303_P03
```

```
## will read Tet12_0303_P04
```

```
## will read Tet12_1_0616
```

```
## will read Tet12_2_0616
```

```
## will read Tet12DKO_10
```

```
## will read Tet12DKO_11
```

```
## will read Tet12DKO_12
```

```
## will read Tet12DKO_13
```

```
## will read Tet12DKO_14
```

```
## will read Tet12DKO_15
```

```
## will read Tet12DKO_7
```

```
## will read Tet12DKO_8
```

```
## will read Tet12DKO_9
```

```
## [1] TRUE
```

### Import Tet 1/3 DKO chimera embryos


```r
mat_nm <- "dko13_chim"

metadata <- read.csv(file = "data/plate_metadata.csv", stringsAsFactors = F)
genotype <- "Tet1/3 DKO"
treatment <- c("Chimera assay")
metadata <- metadata[metadata$treatment %in% treatment, ]
metadata <- metadata[grep(pattern = genotype, x = metadata$genotype), ]

import_plates(mat_nm = mat_nm, metadata = metadata)
```

```
## will read mTet13_P01
```

```
## will read mTet13_P03
```

```
## will read mTet13_P07
```

```
## will read mTet13_P08
```

```
## will read Tet13_0309_P02
```

```
## will read Tet13_0309_P03
```

```
## will read Tet13_0309_P04
```

```
## will read Tet13_0309_P05
```

```
## will read Tet13_2_0613
```

```
## will read Tet13_DKO_3
```

```
## will read Tet13_DKO_4
```

```
## will read Tet13_P05
```

```
## will read Tet13_P06
```

```
## [1] TRUE
```

### Import Tet 2/3 DKO chimera embryos


```r
mat_nm <- "dko23_chim"

metadata <- read.csv(file = "data/plate_metadata.csv", stringsAsFactors = F)
genotype <- "Tet2/3 DKO"
treatment <- c("Chimera assay")
metadata <- metadata[metadata$treatment %in% treatment, ]
metadata <- metadata[grep(pattern = genotype, x = metadata$genotype), ]

import_plates(mat_nm = mat_nm, metadata = metadata)
```

```
## will read Tet23_0303_P01
```

```
## will read Tet23_0303_P02
```

```
## will read Tet23_0303_P03
```

```
## will read Tet23_1
```

```
## will read Tet23_2
```

```
## will read Tet23_4
```

```
## will read Tet23_5
```

```
## will read Tet23_7
```

```
## will read Tet23_8
```

```
## will read Tet23_DKO_1_0616
```

```
## will read Tet23_DKO_2_0616
```

```
## will read Tet23DKO_3_0616
```

```
## will read Tet23DKO_4_0616
```

```
## [1] TRUE
```

### Import AAV-Cre delivered whole-embryo TKOs


```r
mat_nm <- "tko_germline"

metadata <- read.csv(file = "data/plate_metadata.csv", stringsAsFactors = F)
genotype <- "TKO"
treatment <- c("AAV Cre-delivered Tet1/2/3 triple KO")
metadata <- metadata[metadata$treatment %in% treatment, ]
metadata <- metadata[grep(pattern = genotype, x = metadata$genotype), ]

import_plates(mat_nm = mat_nm, metadata = metadata)
```

```
## will read TetAAV_P01
```

```
## will read TetAAV_P07
```

```
## will read TetAAV_P08
```

```
## will read TetAAV_P10
```

```
## will read TetAAV_P11
```

```
## will read TetAAV_P12
```

```
## will read TetAAV_P16
```

```
## will read TetAAV_P18
```

```
## will read TetAAV_P13_2
```

```
## will read TetAAV_P19_2
```

```
## [1] TRUE
```


```r
#### next step: gating of populations
```


