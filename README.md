
# The intrinsic and extrinsic effects of Tet proteins during gastrulation

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/500408241.svg)](https://zenodo.org/badge/latestdoi/500408241)
<!-- badges: end -->

This is accompanying code that generates the figures of our paper on TET function during gastrulation. The code is splitted into jupyter notebooks that can be found in the analysis folder.

## Running the notebooks

Prior to any analysis, after cloning the repository, please download first the necessary data by running (in the root directory of the cloned repository):


```bash
R -e "source('scripts/download_data.R'); download_full_data()"
```

This will download all the necessary data including processed metacell objects necessary for generating the figures. You can also download a minimal version ( `download_minimal_scrna_data()` ) in which case you would need to rerun all the notebooks prior to _generate-figures_ (see section *Notebook order* below.)

## Necessary R packages

You can install the necessary R packages by running:

```r
remotes::install_github("tanaylab/tet-gastrulation", subdir = "scripts/tet.gastru")
```

The initialization script (`scripts/init.R`) loads automatically the necessary R packages to run the notebooks. The analysis was done using R 4.0.5 and the following packages:

- gridExtra 2.3
- qvalue 2.22.0
- ggrepel 0.9.1
- ggpubr 0.4.0   
- ggplot2 3.3.5
- zoo 1.8-9
- qlcMatrix 0.9.7
- sparsesvd 0.2  
- slam 0.1-49
- tidyr 1.2.0
- dplyr 1.0.9
- tgutil 0.1.13  
- tgstat 2.3.17
- metacell 0.3.7
- Matrix 1.3-4

## Notebook order 

You should run the notebooks _import_plates_ and _gating_ first and _generate-figures_ last.

1. import_plates
2. gating
3. TKO-chimera-embryos-processing
4. TKO-tetraploid-complemented-embryos-processing
5. Control-tetraploid-complemented-embryos-processing
6. DKO-chimera-embryos-processing
7. AAV_Cre-TKO-embryos-processing
8. Methylation
9. generate-figures

## Docker

We also provide a docker image which contains all the needed dependencies, to use it run:

Download the analysis files: 

```bash
R -e "source('scripts/download_data.R'); download_full_data()"
```

Change permissions for the analysis files:

```bash
chmod a+wx data/
chmod a+rw -R data/
chmod a+rx -R db/
chmod a+rx scrna_db/
chmod a+rw -R scrna_db/
chmod a+w figs/
chmod a+w output/
```

Run the container:

```bash
docker run \
    -v $(pwd)/db:/workdir/db \
    -v $(pwd)/data:/workdir/data \
    -v $(pwd)/scrna_db:/workdir/scrna_db \
    -v $(pwd)/output:/workdir/output \
    -v $(pwd)/figs:/workdir/figs  \
    -ti \
    -p 8888:8888 \
    tanaylab/tet-gastrulation
```

Connect to the jupyter server running at port 8888.



