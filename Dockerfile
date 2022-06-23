FROM jupyter/datascience-notebook:r-4.1.3

# Install rpm dependencies
USER root
RUN apt-get update && apt-get install -y cmake git-core libcurl4-openssl-dev libgit2-dev libicu-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev libgtk2.0-dev libhiredis-dev libcairo2-dev libxt-dev xvfb xauth xfonts-base vim && rm -rf /var/lib/apt/lists/*
USER ${NB_UID}

# Install reticulate
RUN R -e 'install.packages("reticulate", repos = "http://cran.us.r-project.org")'

ENV JUPYTER_ENABLE_LAB=yes

USER root

COPY ./scripts /workdir/scripts
COPY ./config /workdir/config
COPY ./analysis_notebooks /workdir/analysis_notebooks

RUN R -e 'remotes::install_local("/workdir/scripts/tet.gastru", Ncpus = 40)'

RUN chown -R ${NB_UID} /workdir
RUN chgrp -R ${NB_GID} /workdir

USER ${NB_UID}

RUN touch /workdir/.here

WORKDIR /workdir

