FROM bioconductor/bioconductor_docker:devel

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install(ask=FALSE)"

RUN R -e "install.packages(c('magick', 'remotes'))"

RUN R -e "BiocManager::install(c('BiocStyle', 'MsDataHub', 'QFeatures', 'mzR', 'mzID'))"

RUN R -e "remotes::install_github(c('rformassspectrometry/PTMods', 'rformassspectrometry/PSMatch@validatePSM'), build_vignettes = TRUE, dependencies = TRUE)"

RUN Rscript -e "options(repos = BiocManager::repositories()); devtools::install('.', dependencies=TRUE, build_vignettes=TRUE, force=TRUE)"
