FROM bioconductor/bioconductor_docker:devel

WORKDIR /home/rstudio

COPY --chown=rstudio:rstudio . /home/rstudio/

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install(ask=FALSE)"

RUN R -e "install.packages('remotes')"

RUN R -e "remotes::install_github(c('rformassspectrometry/PTMods', 'rformassspectrometry/PSMatch@validatePSM'), build_vignettes = TRUE)"

RUN Rscript -e "options(repos = BiocManager::repositories()); devtools::install('.', dependencies=TRUE, build_vignettes=TRUE)"
