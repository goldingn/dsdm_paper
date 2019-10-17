FROM rocker/binder:3.6.0
LABEL maintainer='Nick Golding'
USER root
COPY . ${HOME}
RUN chown -R ${NB_USER} ${HOME}
USER ${NB_USER}



RUN wget https://github.com/goldingn/dsdm_paper/raw/master/DESCRIPTION && R -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-10-17/')); devtools::install_deps()"

RUN rm DESCRIPTION.1; exit 0
