FROM rocker/r-ver:4.0.0

RUN apt-get update && apt-get install -y \
    sudo \
    gdebi-core \
    libhdf5-dev \
    libssl-dev \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libxml2-dev \
    libglpk-dev \
    xtail \
    wget


# Download and install shiny server
RUN apt-get -y update && apt-get -y install \
    python3-pip

RUN pip3 install \
    umap-learn

RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    . /etc/environment && \
    R -e "install.packages(c('shiny', 'rmarkdown', 'BiocManager', 'shinydashboard', 'ggplot2', 'dplyr', 'deldir', 'Rcpp', 'Seurat', 'DT', 'umap'))" && \
    R -e "BiocManager::install(c('SingleCellExperiment','scater'))"

EXPOSE 3838

COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh
COPY /app/app.R /srv/shiny-server/app.R

CMD ["/usr/bin/shiny-server.sh"]
