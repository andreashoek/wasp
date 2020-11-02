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
    xtail \
    wget


# Install required packages for Shiny App
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyFiles', 'shinyjs', 'scales', 'rjson', 'ggplot2', 'DT', 'readr', 'plotly'))"

# Expose Port for connection  
EXPOSE 3838

# Copy Shiny app & required files
COPY app/ ./app

CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
