# Use official R base image with Ubuntu
FROM rocker/r-ver:4.3.3

# Set maintainer
LABEL maintainer="your-email@domain.com"
LABEL description="ecDNA Repeat Element Enrichment Analysis using regioneR"
LABEL version="1.0"

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV R_REPOS=https://cran.rstudio.com/

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libgit2-dev \
    bedtools \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages('BiocManager', repos='$R_REPOS')" && \
    R -e "BiocManager::install(c('regioneR', 'GenomicRanges', 'rtracklayer', 'BSgenome.Hsapiens.UCSC.hg38', 'BSgenome.Hsapiens.UCSC.hg19'))" && \
    R -e "install.packages(c('ggplot2', 'gridExtra', 'dplyr', 'optparse', 'readr'), repos='$R_REPOS')"

# Create working directory
WORKDIR /app

# Copy R scripts
COPY ecdna_analysis.R /app/
COPY run_analysis.R /app/

# Make scripts executable
RUN chmod +x /app/*.R

# Create output directory
RUN mkdir -p /app/output

# Set default command
ENTRYPOINT ["Rscript", "/app/run_analysis.R"]

# Default help message if no arguments provided
CMD ["--help"]
