FROM continuumio/miniconda3:4.12.0

LABEL maintainer="Jordan Brown"
LABEL maintainer.email="jordbrown.10@gmail.com"

RUN mkdir -p /usr/share/man/man1

WORKDIR /project

RUN apt-get --allow-releaseinfo-change update -y && \
    apt-get install -y unzip \
    build-essential \
    libz-dev \ 
    r-base-core libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev pandoc \
    r-base-dev \
    default-jre \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    r-cran-xml \
    pandoc \
    cutadapt \
    bedtools


RUN Rscript -e """\
install.packages(c(\"languageserver\", \"systemfonts\", \"devtools\", \"ggplot2\", \"dplyr\", \"reshape2\", \"ggsci\", \"BiocManager\", \"openxlsx\", \"rmarkdown\")) \
"""

RUN Rscript -e """\
BiocManager::install(c(\"Biostrings\"), update = TRUE, ask = FALSE) \
"""

RUN pip install radian

RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    { conda config --remove repodata_fns current_repodata.json || true ; } && \
    conda config --prepend repodata_fns repodata.json && \
    conda config --set auto_update_conda False 

#Setup File System
RUN mkdir data
VOLUME ["/data"]


CMD ["/bin/bash"]
