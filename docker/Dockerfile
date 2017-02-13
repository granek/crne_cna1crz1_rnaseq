# Distributed under the terms of the Modified BSD License.

FROM rocker/rstudio:3.3.2

MAINTAINER Josh Granek

USER root

RUN apt-get update &&  \
    apt-get install -y --no-install-recommends \
    fastqc default-jre \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN echo "deb http://ftp.debian.org/debian jessie-backports main" >  /etc/apt/sources.list.d/backports.list && \
    apt-get update &&  \
    apt-get -t jessie-backports install -y --no-install-recommends \
    libxml2-dev \
    less \
    make \
    git \
    samtools \
    tophat \
    libgsl0-dev \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# download and compile ea-utils, since version from jessie-backports depends on R,
# which fails
# libgsl0-dev is for compiling ea-utils
RUN git clone https://github.com/ExpressionAnalysis/ea-utils.git && \
    cd ea-utils/clipper && \
    make && \
    make install && \
    rm -rf ../../ea-utils/

RUN Rscript -e "install.packages(pkgs = c('optparse','RColorBrewer','gplots','formatR','rprojroot','rmarkdown','dplyr','tidyr'), \
    repos='https://cran.revolutionanalytics.com/', \
    dependencies=TRUE)" && \
    Rscript -e "source('https://bioconductor.org/biocLite.R'); \
    biocLite(pkgs=c('DESeq2'))"

##------------------------------------------------------------
# USER $RSTUDIO_USER

# ##------------------------------------------------------------
# USER root

CMD ["/init"]
