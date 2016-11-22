# Distributed under the terms of the Modified BSD License.

FROM rocker/rstudio

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
    ea-utils \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(pkgs = c('optparse','RColorBrewer','gplots'), \
    repos='https://cran.revolutionanalytics.com/', \
    dependencies=TRUE)" && \
    Rscript -e "source('https://bioconductor.org/biocLite.R'); \
    biocLite(pkgs=c('DESeq2'))"

##------------------------------------------------------------
USER $RSTUDIO_USER

##------------------------------------------------------------
USER root

CMD ["/init"]
