#!/usr/bin/env bash

# https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image

# https://support.rstudio.com/hc/en-us/community/posts/202827628-Using-HTTPS-instead-of-HTTP
# https://www.digitalocean.com/community/tutorials/how-to-configure-nginx-with-ssl-as-a-reverse-proxy-for-jenkins

PASS="8_juggleD_albiNo_12_eleVens_cRush"
DOCKER_MNTPOINT="/home/rstudio/crne_cna1crz1_rnaseq"
HOST_BASE="$HOME/crne_cna1crz1_rnaseq"
HOST_SCRATCH="/mnt/ibiem_scratch/Members/josh/crne_cna1crz1_rnaseq"
WORKSPACE="$HOST_SCRATCH/workspace"
RAW_DATA="$HOST_SCRATCH/raw_data"

if [ "$1" == "shell" ]; then
    DOCKER_COMMAND="/bin/bash"
    CONTAINER_NAME="rstudio_shell"
    DOCKER_ARGS="--rm --interactive --tty --user rstudio"
    echo "------------------------------"
    echo "In docker run the following:"
    echo "cd $DOCKER_MNTPOINT"
    echo "make --dry-run fastqs -j4 -f analyze_culture_sequences.mk --warn-undefined-variables NUMTHREADS=4 RAW_DATA_DIR=$DOCKER_MNTPOINT/raw_data"
    echo "------------------------------"
elif [ "$1" == "rstudio" ]; then
    DOCKER_COMMAND=""
    CONTAINER_NAME="rstudio_web"
    DOCKER_ARGS="--detach --publish 8787:8787 -e PASSWORD=$PASS"
else
   echo "Must supply command line argument! Should be run as one of the following commands:"
   echo "'$0 shell'"
   echo "'$0 rstudio'"
   exit 1
fi



if [ ! -d "$WORKSPACE" ]; then
    echo "NOT MOUNTED: $WORKSPACE"
    echo "REFUSING TO START DOCKER"
    exit 1
fi

if [ ! -d "$RAW_DATA" ]; then
    echo "NOT MOUNTED: $RAW_DATA"
    echo "REFUSING TO START DOCKER"
    exit 1
fi

docker run $DOCKER_ARGS \
       --name $CONTAINER_NAME \
       -v $HOST_BASE:$DOCKER_MNTPOINT \
       -v $WORKSPACE:$DOCKER_MNTPOINT/workspace \
       -v $RAW_DATA:$DOCKER_MNTPOINT/raw_data \
       granek/rnaseq:v0rc0 \
       $DOCKER_COMMAND
