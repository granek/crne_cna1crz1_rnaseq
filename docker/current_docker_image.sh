DOCKER_IMAGE_TAG="DIv0rc0"
DOCKER_IMAGE_NAME="granek/rnaseq"
DOCKER_IMAGE="${DOCKER_IMAGE_NAME}:${DOCKER_IMAGE_TAG}"

DOCKER_MNTPOINT="/home/rstudio/crne_cna1crz1_rnaseq"
WORKSPACE_MNTPOINT="$DOCKER_MNTPOINT/workspace"
DATA_MNTPOINT="$DOCKER_MNTPOINT/raw_data"

HOST_BASE="$HOME/crne_cna1crz1_rnaseq"
HOST_SCRATCH="/mnt/ibiem_scratch/Members/josh/crne_cna1crz1_rnaseq"
WORKSPACE="$HOST_SCRATCH/workspace"
RAW_DATA="$HOST_SCRATCH/raw_data"


TIME_ZONE="-e TZ=America/New_York -v /etc/timezone:/etc/timezone"

PORT_NUMBER=8784
