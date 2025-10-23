#!/bin/bash
#SBATCH -c 8 # number of cores
#SBATCH --mem 128000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"

workflow_dir="${VIMOP_DB_UPDATE_SRC}/virus"
cd $workflow_dir

config_path="${workflow_dir}/configs/db${VIMOP_DB_UPDATE_VIRUSDB_VERSION}/${VIMOP_DB_UPDATE_VIRUSDB_VERSION}.config"

nextflow main.nf -c $config_path -profile conda
