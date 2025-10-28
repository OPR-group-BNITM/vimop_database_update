#!/usr/bin/env bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail

outdir="$VIMUPDATE_GENOMES/rvdb"
mkdir -p $outdir
cd $outdir
wget $VIMUPDATE_RVDB_URL
