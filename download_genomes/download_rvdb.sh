#!/usr/bin/env bash
set -euo pipefail

outdir="$VIMUPDATE_GENOMES/rvdb"
mkdir -p $outdir
cd $outdir
wget $VIMUPDATE_RVDB_URL
