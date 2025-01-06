#!/usr/bin/env bash

datadir="data/lassa/filtered"
outdir="$datadir/clustered"
thresh=0.98

../shared/cluster_and_msa.sh $outdir $thresh $datadir/S.fasta $datadir/L.fasta
