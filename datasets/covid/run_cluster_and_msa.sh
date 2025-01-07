#!/usr/bin/env bash

datadir="data/filtered"
outdir="data/clustered"
thresh=0.98

# Add the fasta files for all segments
../../shared/cluster_and_msa.sh $outdir $thresh $datadir/*.fasta
