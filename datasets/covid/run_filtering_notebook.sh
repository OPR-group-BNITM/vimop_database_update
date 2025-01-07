#!/usr/bin/env bash

set -x

fname_config=align_orient.nextflow.config
fname_sequences=$(grep -E "fasta_seqs=" "$fname_config" | awk -F= '{print $2}' | tr -d "'" | tr -d '"')
fname_stats="data/all/collected_stats.tsv"
outdir="data/filtered"


papermill ../../shared/filter.ipynb filtering.ipynb \
    -p fname_stats "$fname_stats" \
    -p fname_sequences "$fname_sequences" \
    -p outdir "$outdir"

#jupyter nbconvert --to html lassa_filter.ipynb
