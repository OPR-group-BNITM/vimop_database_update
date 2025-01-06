#!/usr/bin/env bash


fname_stats="data/lassa/all/collected_stats.tsv"
fname_sequences="/Volumes/DataCurate/CurationDatasets/Lassa/lassa_ncbi_20241024.fasta"
outdir="data/lassa/filtered"


papermill ../shared/filter.ipynb lassa_filter.ipynb \
    -p fname_stats "$fname_stats" \
    -p fname_sequences "$fname_sequences" \
    -p outdir "$outdir"

#jupyter nbconvert --to html lassa_filter.ipynb
