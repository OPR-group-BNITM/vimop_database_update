#!/usr/bin/env bash


fname_sequences="/Volumes/DataCurate/CurationDatasets/Dengue/dengue_ncbi_20240105.fasta"
fname_stats="data/all/collected_stats.tsv"
outdir="data/filtered"


papermill ../../shared/filter.ipynb filtering.ipynb \
    -p fname_stats "$fname_stats" \
    -p fname_sequences "$fname_sequences" \
    -p outdir "$outdir"

#jupyter nbconvert --to html lassa_filter.ipynb
