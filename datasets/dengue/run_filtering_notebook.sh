#!/usr/bin/env bash


fname_sequences="path_to_fasta_file_with_all_sequence_to_filter"
fname_stats="data/all/collected_stats.tsv"
outdir="data/filtered"


papermill ../../shared/filter.ipynb filtering.ipynb \
    -p fname_stats "$fname_stats" \
    -p fname_sequences "$fname_sequences" \
    -p outdir "$outdir"

#jupyter nbconvert --to html lassa_filter.ipynb
