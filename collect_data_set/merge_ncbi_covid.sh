#!/usr/bin/env bash

fname_covid="/data2/le-petersen/nils/CurationDatasets/Covid/CovidRVDB/covid_rvdb_v29.fasta"
fname_ncbi="/data2/le-petersen/nils/CurationDatasets/NCBI_all/ncbi_all_NoCovid_20241026.fasta"
fname_out="/data2/le-petersen/nils/base_datasets/merged_datasets/ncbi_all_NoCovid_20241026_rvdb29covid.fasta"

# merge ncbi virus without covid with covid sequences
cat $fname_covid $fname_ncbi > $fname_out
