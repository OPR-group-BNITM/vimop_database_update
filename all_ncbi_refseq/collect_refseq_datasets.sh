#!/usr/bin/env bash

set -x

# Input
fname_all_ncbi="/Volumes/DataCurate/CurationDatasets/NCBI_all/ncbi_all_NoCovid_20241026.fasta"

# Output
outdir="/Volumes/DataCurate/CurationDatasets/NCBI_all/grouped_by_species_with_refseqs"
fname_ids_organisms="$outdir/fasta_header_content.tsv"

mkdir -p $outdir

echo "extracting fasta headers"
grep "^>" $fname_all_ncbi | awk -F'[>|]' '{print $2 "\t" $3 "\t" $4 "\t" $5}' > "$fname_ids_organisms"

python make_refseq_datasets.py "$fname_ids_organisms" "$fname_all_ncbi" "$outdir/species"

# echo "Extracting fasta sequences"

# for species_dir in $outdir/species/*
# do
#     seqtk subseq "$fname_all_ncbi" "$species_dir/ids.txt" > "$species_dir/sequences.fasta"
# done

# echo "done"

# TODO in a separate script:
# - for each organism
#   - make list with reference IDs -> create the yaml/nextflow config
#   - copy a shell script to run nextflow with the config
