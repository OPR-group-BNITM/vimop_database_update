#!/usr/bin/env bash

tasks=(
    "rs_setup download_genomes"
    "virus_taxinfo download_genomes"
    "virus_sequences download_genomes"
    "rvdb download_genomes"
    "merge_ncbi_rvdb_covid download_genomes"
    "rna download_genomes"
    "contaminants download_genomes"
    "virus_db virus"
    "contaminant_db host"
    "build_centrifuge_index centrifuge"
    "ship shipit"
)

for task in "${tasks[@]}"
do
    set -- $task   # split
    task_name=$1
    directory=$2

    echo "Running ${directory}/${task_name}.sh"

    cd "$directory"
    ./${task_name}.sh
    cd ..
done
