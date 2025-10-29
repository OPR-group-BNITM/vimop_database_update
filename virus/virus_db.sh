#!/bin/bash
#SBATCH -c 8 # number of cores
#SBATCH --mem 128000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"

workflow_dir="${VIMUPDATE_SRC}/virus"
cd $workflow_dir

fname_taxa_config="$VIMUPDATE_VIRUS/groups_refs_and_organisms.yaml"

mkdir -p "$VIMUPDATE_VIRUS"


conda activate datasets

python scripts/fetch_organisms_for_taxids.py \
    --email "$VIMUPDATE_MAIL" \
    --config "configs/groups_and_refs.yaml" \
    --out $fname_taxa_config

conda deactivate


nextflow main.nf \
    --taxa_config $fname_taxa_config \
    --remove_set configs/remove.txt \
    --fasta_sequences "${VIMUPDATE_GENOMES}/merged_ncbi_rvdb_virus/virus.fasta" \
    --out_dir "$VIMUPDATE_VIRUS" \
    --output_version "$VIMUPDATE_VIRUSDB_VERSION" \
    --output_description "$VIMUPDATE_VIRUSDB_DESCRIPTION" \
    -profile conda


mv "${VIMUPDATE_VIRUS}/virus" "${VIMUPDATE_DB}/virus"
