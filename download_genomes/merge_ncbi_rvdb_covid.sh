#!/usr/bin/env bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate datasets

mkdir -p "${VIMUPDATE_GENOMES}/merged_ncbi_rvdb_virus"

python ${VIMUPDATE_SRC}/download_genomes/merge_ncbi_rvdb_covid.py \
    --email "$VIMUPDATE_MAIL" \
    --virus-taxid "$VIMUPDATE_VIRUS_TAXID" \
    --covid-taxid "$VIMUPDATE_COVID_TAXID" \
    --covid-refseq-id "$VIMUPDATE_COVID_REFSEQ" \
    --rvdb "${VIMUPDATE_GENOMES}/rvdb/C-RVDBvCurrent.fasta.gz" \
    --ncbi "${VIMUPDATE_GENOMES}/ncbi_virus/ncbi_virusviruses.fasta.gz" \
    --taxinfo "${VIMUPDATE_GENOMES}/ncbi_virus/viruses.nocovid.tsv" \
    --fasta-out "${VIMUPDATE_GENOMES}/merged_ncbi_rvdb_virus/virus.fasta" \
    --seqid-to-taxid-out "${VIMUPDATE_GENOMES}/merged_ncbi_rvdb_virus/virus.seqid2taxid.map"
