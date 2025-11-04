#!/usr/bin/env bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 8000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate datasets

# Force EDirect away from curl
export EDIRECT_PREFER_CURL=0
export EDIRECT_PREFER_WGET=1

# Silence locale warnings
export LC_ALL=C.UTF-8
export LANG=C.UTF-8
export LANGUAGE=en_US:en

export NCBI_EMAIL=$VIMUPDATE_MAIL

fname_genbank_tsv="${VIMUPDATE_SRC}/host/reagent/sequences.tsv"
fname_refseq_tsv="${VIMUPDATE_SRC}/host/reagent/assemblies.tsv"

OUT_DIR="${VIMUPDATE_GENOMES}/contaminant_list"
fname_out="${OUT_DIR}/contaminants.fasta.gz"
mkdir -p "${OUT_DIR}"

{
    awk '{print $1}' "$fname_genbank_tsv" \
    | epost -db nuccore \
    | efetch -format fasta

    awk '{print $1}' "$fname_refseq_tsv" \
    | epost -db assembly \
    | elink -target nuccore -name assembly_nuccore_refseq \
    | efetch -format fasta
} | gzip > "$fname_out"
