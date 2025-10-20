#!/bin/bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL


set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate env

outdir=${VIMOP_DB_UPDATE_OUTPUT_REFSEQ_DATA}

#rm -rf $outdir
mkdir -p $outdir
cd $outdir

taxa=(
    #"archaea 2157"
    "bacteria 2"
    "human 9606"  # homo sapiens
    "mouse 10090"  # mus musculus
)

for tax in "${taxa[@]}"
do
    set -- $tax   # split
    taxon=$1
    taxid=$2

    datasets summary genome taxon "$taxid" \
        --reference --assembly-source refseq --as-json-lines \
        | dataformat tsv genome --fields accession,organism-tax-id,organism-name \
        > "${taxon}_tax.tsv"

    awk '{print $1}' "${taxon}_tax.tsv" | tail -n+2 > "${taxon}_reference_acc.txt"

    datasets download genome accession \
        --inputfile "${taxon}_reference_acc.txt" \
        --assembly-source refseq \
        --dehydrated \
        --include genome \
        --filename "${taxon}_refseq.zip"

    unzip -q ${taxon}_refseq.zip -d "${taxon}_pkg"
    datasets rehydrate --directory "${taxon}_pkg"

    python ${VIMOP_DB_UPDATE_SRC}/centrifuge/merge_refseq_seqs.py \
        --sequence-directory "${taxon}_pkg/ncbi_dataset/data" \
        --assembly-to-tax "${taxon}_tax.tsv" \
        --kingdom-taxid "$taxid" \
        --prefix "${PWD}/${taxon}"

    rm -r "${taxon}_pkg/"
    rm "${taxon}_refseq.zip"

done
