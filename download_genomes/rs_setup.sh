#!/bin/bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL


set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate datasets

outdir="${VIMUPDATE_GENOMES}/refseq"

mkdir -p $outdir
cd $outdir

all_levels="complete,chromosome,contig,scaffold"

taxa=(
    "archaea 2157 refseq ${all_levels}"
    "homo_sapiens 9606 refseq ${all_levels}"
    "mus_musculus 10090 refseq ${all_levels}"
    "mastomys_natalensis 10112 genbank ${all_levels}"
    "aedes_aegypti 7159 refseq ${all_levels}"
    "bacteria 2 refseq complete"
)

for tax in "${taxa[@]}"
do
    set -- $tax   # split
    taxon=$1
    taxid=$2
    database=$3
    assembly_level=$4

    if [[ -e "${taxon}_tax.tsv" ]]
    then
        echo "${taxon}_tax.tsv exists, exiting" >&2
        exit 1
    fi

    datasets summary genome taxon "$taxid" \
        --assembly-source "$database" \
        --assembly-level $assembly_level \
        --reference \
        --as-json-lines \
        | dataformat tsv genome --fields accession,organism-tax-id,organism-name,assminfo-level,assminfo-refseq-category \
        > "${taxon}_tax.tsv"

    awk '{print $1}' "${taxon}_tax.tsv" | tail -n+2 > "${taxon}_reference_acc.txt"

    datasets download genome accession \
        --inputfile "${taxon}_reference_acc.txt" \
        --assembly-source "$database" \
        --dehydrated \
        --include genome \
        --filename "${taxon}_refseq.zip"

    unzip -q ${taxon}_refseq.zip -d "${taxon}_pkg"
    datasets rehydrate --directory "${taxon}_pkg"

    python "${VIMUPDATE_SRC}/download_genomes/merge_refseq_seqs.py" \
        --sequence-directory "${taxon}_pkg/ncbi_dataset/data" \
        --assembly-to-tax "${taxon}_tax.tsv" \
        --kingdom-taxid "$taxid" \
        --prefix "${PWD}/${taxon}"

    rm -r "${taxon}_pkg/"
    rm "${taxon}_refseq.zip"

done
