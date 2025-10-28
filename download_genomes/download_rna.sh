#!/usr/bin/env bash
set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate datasets

OUT_DIR="${VIMUPDATE_GENOMES}/refseq"
TMP_DIR="${OUT_DIR}/out_rna"

datasets download genome taxon human --reference --assembly-source refseq --include rna
unzip ncbi_dataset.zip -d "${TMP_DIR}"

{
    # This is usually one file, but we do not know the name of the latest assembly
    for fn in "${TMP_DIR}"/ncbi_dataset/data/GCF_*/rna.fna
    do
        cat "$fn"
    done
} | gzip > "${OUT_DIR}/human_rna.fasta.gz"

rm -r "${OUT_DIR}/out_rna"
