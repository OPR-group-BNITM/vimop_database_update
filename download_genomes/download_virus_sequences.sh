#!/usr/bin/env bash
set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate datasets

source "${VIMUPDATE_SRC}/download_genomes/utils.sh"

OUT_DIR="${VIMUPDATE_GENOMES}/ncbi_virus"

VIRUSES="${OUT_DIR}/viruses.nocovid.tsv"
SEQS="${OUT_DIR}/ncbi_virusviruses.fasta.gz"

TMP_DIR="${OUT_DIR}/tmp_virus_genomes_download"

mkdir -p $TMP_DIR

# get the sequence IDs 
SEQIDS="${TMP_DIR}/sequences.ids.txt"
awk 'NR>1{print $1}' "$VIRUSES" | LC_ALL=C sort -u > "$SEQIDS"

echo "Sequences to download: $(wc -l < "$SEQIDS")"

# split across multiple files
BATCH_DIR="$TMP_DIR/batches"
mkdir -p "$BATCH_DIR"

split -l 20000 -d -a 4 "$SEQIDS" "${BATCH_DIR}/seqids_"

{
    i=0
    for f in "$BATCH_DIR"/seqids_*
    do
        ((++i))
        zipfile="${TMP_DIR}/viruses_part${i}.zip"
        echo "Downloading batch $i ..."  >&2
        run_with_retries 8 datasets download virus genome accession \
            --inputfile "$f" \
            --include genome \
            --filename "$zipfile"
        unzip -qq -p "$zipfile" "ncbi_dataset/data/genomic.fna" 
        rm "$zipfile"
    done
} | gzip > "$SEQS"

rm -r "${TMP_DIR}"

echo "DONE. Sequences written to ${SEQS}"