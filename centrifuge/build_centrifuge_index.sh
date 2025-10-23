#!/bin/bash
#SBATCH -c 8 # number of cores
#SBATCH --mem 100000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate centrifuge

index_build_dir="${VIMOP_DB_UPDATE_OUTPUT_DIR}/centrifuge/centrifuge"

if [[ -d $index_build_dir ]]
then
    rm -rf $index_build_dir
fi

mkdir -p $index_build_dir
cd $index_build_dir

refseq_taxons=(
    "archaea"
    "bacteria"
    "human"
    "mouse"
)

fname_fasta_merged=merged.fasta
fname_seqid_to_taxid_merged=seqid2taxid.tsv

for taxon in "${refseq_taxons[@]}"
do
    fname_fasta="${VIMOP_DB_UPDATE_OUTPUT_REFSEQ_DATA}/${taxon}.fasta.gz"
    gunzip -c $fname_fasta >> $fname_fasta_merged
    fname_seqid_to_taxid="${VIMOP_DB_UPDATE_OUTPUT_REFSEQ_DATA}/${taxon}.seqid2taxid.tsv"
    cat $fname_seqid_to_taxid >> $fname_seqid_to_taxid_merged
done

fname_seqid_to_taxid_virus="${VIMOP_DB_UPDATE_VIRUS_TAXONDATA_DIR}/virus.seqid2taxid.tsv"
fname_fasta_virus="${VIMOP_DB_UPDATE_VIRUS_ALL_FILE}"

cat $fname_seqid_to_taxid_virus >> $fname_seqid_to_taxid_merged
cat $fname_fasta_virus >> $fname_fasta_merged


cp "${VIMOP_DB_UPDATE_VIRUS_TAXONDATA_DIR}/virus_taxids.txt" "virus_taxids.txt"
ktUpdateTaxonomy.sh .

centrifuge-build \
    -p 8 \
    --conversion-table $fname_seqid_to_taxid_merged \
    --taxonomy-tree "${VIMOP_DB_UPDATE_VIRUS_TAXONDATA_DIR}/taxonkit/nodes.dmp" \
    --name-table "${VIMOP_DB_UPDATE_VIRUS_TAXONDATA_DIR}/taxonkit/names.dmp" \
    $fname_fasta_merged \
    all

{
    echo "version: ${VIMOP_DB_UPDATE_CENTRIFUGEDB_VERSION}"
    echo "description: \"${VIMOP_DB_UPDATE_CENTRIFUGEDB_DESCRIPTION}\""
    echo "virus_taxid_file: virus_taxids.txt"
    echo "index_name: all"
    echo "files:"
    for f in all.*.cf
    do
        echo "- \$f"
    done
} > centrifuge.yaml

mkdir -p taxonomy
mv taxonomy.tab taxonomy/taxonomy.tab

rm images.dmp
rm $fname_seqid_to_taxid_merged
rm $fname_fasta_merged
