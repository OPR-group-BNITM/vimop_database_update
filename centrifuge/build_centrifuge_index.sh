#!/bin/bash
#SBATCH -c 16 # number of cores
#SBATCH --mem 128000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate centrifuge

index_build_dir="${VIMUPDATE_CENTRIFUGE}/centrifuge"

if [[ -d $index_build_dir ]]
then
    rm -rf $index_build_dir
fi

mkdir -p $index_build_dir
cd $index_build_dir

refseq_taxons=(
    "archaea"
    "bacteria"
    "homo_sapiens"
    "mus_musculus"
)
# TODO: add aedes aegypti?

fname_fasta_merged=merged.fasta
fname_seqid_to_taxid_merged=seqid2taxid.tsv

: > $fname_fasta_merged
: > $fname_seqid_to_taxid_merged

for taxon in "${refseq_taxons[@]}"
do
    fname_fasta="${VIMUPDATE_GENOMES}/refseq/${taxon}.fasta.gz"
    gunzip -c $fname_fasta >> $fname_fasta_merged
    fname_seqid_to_taxid="${VIMUPDATE_GENOMES}/refseq/${taxon}.seqid2taxid.tsv"
    cat $fname_seqid_to_taxid >> $fname_seqid_to_taxid_merged
done

cat "${VIMUPDATE_GENOMES}/merged_ncbi_rvdb/viruses.seqid2taxid.tsv" >> $fname_seqid_to_taxid_merged
cat "${VIMUPDATE_DB}/virus/ALL.fasta" >> $fname_fasta_merged

cp "${VIMUPDATE_GENOMES}/ncbi_virus/virus_taxids.txt" "virus_taxids.txt"
ktUpdateTaxonomy.sh .

centrifuge-build \
    -p 8 \
    --conversion-table $fname_seqid_to_taxid_merged \
    --taxonomy-tree "${VIMUPDATE_TAXONKIT}/nodes.dmp" \
    --name-table "${VIMUPDATE_TAXONKIT}/names.dmp" \
    $fname_fasta_merged \
    all

{
    echo "version: ${VIMUPDATE_CENTRIFUGEDB_VERSION}"
    echo "description: \"${VIMUPDATE_CENTRIFUGEDB_DESCRIPTION}\""
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

cd ..
mv centrifuge "${VIMUPDATE_DB}/centrifuge"
