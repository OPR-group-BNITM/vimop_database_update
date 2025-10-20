#!/bin/bash
#SBATCH -c 4 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate env

n_threads=4
script_dir="$VIMOP_DB_UPDATE_SRC/centrifuge"
fname_virus_genomes="${VIMOP_DB_UPDATE_VIRUS_ALL_FILE}"
taxdump_url="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
outdir="${VIMOP_DB_UPDATE_VIRUS_TAXONDATA_DIR}"
taxonkit_datadir="${outdir}/taxonkit"
virus_taxid="10239"

workdir=$(pwd)

rm -rf "$taxonkit_datadir"
mkdir -p "$taxonkit_datadir"

cd "$taxonkit_datadir"
wget "${taxdump_url}"
tar -xzf taxdump.tar.gz
cd "$workdir"

echo "$virus_taxid" | taxonkit list --data-dir "$taxonkit_datadir" > "${outdir}/virus_taxids_one.txt"
sed 's/.* //' "${outdir}/virus_taxids_one.txt" > "${outdir}/virus_taxids.txt"

python ${script_dir}/get_family_species.py \
    --genomes $fname_virus_genomes \
    --outdir ${outdir}

taxonkit name2taxid --data-dir $taxonkit_datadir --threads ${n_threads} < "${outdir}/virus_families.txt" > "${outdir}/families.taxids"
taxonkit name2taxid --data-dir $taxonkit_datadir --threads ${n_threads} < "${outdir}/virus_species.txt" > "${outdir}/species.taxids"

echo -e "Viruses\t${virus_taxid}" > "${outdir}/kingdoms.taxids"

python ${script_dir}/seqids_to_taxids.py \
    --taxon-table  "${outdir}/virus_seqid_to_tax.tsv" \
    --kingdom-taxids "${outdir}/kingdoms.taxids" \
    --family-taxids "${outdir}/families.taxids" \
    --species-taxids "${outdir}/species.taxids" \
    --seqid-to-taxid "${outdir}/virus.seqid2taxid.tsv"
