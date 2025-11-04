#!/usr/bin/env bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate datasets

out_dir="${VIMUPDATE_DB}/split"
mkdir -p "$out_dir/files"

cd "${VIMUPDATE_SRC}/shipit"

for db_name in virus centrifuge contaminants 
do
    data_dir="${VIMUPDATE_DB}/${db_name}"
    ./archive_and_split.sh $data_dir $db_name "${out_dir}/files"
done

python merge_config.py \
    --virus "${out_dir}/files/virus.v${VIMUPDATE_VIRUSDB_VERSION}.files.yaml" \
    --contaminants "${out_dir}/files/contaminants.v${VIMUPDATE_CONTAMINANTSDB_VERSION}.files.yaml" \
    --centrifuge "${out_dir}/files/centrifuge.v${VIMUPDATE_CENTRIFUGEDB_VERSION}.files.yaml" \
    --version "${VIMUPDB_VERSION}" \
    --description "${VIMUPDB_DESCRIPTION}" \
    --output "${out_dir}/vimop_db.${VIMUPDB_VERSION}.yaml"
