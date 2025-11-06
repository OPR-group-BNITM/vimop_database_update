#!/usr/bin/env bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate datasets

out_dir="${VIMUPDATE_DB}"
mkdir -p "$out_dir"

cd "${VIMUPDATE_SRC}/shipit"

for db_name in virus centrifuge contaminants 
do
    data_dir="${VIMUPDATE_DB}/${db_name}"
    ./archive_and_split.sh $data_dir $db_name "${out_dir}"
done

python merge_config.py \
    --virus "${out_dir}/split/files/virus.v${VIMUPDATE_VIRUSDB_VERSION}.files.yaml" \
    --contaminants "${out_dir}/split/files/contaminants.v${VIMUPDATE_CONTAMINANTSDB_VERSION}.files.yaml" \
    --centrifuge "${out_dir}/split/files/centrifuge.v${VIMUPDATE_CENTRIFUGEDB_VERSION}.files.yaml" \
    --version "${VIMUPDATE_DB_VERSION}" \
    --description "${VIMUPDATE_DB_DESCRIPTION}" \
    --output "${out_dir}/split/vimop_db.${VIMUPDATE_DB_VERSION}.yaml"
