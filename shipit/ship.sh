#!/usr/bin/env bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL

out_dir="${VIMUPDATE_DB}/split"
mkdir -p $out_dir

cd "${VIMUPDATE_SRC}/shipit"

for db_name in virus centrifuge contaminants 
do
    data_dir="${VIMUPDATE_DB}/${db_name}"
    ./archive_and_split.sh $data_dir $db_name $out_dir
done
