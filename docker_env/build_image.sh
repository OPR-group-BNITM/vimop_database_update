#!/usr/bin/env bash

# build the only the image. Use this if you do not want to use slurm, 
# but want to us run_in_docker.sh to run everything locally instead.

set -euo pipefail

prefix=vimop_db_update
images=(general)
if [[ $# -gt 0 ]]
then
    images=("$@")
fi

for img in "${images[@]}"
do
    image="${prefix}_${img}"
    version=$(cat ${img}/version.txt)
    
    echo $img:$version
    docker build --no-cache -t oprgroup/${img}:$version $img/
    echo ""
done
