#!/usr/bin/env bash

# build the images and a container with the same name in $HOME to use with slurm.

set -euo pipefail

prefix=vimop_db_update
images=(nextflow ncbi_datasets)
if [[ $# -gt 0 ]]
then
  images=("$@")
fi

for img in "${images[@]}"
do
    image="${prefix}_${img}"
    version=$(cat ${img}/version.txt)

    container="$image"
    container_dir="${HOME}/containers/${container}"
    rootfs_dir="${container_dir}/rootfs"

    docker rm -f $container 2>/dev/null || true
    rm -rf $container_dir

    docker build --no-cache -t $image:$version $img/

    mkdir -p "${HOME}/run"
    mkdir -p $rootfs_dir

    docker export $(docker create --name ${container} "${image}:${version}") | tar -C ${rootfs_dir} -xvf -

    /opt/gaia-tools/bin/helper/generate-config -hostname "${container}" -path ${container_dir}
done
