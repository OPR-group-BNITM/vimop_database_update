#!/usr/bin/env bash

# build the images and a container with the same name in $HOME to use with slurm.

set -euo pipefail

set -x

images=(dbupdate)

for image in "${images[@]}"
do
    version=$(cat $image/version.txt)

    container="$image"
    container_dir="${HOME}/containers/${container}"
    rootfs_dir="${container_dir}/rootfs"

    docker rm -f $container 2>/dev/null || true
    rm -rf $container_dir

    docker build --no-cache -t $image:$version $image/

    mkdir -p "${HOME}/run"
    mkdir -p $rootfs_dir

    docker export $(docker create --name ${container} "${image}:${version}") | tar -C ${rootfs_dir} -xvf -

    /opt/gaia-tools/bin/helper/generate-config -hostname "${container}" -path ${container_dir}
done
