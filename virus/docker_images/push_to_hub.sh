#!/usr/bin/env bash

images=(general)
docker_user=oprgroup

for img in "${images[@]}";
do
    version=$(cat $img/version.txt)
    echo $img:$version
    docker push ${docker_user}/vimopdb_${img}:$version
    echo ""
done
