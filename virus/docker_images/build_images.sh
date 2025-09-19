#!/usr/bin/env bash

images=(general)

for img in "${images[@]}";
do
    echo $img
    docker build --no-cache -t db_${img}:0.0.1 $img/
    echo ""
done
