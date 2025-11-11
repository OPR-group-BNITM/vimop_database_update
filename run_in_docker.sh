#!/usr/bin/env bash

version=$(cat docker_env/general/version.txt)

docker run --rm -it \
    -v "$(pwd):/work" \
    -w /work \
    oprgroup/vimop_db_update_general:${version} \
    bash -lc './run_all.sh'  
