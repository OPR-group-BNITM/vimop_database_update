#!/usr/bin/env bash

casename=MARV
container_name=nextflow

sbatch --container=$HOME/containers/${container_name} \
    --job-name=$casename \
    ../../shared/job.sh ALL $casename
