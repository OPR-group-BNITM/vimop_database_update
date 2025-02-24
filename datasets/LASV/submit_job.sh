#!/usr/bin/env bash

casename=LASV
container_name=nextflow

sbatch --container=$HOME/containers/${container_name} \
    --job-name=$casename \
    ../../shared/job.sh MSA $casename
