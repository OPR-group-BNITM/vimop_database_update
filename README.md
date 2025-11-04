# ViMOP database setup

This repository contains scripts and documentation on how to automatically create an updated OPR ViMOP database.
The scripts

- download all required sequences from NCBI GenBank and RefSeq (`download_genomes`)
- filter the virus genomes to build the virus reference dataset (`virus`)
- build the centrifuge index (`centrifuge`)
- collect the host genomes (`host`)
- create zipped packages to make the database downloadable (`shipit`)

The scripts are made to run on our slurm system.
Before running them, build a docker container with the script `docker_env/build_containers.sh`.
Then send the job to the hpc using `./slurmrun_build_vimop_database.sh`.
