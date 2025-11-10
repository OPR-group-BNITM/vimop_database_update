# ViMOP database setup

This repository contains scripts and documentation on how to automatically create an updated database for our pipeline [ViMOP](https://github.com/OPR-group-BNITM/vimop).
The scripts

- download all required sequences from NCBI GenBank and RefSeq (`downloads`)
- filter the virus genomes to build the virus reference dataset (`virus`)
- build the centrifuge index (`centrifuge`)
- collect the host genomes (`host`)
- create zipped packages to make the database downloadable (`shipit`)

The scripts are made to run on our slurm system.

## Build docker container

Before running them, build a docker container with the script `docker_env/build_containers.sh`.

## Set versions and paths

Before building a database, open env.sh and set all the versions and paths required.
Additionally, set the 
```bash
set export VIMUPDATE_FILESHARE_PREFIX="USER@SERVER:/path/to/files" 
```
in your own .bashrc to upload the database to the distribution server.

## Build the database

Send the jobs to the slurm system with
```bash
./slurmrun_all.sh
```

## Support
