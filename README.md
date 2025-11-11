# ViMOP database setup

This repository contains scripts and documentation on how to automatically create an updated database for our pipeline [ViMOP](https://github.com/OPR-group-BNITM/vimop).
The scripts

- download all required sequences from NCBI GenBank and RefSeq (`download_genomes`)
- filter the virus genomes to build the virus reference dataset (`virus`)
- build the centrifuge index (`centrifuge`)
- collect the host genomes (`host`)
- create zipped packages to make the database downloadable (`shipit`)

The scripts are made to run on our slurm system.

## Build docker container

Before running them, build a docker container with the script `build_containers.sh` in `docker_env`.

## Set versions and paths

Before building a database, open env.sh and set all the versions and paths required.
Additionally, set the 
```bash
export VIMUPDATE_FILESHARE_PREFIX="USER@SERVER:/path/to/files" 
```
in your own .bashrc to upload the database to the distribution server.

## Build the database

Send the jobs to the slurm system with
```bash
./slurmrun_all.sh
```

## Build the data base locally (not in a slurm system)

Instead of building a container, only build the image (`build_image.sh` in `docker_env`).
Then run `./run_in_docker.sh`.

## Support

This repository is built for our team to create the database and it is shared for transparency and reproducability.
However, if you want to use these scripts on your system and need help, please feel free to contact us here!
