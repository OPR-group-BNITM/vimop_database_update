#!/bin/bash
#SBATCH -c 8 # number of cores
#SBATCH --mem 128000 # memory pool for all cores
#SBATCH -o /data2/le-petersen/shared_data/slurmout/slurm.db_build.%N.%j.out # STDOUT
#SBATCH -e /data2/le-petersen/shared_data/slurmout/slurm.db_build.%N.%j.err # STDERR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nils.petersen@bnitm.de
#SBATCH --job-name=db_build

workflow_dir=/data2/le-petersen/shared_data/db_different_versions/VimopDatabaseUpdate/virus
cd $workflow_dir
config_path=configs/db2.2/db2.2.config
nextflow main.nf -c $config_path -profile conda