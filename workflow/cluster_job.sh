#!/bin/bash
#SBATCH -c 16 # number of cores
#SBATCH --mem 128000 # memory pool for all cores
#SBATCH -o /data2/le-petersen/nils/slurmout/slurm.db_build.%N.%j.out # STDOUT
#SBATCH -e /data2/le-petersen/nils/slurmout/slurm.db_build.%N.%j.err # STDERR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nils.petersen@bnitm.de
#SBATCH --job-name=db_build

source /opt/conda/etc/profile.d/conda.sh
export PATH="/opt/conda/bin:$PATH"
conda activate aligner

cd /data/home/nils.petersen/dev/VirusDatasetCuration/workflow
nextflow main.nf -c /data/home/nils.petersen/dev/VirusDatasetCuration/workflow/configs/db2.0.config
# nextflow main.nf -c /data/home/nils.petersen/dev/VirusDatasetCuration/workflow/testset/test.config
