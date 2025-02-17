#!/bin/bash
#SBATCH -c 16 # number of cores
#SBATCH --mem 64000 # memory pool for all cores
#SBATCH -o /data2/le-petersen/nils/slurmout/slurm.%N.%j.out # STDOUT
#SBATCH -e /data2/le-petersen/nils/slurmout/slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nils.petersen@bnitm.de

source /opt/conda/etc/profile.d/conda.sh
export PATH="/opt/conda/bin:$PATH"

conda activate yaml

python /data/home/nils.petersen/dev/VirusDatasetCuration/build_data_base/build_virus_db.py \
    /data/home/nils.petersen/dev/VirusDatasetCuration/build_data_base/virus_db.yaml \
    test_version_0.0.1
