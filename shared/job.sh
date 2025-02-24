#!/bin/bash
#SBATCH -c 16 # number of cores
#SBATCH --mem 128000 # memory pool for all cores
#SBATCH -o /data2/le-petersen/nils/slurmout/slurm.%N.%j.out # STDOUT
#SBATCH -e /data2/le-petersen/nils/slurmout/slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nils.petersen@bnitm.de

set -x

threads=16

usage() {
    echo "Usage: $0 {ALIGN|FILTER|MSA|ALL} casename"
    exit 1
}

if [ $# -lt 2 ]; then
    usage
fi
option=$1
casename=$2

if [[ "$option" != "ALIGN" && "$option" != "FILTER" && "$option" != "MSA" && "$option" != "ALL" ]]; then
    echo "Invalid option: $option"
    usage
fi

source /opt/conda/etc/profile.d/conda.sh
export PATH="/opt/conda/bin:$PATH"

VIRUS_DATASET_CURATION_SRC=/data/home/nils.petersen/dev/VirusDatasetCuration
fname_config=$VIRUS_DATASET_CURATION_SRC/datasets/$casename/align_orient.nextflow.config

if [[ "$option" == "ALIGN" || "$option" == "ALL" ]]; then
    conda activate aligner
    nextflow $VIRUS_DATASET_CURATION_SRC/shared/workflow/align_to_refs.nf -c $fname_config -with-conda
fi

fname_sequences=$(grep -E "fasta_seqs=" "$fname_config" | awk -F= '{print $2}' | tr -d "'" | tr -d '"')
datadir=$(grep -E "out_dir=" "$fname_config" | awk -F= '{print $2}' | tr -d "'" | tr -d '"')

dir_notebook="$datadir/notebook"
dir_filt="$datadir/filtered"
dir_clust="$datadir/clustered"

if [[ "$option" == "FILTER" || "$option" == "ALL" ]]; then

    conda activate jupyter
    echo "Running script in conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

    mkdir -p $dir_notebook

    fname_stats="$datadir/all/collected_stats.tsv"

    papermill $VIRUS_DATASET_CURATION_SRC/shared/filter.ipynb $dir_notebook/filtering.ipynb \
        -p fname_stats "$fname_stats" \
        -p fname_sequences "$fname_sequences" \
        -p outdir "$dir_filt"

    jupyter nbconvert --to html  $dir_notebook/filter.ipynb
fi

if [[ "$option" == "MSA" || "$option" == "ALL" ]]; then

    conda activate msa
    echo "Running script in conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"

    thresh=0.98
    $VIRUS_DATASET_CURATION_SRC/shared/cluster_and_msa.sh $dir_clust $thresh $threads $dir_filt/*.fasta
fi
