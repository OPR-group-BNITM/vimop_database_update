#!/usr/bin/env bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 8000 # memory pool for all cores
#SBATCH --mail-type=ALL

outdir="${VIMUPDATE_DB}/contaminants"

fname_human_dna="${VIMUPDATE_GENOMES}/refseq/homo_sapiens.fasta.gz"
fname_human_rna="${VIMUPDATE_GENOMES}/refseq/human_rna.fasta.gz"
fname_aedes_aegypti="${VIMUPDATE_GENOMES}/refseq/aedes_aegypti.fasta.gz"
fname_mus_musculus="${VIMUPDATE_GENOMES}/refseq/mus_musculus.fasta.gz"
fname_mastomys_natalensis="${VIMUPDATE_GENOMES}/refseq/mastomys_natalensis.fasta.gz"
fname_contaminants="${VIMUPDATE_GENOMES}/contaminant_list/contaminants.fasta.gz"

fname_human_dna_out=human_genome.fasta.gz
fname_human_rna_out=human_transcriptome.fasta.gz
fname_aedes_aegypti_out=aedis_aegyptis_genome.fasta.gz
fname_mus_musculus_out=mus_musculus_genome.fasta.gz
fname_mastomys_natalensis_out=mastomys_natalensis_genome.fasta.gz
fname_contaminants_out=contaminants.fasta.gz

mkdir -p $outdir
cd $outdir

cp "$fname_human_dna" "$fname_human_dna_out"
cp "$fname_human_rna" "$fname_human_rna_out"
cp "$fname_aedes_aegypti" "$fname_aedes_aegypti_out" 
cp "$fname_mus_musculus" "$fname_mus_musculus_out"
cp "$fname_mastomys_natalensis" "$fname_mastomys_natalensis_out"
cp "$fname_contaminants" "$fname_contaminants_out"

{
    echo "filters:"
    echo "  human_dna: \"${fname_human_dna_out}\""
    echo "  human_rna: \"${fname_human_rna_out}\""
    echo "  mouse: \"${fname_mus_musculus_out}\""
    echo "  mastomys: \"${fname_mastomys_natalensis_out}\""
    echo "  aedes_aegypti: \"${fname_aedes_aegypti_out}\""
    echo "  reagent: \"${fname_contaminants_out}\""
    echo "version: \"${VIMUPDATE_CONTAMINANTSDB_VERSION}\""
    echo "description: \"${VIMUPDATE_CONTAMINANTSDB_DESCRIPTION}\""
} > contaminants.yaml
