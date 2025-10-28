#!/usr/bin/env bash

# basic settings
export VIMUPDATE_MAIL="nils.petersen@bnitm.de"
export VIMUPDATE_SRC="/data/home/nils.petersen/dev/VimopDatabaseUpdate"
export VIMUPDATE_BASEDIR="/data/home/nils.petersen/data/vimop_db_update/"

# output sub-directories
export VIMUPDATE_SLURM="${VIMUPDATE_BASEDIR}/slurm"
export VIMUPDATE_GENOMES="${VIMUPDATE_BASEDIR}/genomes"
export VIMUPDATE_TAXONKIT="${VIMUPDATE_BASEDIR}/taxonkit"
export VIMUPDATE_CENTRIFUGE="${VIMUPDATE_BASEDIR}/centrifuge"

# URLs
export VIMUPDATE_TAXDUMP_URL="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
export VIMUPDATE_RVDB_URL="https://rvdb.dbi.udel.edu/download/C-RVDBvCurrent.fasta.gz"

# Output DB versions and descriptions
export VIMUPDATE_CENTRIFUGEDB_VERSION="1.1"
export VIMUPDATE_CENTRIFUGEDB_DESCRIPTION="Viruses sequences from ViMOP virus DB version ${VIMUPDATE_CENTRIFUGEDB_VERSION}, RefSeq files for Archea, Bacteria and Human downloaded in Octobre 2025"

# Taxonomy
export VIMUPDATE_COVID_REFSEQ="NC_045512.2"
export VIMUPDATE_VIRUS_TAXID="10239"
export VIMUPDATE_COVID_TAXID="2697049"




# export VIMOP_DB_UPDATE_VIRUS_DB_INPUT="/data2/le-petersen/shared_data/db_different_versions/update_data/virus_input_genomes"

# export VIMOP_DB_UPDATE_VIRUSDB_VERSION="2.3"  # used to build centrifuge index from virus DB file
# export VIMOP_DB_UPDATE_VIRUS_ALL_FILE="${VIMOP_DB_UPDATE_OUTPUT_DIR}/virus_db/db${VIMOP_DB_UPDATE_VIRUSDB_VERSION}/virus/ALL.fasta"

# export VIMOP_DB_UPDATE_OUTPUT_REFSEQ_DATA="${VIMOP_DB_UPDATE_OUTPUT_DIR}/centrifuge/refseq_data"
# export VIMOP_DB_UPDATE_VIRUS_TAXONDATA_DIR="${VIMOP_DB_UPDATE_OUTPUT_DIR}/centrifuge/taxon_data"

# export VIMOP_DB_UPDATE_CENTRIFUGEDB_VERSION=
# export VIMOP_DB_UPDATE_CENTRIFUGEDB_DESCRIPTION=