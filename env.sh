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
export VIMUPDATE_DB="${VIMUPDATE_BASEDIR}/db"

# URLs
export VIMUPDATE_TAXDUMP_URL="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
export VIMUPDATE_RVDB_URL="https://rvdb.dbi.udel.edu/download/C-RVDBvCurrent.fasta.gz"

# Output DB versions and descriptions
export VIMUPDATE_VIRUSDB_VERSION="2.4"
export VIMUPDATE_VIRUSDB_DESCRIPTION="NCBI virus genomes combined with covid sequences from RVDB"

export VIMUPDATE_CENTRIFUGEDB_VERSION="1.1"
export VIMUPDATE_CENTRIFUGEDB_DESCRIPTION="Centrifuge index with virus Genbank sequences and RefSeq files for Archea, Bacteria and Human"

export VIMUPDATE_CONTAMINANTSDB_VERSION="1.2"
export VIMUPDATE_CONTAMINANTSDB_DESCRIPTION="Human, mouse, mastomys, aedes aegypti and contaminant filter set"

# Taxonomy
export VIMUPDATE_COVID_REFSEQ="NC_045512.2"
export VIMUPDATE_VIRUS_TAXID="10239"
export VIMUPDATE_COVID_TAXID="2697049"




# export VIMOP_DB_UPDATE_VIRUS_DB_INPUT="/data2/le-petersen/shared_data/db_different_versions/update_data/virus_input_genomes"

# export VIMOP_DB_UPDATE_VIRUSDB_VERSION="2.3"  # used to build centrifuge index from virus DB file
# export VIMOP_DB_UPDATE_VIRUS_ALL_FILE="${VIMOP_DB_UPDATE_OUTPUT_DIR}/virus_db/db${VIMOP_DB_UPDATE_VIRUSDB_VERSION}/virus/ALL.fasta"

# export VIMOP_DB_UPDATE_OUTPUT_REFSEQ_DATA="${VIMOP_DB_UPDATE_OUTPUT_DIR}/centrifuge/refseq_data"
# export VIMOP_DB_UPDATE_VIRUS_TAXONDATA_DIR="${VIMOP_DB_UPDATE_OUTPUT_DIR}/centrifuge/taxon_data"
