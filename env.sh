#!/usr/bin/env bash

export VIMOP_DB_UPDATE_VIRUS_DB_INPUT="/data2/le-petersen/shared_data/db_different_versions/update_data/virus_input_genomes"
export VIMOP_DB_UPDATE_OUTPUT_DIR="/data2/le-petersen/shared_data/db_different_versions/update_data"
export VIMOP_DB_UPDATE_SRC="/data/home/nils.petersen/dev/VimopDatabaseUpdate"
export VIMOP_DB_UPDATE_SRC_SLURM_OUTPUT_DIR="/data2/le-petersen/shared_data/db_different_versions/update_data/slurm"
export VIMOP_DB_UPDATE_USER_MAIL="nils.petersen@bnitm.de"

export VIMOP_DB_UPDATE_VIRUSDB_VERSION="2.3"  # used to build centrifuge index from virus DB file
export VIMOP_DB_UPDATE_VIRUS_ALL_FILE="${VIMOP_DB_UPDATE_OUTPUT_DIR}/virus_db/db${VIMOP_DB_UPDATE_VIRUSDB_VERSION}/virus/ALL.fasta"

export VIMOP_DB_UPDATE_OUTPUT_REFSEQ_DATA="${VIMOP_DB_UPDATE_OUTPUT_DIR}/centrifuge/refseq_data"
export VIMOP_DB_UPDATE_VIRUS_TAXONDATA_DIR="${VIMOP_DB_UPDATE_OUTPUT_DIR}/centrifuge/taxon_data"

export VIMOP_DB_UPDATE_CENTRIFUGEDB_VERSION="1.1"
export VIMOP_DB_UPDATE_CENTRIFUGEDB_DESCRIPTION="Viruses sequences from ViMOP virus DB version ${VIMOP_DB_UPDATE_VIRUSDB_VERSION}, RefSeq files for Archea, Bacteria and Human downloaded in Octobre 2025"
