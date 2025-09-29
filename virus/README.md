# Virus data set curation for metagenomic sequencing analysis

This repository contains scripts to create data sets for different virus species.
The data sets are used as references in our metagenomics pipeline (currently named vimop).

## Overview

The workflow to create a new data base version requires several steps.

- Download sequence data from NCBI virus and RVDB and merge them
- set up config files to create the data base
  - get taxonomy information
- run the data set curation pipeline
- review the outcome (using a jupyter notebook)

## Dowload virus genomes

Read the following to learn how to download and merge data from RVDB and NCBI virus.

### Download RVDB

Go to the [RVDB-website](https://rvdb.dbi.udel.edu/).
Download the latest Unclustered DB data set and save it in `data/input_genomes/U-RVDBvX.fasta.gz` where you replace X with the version number (e.g. 30.0).

### Download NCBI virus genomes without Covid

Open this [NCBI virus link](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Viruses,%20taxid:10239&VirusLineage_ss=taxid:NOT%202697049). It's all viruses with covid excluded. Click the following

- Download All Results
- Nucleotide
- Donwload All Records
- Build Custom
  - Accession
  - GenBankTitle
  - Accession
  - Species
- Download

Save this as `data/input_genomes/ncbi_nocovid_YYYYMMDD` and replace date in the format year, month, data (e.g. 20251029 for the 29th octobre 2025).

### Merge the NCBI and RVDB-Covid data

Run the script
```
python scripts/merge_ncbi_rvdb_covid.py \
--email user@host.de \
--rvdb data/input_genomes/C-RVDBv30.0.fasta.gz \
--ncbi data/input_genomes/ncbi_nocovid_20250922.fasta \
--output data/input_genomes/ncbi_nocovid_20250922_crvdbv30.0_covid.fasta
```
Replace the email with yours and set the versions in the fasta filenames according to what you have downloaded. 
It will merge the NCBI genomes with the covid genomes from RVDB and the Covid reference genome.
The Covid reference genome is downloaded from NCBI using Entrez, that is why you have to add your email adress.  

## Create the configs

### Define the curated data set

Create a configuration file to define curated data sets and filters.
Save this in `configs/dbX.Y/groups_and_refs.yaml` (replace X and Y with your version).
Typically you want to copy this from a previous version and modify it if needed.

Example:
```
curated:
  COVID:
    name: Severe acute respiratory syndrome coronavirus 2
    taxa:
    - taxid: 2697049
    segments:
      Unsegmented:
        refs:
        - NC_045512.2
        seqs:
        - NC_045512.2
  LASV:
    name: Mammarenavirus lassense
    taxa:
    - taxid: 3052310
    segments:
      S:
        refs:
        - NC_004296.1
        seqs:
        - NC_004296.1
        - KM822128.1
        - GU481068.1
        - OL774861.1
      L:
        refs:
        - NC_004297.1
        seqs:
        - NC_004297.1
        - KM822127.1
        - GU481069.1
        - OL774860.1
filters:
  ARENA:
    name: Arenaviridae
    taxa:
    - taxid: 11617
```

The first section `curated` holds the curated species.
`refs` is the place to put one reference that defines the orientation in which all genomes will be deposited.
`seqs` are more reference genomes.
Note, that all genomes in the input data set will be considered for this data set.
The sequences here are used to compare all sequences to in order to filter out genomes in the data set, that are assigned falsely to this taxon.

### Add organism names

Run 
```
python scripts/fetch_organisms_for_taxids.py \
  --config configs/db2.2/groups_and_refs.yaml \
  --out configs/db2.2/groups_refs_and_organisms.yaml
```
to add organism names.

## Run a data set curation

There is a test run in `testset`. You can run it by typing `nextflow main.nf -c testset/test.config`. The results will be written to `data/output_databases/test`. The reference virus genome data base is then found in `data/output_databases/test/db`. It contains files with genomes for the respective filters, a yaml config file and the blast data base.

By default this pipeline runs with docker, but you can also choose other profiles using `-profile conda` or `-profile noenv`.
The latter requires the dependencies to be installed in the sourrounding system.  
This can be useful, if you want to run this in a slurm container. 

## Review the outcome

To have a look into the statistics of the curated data sets use the jupyter notebook workflow/utils/evaluate.ipynb and change the paths and data set name to the one you are interested in.

## Virus with taxid

The following viruses and families are currently in our data sets.

### Viruses

| Virus                              | Abbreviation | TaxId     |
| ---------------------------------- | ------------ | --------- |
| Mammarenavirus lassense            | LASV         | 3052310   |
| Mammarenavirus choriomeningitidis  | LCMV         | 305230    |
| Mammarenavirus juninense           | JUNV         | 2169991   |
| Orthoebolavirus                    | EBOV         | 3044781   |
| Orthomarburgvirus                  | MARV         | 3044783   |
| Orthoflavivirus denguei            | DENV         | 3052464   |
| Orthoflavivirus zikaense           | ZIKA         | 3048459   |
| Emesvirus zinderi                  | MS2          | 329852    |
| Yellow fever virus                 | YFV          | 3046277   |
| West nile virus                    | WNV          | 3048448   |
| Orthonairovirus hazaraense         | HAZV         | 3052519   |


### Families

| Virus family            | Abbreviation | TaxId     |
| ----------------------- | ------------ | --------- |
| Arenaviridae            | ARENA        | 11617     |
| Filoviridae             | FILO         | 11266     |
| Hantaviridae            | HANTA        | 1980413   |
| Nairoviridae            | NAIRO        | 1980415   |

## Container setup for slurm

To run this inside a container (e.g. in slurm system) you need to install the following inside the container:
- conda (best to use an image with conda and mamba already installed)
- java (apt-get install -y openjdk-11-jdk-headless)
- nextflow:
    - wget -qO- https://get.nextflow.io | bash && \
      chmod +x nextflow && \
      mv nextflow /usr/local/bin/
- the conda environment in docker_images/general/env.yaml
