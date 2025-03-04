# Virus data set curation for metagenomic sequencing analysis

This repository contains scripts to create data sets for different virus species. The data sets are used as references in our metagenomics pipeline (currently named nanoflow).

Each data set should have an own readme which documents the data set and it's creation.

## Virus with taxid

### Viruses

| Virus                              | Abbreviation | TaxId     |
| ---------------------------------- | ------------ | --------- |
| Mammarenavirus lassense            | LASV         | 3052310   |
| Mammarenavirus choriomeningitidis  | LCMV         | 305230    |
| Mammarenavirus juninense           | JUNV         | 2169991   |
| Orthoebolavirus                    | EBOV         | 3044781   |
| Orthomarburgvirus                  | EBOV         | 3044783   |
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

## Gaia installation

- created a container using init_container.sh. The container is called nextflow.
- started the container with run_instance.sh.
  - installed java (apt-get install -y openjdk-11-jdk-headless)
  - installed nextflow 
    - wget -qO- https://get.nextflow.io | bash && \
      chmod +x nextflow && \
      mv nextflow /usr/local/bin/
  - created conda environments
    - conda create -n aligner
      - python=3.11
      - seqtk=1.4
      - minimap2=2.1.1
      - samtools=1.21
      - pysam=0.22.1
      - biopython=1.85
      - pandas=2.2.3
    - conda create -n jupyter
      - python=3.12
      - jupyter=1.1.1
      - papermill=2.6.0
      - pandas=2.2.2.
      - matplotlib=3.10.0
      - biopython=1.85
    - conda create -n msa
      - cd-hit=4.8.1
      - muscle=5.3
