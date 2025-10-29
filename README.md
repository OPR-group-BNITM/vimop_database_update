# ViMOP database setup

This repository contains scripts and documentation on how to setup automatically setup the OPR ViMOP database.
The datatabase has three independent compentents

- the virus reference genomes
- the centifuge index
- the host genomes

This pipeline will build a new data base based on the most recent genomes available, which will be downloaded automatically.

## Downloading the data

First, sequence data need to be downloaded,
All scripts to download the necessary sequences are found in the directory `download_genomes`.

## Building the data base components

The three component are build with the scripts found in `virus`, `centrifuge` and `host` respectively.
Building the centrifuge index requires to build the virus reference data set first.
