# ViMOP DB centrifuge index creation

This repository contains scripts to build the vimop database centrifuge index.

Before you need to
- set all paths and configs in `env.sh`. Be sure to set version and description of the centrifuge part.
- download the data (see `download_genomes`)
- build the virus data base (see `virus`)

Then run `./slurmrun_build_centrifuge_index.sh`
