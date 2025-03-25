#!/usr/bin/env bash

# Setup
# after installation (conda install taxonkit) go to /Users/nils.petersen/.taxonkit
# and download taxonomy data (wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) and uncompress  

taxonkit list --ids 10239 | sed 's/.* //' > virus_taxids.txt
