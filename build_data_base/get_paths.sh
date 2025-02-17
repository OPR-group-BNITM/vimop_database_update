#!/usr/bin/env bash

echo "curated:"

for f in /data2/le-petersen/nils/curated_virus_reference_datasets/*/clustered/*.clust_0.98.fasta
do
    echo $f | awk -F"/" '{print "  "$6 ":"}'
    echo "    organisms:"
    grep ">" $f | awk -F"|" '{print $4}' | sort | uniq | awk '{print "    - "$0}'
        
    segment=$(basename "${f%.clust_0.98.fasta}")
    echo "    segments:"
    echo "    - $segment"
    echo "      file: '$f'"
    msa="${f%.fasta}.muscle.fasta"
    if [ -e $msa ]
    then
        echo "      msa: '$msa'"
    else
        echo "      msa: ''"
    fi
done