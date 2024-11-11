#!/usr/bin/env bash

set -x
outdir="../data/lassa/filtered/clustered"
mkdir -p "$outdir"

for segment in S L
do
    fname_in="../data/lassa/filtered/filtered_segment_${segment}.fasta"
    for thresh in all 0.99
    do
        fname_clust="$outdir/${segment}_${thresh}.fasta"
        fname_msa="$outdir/${segment}_${thresh}.msa.fasta"
        if [ "$thresh" != "all" ]; then
            cd-hit-est -i $fname_in -o $fname_clust -c $thresh
        else
            cp $fname_in $fname_clust
        fi
        muscle -super5 $fname_clust -output $fname_msa
    done
done
