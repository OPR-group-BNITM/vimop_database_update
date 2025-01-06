#!/usr/bin/env bash


# Check if at least three arguments are provided (output directory, threshold, and at least one input file)
if [[ $# -lt 3 ]]; then
    echo "Usage: $0 <output_directory> <threshold> <input_file1> [input_file2 ...]"
    exit 1
fi

# Assign the first argument to the output directory variable
outdir="$1"

# Assign the second argument to the threshold variable
thresh="$2"

# Validate that the threshold is a positive number
if [ "$thresh" != "all" ]; then
    if ! [[ "$thresh" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
        echo "Error: Threshold must be a positive number."
        exit 1
    fi
fi

# Create the output directory if it doesn't exist
mkdir -p "$outdir"

# Shift the first two arguments to get the input files
shift 2

# Iterate over the input files
for fname_in in "$@"; do
    echo "Processing file: $fname_in with threshold: $thresh"

    base_name=$(basename "$fname_in" | sed -E 's/\.[^.]+$//')
    fname_clust="$outdir/${base_name}.clust_${thresh}.fasta"
    fname_msa="$outdir/${base_name}.clust_${thresh}.msa.fasta"

    if [ "$thresh" != "all" ]; then
        cd-hit-est -i $fname_in -o $fname_clust -c $thresh
    else
        cp $fname_in $fname_clust
    fi
    muscle -super5 $fname_clust -output $fname_msa
done

echo "All files have been processed and the output is written to $outdir."
