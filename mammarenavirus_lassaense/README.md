# Lassa virus data set creation

There are two steps. The first is to align and orient the reference sequences. The second is to align all sequences in our raw data set


## Raw data set

This is typically all sequences downloaded from NCBI.

## Orient the reference sequences

In the directory refset.

## Orient and filter the whole data set

TODO

## Input data requirements

- fasta header need to be formatted like this: ID|description|family|species

## TODO

- orient the sequences and add segment information to output file
- filter the sequences in notebook and create final fasta (add forward and reverse tags to fasta headers + orient the sequences)

- CD-hit in the main workflow?

- MSA?
    - add msa for the referece sequences?
    - add visualization of the msa of the references to the jupyter notebook?
    - add msa for all final sequences
    - msa integritiy test
        - do pairwise assignments change with the msa creation?

- add test for integrity of the input data (e.g. formatting of the fasta header)
