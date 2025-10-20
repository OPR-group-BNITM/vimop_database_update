#!/usr/bin/env python


import os
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genomes', required=True, help='')
    parser.add_argument('--outdir', required=True, help='')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    fname_seqid_to_tax = os.path.join(args.outdir, 'virus_seqid_to_tax.tsv')
    fams_set = set()
    species_set = set()

    with open(fname_seqid_to_tax, 'w') as f_out:
        for record in SeqIO.parse(args.genomes, 'fasta'):
            seqid = record.id
            descr_split = record.description.split('|')
            family = descr_split[-4]
            fams_set.add(family)
            species = descr_split[-3]
            species_set.add(species)
            f_out.write(f'{seqid}\tViruses\t{family}\t{species}\n')

    with open(os.path.join(args.outdir, 'virus_families.txt'), 'w') as f_fam:
        f_fam.write('\n'.join(fams_set))

    with open(os.path.join(args.outdir, 'virus_species.txt'), 'w') as f_species:
        f_species.write('\n'.join(species_set))


if __name__ == '__main__':
    main()
