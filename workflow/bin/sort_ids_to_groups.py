#!/usr/bin/env python3
import os
import argparse
import yaml
from Bio import SeqIO


def read_remove_set(fname):
    with open(fname) as f_in:
        return {
            line.strip()
            for line in f_in
            if not line.strip().startswith('#')
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--taxgroups', required=True, help='Groups with the organism names')
    parser.add_argument('--sequences', required=True, help='Fasta file with sequences')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--category', required=True, help='Category like "curated" or "filters"')
    parser.add_argument('--remove_set', required=True, help='Text file with sequence IDs to remove.')
    parser.add_argument(
        '--organism_label_position',
        default=1,
        help='Position of organism name in segments separated by "|" in fasta header. Counted from the end.',
        type=int,
    )
    args = parser.parse_args()

    fname_groups = args.taxgroups
    fname_seqs = args.sequences
    dir_out = args.outdir
    category = args.category  # curated or filter

    with open(fname_groups) as f_in:
        groups = yaml.safe_load(f_in)

    remove_set = read_remove_set(args.remove_set)

    # build a dictionary organism name -> groupkey
    # (e.g. mammarenavirus lassaense -> LASV)
    organism_to_group = {}
    for groupkey, group_dict in groups[category].items():
        for taxon in group_dict['taxa']:
            for organism in taxon['organisms']:
                org_key = organism.lower()
                has_non_matching_entry = (
                    category == 'curated'
                    and org_key in organism_to_group
                    and organism_to_group[org_key] != groupkey
                )
                if has_non_matching_entry:
                    raise RuntimeError(
                        f'Organism {org_key} found in two groups '
                        f'{organism_to_group[org_key]} and {groupkey}'
                    )
                organism_to_group.setdefault(organism, []).append(groupkey)

    # Extract the sequence IDs and write them to files
    os.makedirs(dir_out, exist_ok=True)
    outfiles = {
        groupkey: open(os.path.join(dir_out, f'{groupkey}.txt'), 'w')
        for groupkey in list(groups[category]) + ['NOGROUP']
    }

    for seq in SeqIO.parse(fname_seqs, 'fasta'):
        if seq.id in remove_set:
            continue
        organism = seq.description.rsplit('|')[-args.organism_label_position].strip()
        for groupkey in organism_to_group.get(organism, ['NOGROUP']):
            outfiles[groupkey].write(seq.id + '\n')

    for f_out in outfiles.values():
        f_out.close()


if __name__ == '__main__':
    main()
