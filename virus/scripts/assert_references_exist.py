"""Check if the sequence-IDs that are used as reference sequences for the curated datasets exist in our virus genomes.
"""
import argparse
import gzip
import yaml
import sys

from Bio import SeqIO


def gzip_open(fname, mode='r'):
    if fname.endswith('.gz'):
        return gzip.open(fname, f'{mode}t')
    return open(fname, mode)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--taxgroups', required=True, help='Groups with the organism names')
    parser.add_argument('--fasta', required=True, help='Fasta file to check')
    args = parser.parse_args()

    with open(args.taxgroups) as f_in:
        groups = yaml.safe_load(f_in)

    ids_to_search = {
        seqid
        for curated_set in groups.get('curated', {}).values()
        for ref_and_seqs in curated_set.get('segments', {}).values()
        for seqlist in ref_and_seqs.values()
        for seqid in seqlist
    }

    with gzip_open(args.fasta) as f_in:
        ids_found = {record.id for record in SeqIO.parse(f_in, 'fasta')} 

    ids_to_add = sorted(ids_to_search - ids_found)
    if ids_to_add:
        print('Missing sequences:')
        print('\n'.join(ids_to_add))
        sys.exit(1)


if __name__ == '__main__':
    main()
