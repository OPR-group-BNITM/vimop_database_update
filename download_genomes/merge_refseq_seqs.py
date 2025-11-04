#!/usr/bin/env python3


import argparse
import gzip
from pathlib import Path
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sequence-directory', required=True, help='')
    parser.add_argument('--assembly-to-tax', required=True, help='')
    parser.add_argument('--kingdom-taxid', required=True, help='')
    parser.add_argument('--prefix', required=True, help='')
    args = parser.parse_args()

    seq_path = Path(args.sequence_directory)

    assembly_to_tax = {}
    with open(args.assembly_to_tax) as f_in:

        # check the header validity
        head_assembly_id, head_tax_id, head_organism, *_ = next(f_in).strip().split('\t')
        assert head_assembly_id == 'Assembly Accession'
        assert head_tax_id == 'Organism Taxonomic ID' 
        assert head_organism == 'Organism Name'

        for line in f_in:
            assembly_id, taxid, *_ = line.split('\t')  # not using organism-name,assminfo-level,assminfo-refseq-category
            assembly_to_tax[assembly_id] = taxid

    path_fasta_out = Path(args.prefix + '.fasta.gz')
    path_seqid_to_taxid = Path(args.prefix + '.seqid2taxid.tsv')

    path_fasta_out.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(path_fasta_out, 'wt') as f_fasta_out, path_seqid_to_taxid.open('w') as f_seqid2taxid:
        for assembly_path in seq_path.iterdir():
            if assembly_path.is_dir():
                assembly_id = assembly_path.name
                taxid = assembly_to_tax.get(assembly_id, args.kingdom_taxid)
                for seqfile in assembly_path.glob('*.fna'):
                    for record in SeqIO.parse(seqfile, 'fasta'):
                        SeqIO.write(record, f_fasta_out, 'fasta')
                        f_seqid2taxid.write(f'{record.id}\t{taxid}\n')


if __name__ == '__main__':
    main()
