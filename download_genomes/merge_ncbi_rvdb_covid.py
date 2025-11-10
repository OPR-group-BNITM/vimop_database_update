"""Extract Covid from RVDB

From a downloaded RVDB data set extract only the covid sequences and store them in 
a fasta file with header according to vimop input requirements.
"""


import argparse
import gzip
import io

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Entrez


def gzip_open(fname):
    if fname.endswith('.gz'):
        return gzip.open(fname, 'rt')
    return open(fname)


def read_taxinfo(fname):
    acc2tax = {}
    acc2fam = {}
    acc2spec = {}
    with open(fname, 'r') as f:
        for line in f:
            if not line.strip() or line.strip().startswith('accession'):  # skip header/empty
                continue
            cols = line.rstrip('\n').split('\t')
            acc = cols[0].strip()
            tax = cols[1].strip()
            family = cols[4].strip()
            species = cols[6].strip()
            if acc and tax:
                acc2tax[acc] = tax
                acc2fam[acc] = family
                acc2spec[acc] = species
    return acc2tax, acc2fam, acc2spec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--email', required=True, help='email used for entrez to download refseq covid sequence')
    parser.add_argument('--covid-refseq-id', help='Covid Refseq sequence is included', default='NC_045512.2')
    parser.add_argument('--taxinfo', required=True, help="TSV file with taxonomic info.")
    parser.add_argument('--rvdb', required=True, help="RVDB input fasta file.")
    parser.add_argument('--ncbi', required=True, help="NCBI virus input fasta file.")
    parser.add_argument('--fasta-out', required=True, help="Output fasta file.")
    parser.add_argument('--seqid-to-taxid-out', required=True, help="Output seqid to taxid map (.tsv file).")
    parser.add_argument('--virus-taxid', help='TaxID to use as default for all missing entries.', type=int, default=10_239)
    parser.add_argument('--covid-taxid', help='TaxID for covid sequences', type=int, default=2_697_049) 
    args = parser.parse_args()

    family = 'Coronaviridae'
    covid = 'Severe acute respiratory syndrome coronavirus 2'
    covid_taxid = args.covid_taxid
    viruses_taxid = args.virus_taxid

    taxids, acc2families, acc2species = read_taxinfo(args.taxinfo)

    with open(args.fasta_out, 'w') as f_fasta, open(args.seqid_to_taxid_out, 'w') as f_seqidtotaxid:

        # Step 1: Download the RefSeq covid entry

        Entrez.email = args.email
        with Entrez.efetch(db='nuccore', id=args.covid_refseq_id, rettype='fasta', retmode='text') as handle:
            record = SeqIO.read(io.StringIO(handle.read()), 'fasta')
            refseq_out = SeqRecord(
                id=record.id,
                description= '|'.join(['', record.description, family, covid]),
                seq=record.seq
            )
            SeqIO.write(refseq_out, f_fasta, 'fasta')
            f_seqidtotaxid.write(f'{record.id}\t{covid_taxid}\n')

        # Step 2: Add all RVDB covid genomes

        with gzip_open(args.rvdb) as f_rvdb:
            for record in SeqIO.parse(f_rvdb, 'fasta'):
                descr_split = record.description.split('|')
                if len(descr_split) < 5:
                    continue
                seqid = descr_split[2]
                description = descr_split[3]
                species = descr_split[4]
                if species == covid and seqid != args.covid_refseq_id:
                    info_string = '|'.join(['', description, family, covid])
                    seq_out = SeqRecord(
                        id=seqid,
                        description=info_string,
                        seq=record.seq
                    )
                    SeqIO.write(seq_out, f_fasta, 'fasta')
                    f_seqidtotaxid.write(f'{seqid}\t{covid_taxid}\n')

        # Step 3: merge with all other species genomes from NCBI virus

        with gzip_open(args.ncbi) as f_ncbi:
            for ncbi_record in SeqIO.parse(f_ncbi, 'fasta'):
                taxid = taxids.get(ncbi_record.id, viruses_taxid)
                fam = acc2families.get(ncbi_record.id, '')
                species = acc2species.get(ncbi_record.id, '')
                new_record = SeqRecord(
                    seq=ncbi_record.seq,
                    id=ncbi_record.id,
                    description=f'|{ncbi_record.description}|{fam}|{species}'
                )
                SeqIO.write(new_record, f_fasta, 'fasta')
                f_seqidtotaxid.write(f'{ncbi_record.id}\t{taxid}\n')


if __name__ == '__main__':
    main()
