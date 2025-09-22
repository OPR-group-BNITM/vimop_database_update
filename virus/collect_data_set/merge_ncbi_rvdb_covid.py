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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--email', help='email used for entrez to download refseq covid sequence')
    parser.add_argument('--covid-refseq-id', help='Covid Refseq sequence is included', default='NC_045512.2')
    parser.add_argument('-i', '--input', required=True, help="RVDB input fasta file.")
    parser.add_argument('-o', '--output', required=True, help="Output fasta file.")
    args = parser.parse_args()

    family = 'Coronaviridae'
    covid = 'Severe acute respiratory syndrome coronavirus 2'

    with open(args.output, 'w') as f_out:
        Entrez.email = args.email
        with Entrez.efetch(db='nuccore', id=args.covid_refseq_id, rettype='fasta', retmode='text') as handle:
            record = SeqIO.read(io.StringIO(handle.read()), 'fasta')
            refseq_out = SeqRecord(
                id=record.id,
                description= '|'.join(['', record.description, family, covid]),
                seq=record.seq
            )
            SeqIO.write(refseq_out, f_out, 'fasta')

        with gzip_open(args.input) as f_in:
            for record in SeqIO.parse(f_in, 'fasta'):
                descr_split = record.id.split('|')
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
                    SeqIO.write(seq_out, f_out, 'fasta')


if __name__ == '__main__':
    main()
