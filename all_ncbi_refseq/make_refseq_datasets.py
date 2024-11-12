import pandas as pd
import argparse
import os
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('headers_tsv')
    parser.add_argument('fasta')
    parser.add_argument('outdir')
    args = parser.parse_args()

    headers = pd.read_csv(
        args.headers_tsv,
        sep='\t',
        header=None,
        names=['id', 'description', 'family', 'species']
    )

    print("Loading data and extracting species with refseqs.")

    headers = headers[['id', 'species']]
    refseqs = headers[headers['id'].str.startswith(('AC_', 'NC_'))][['species', 'id']]

    species_found = set(refseqs['species'])

    sequences_with_a_refseq = headers[headers['species'].isin(species_found)]

    print('Total entries', headers.shape[0])
    print('Refseq entries:', refseqs.shape[0])
    print('Entries of species that do have at least one RefSeq:', sequences_with_a_refseq.shape[0])

    os.makedirs(args.outdir, exist_ok=True)

    species_count = []

    # Iterate over unique species
    for i, (species, group) in enumerate(sequences_with_a_refseq.groupby('species'), 0):

        if i % 1000 == 0:
            print(f'Writing output for species number {i+1} of ({len(species_found)}) ({species})')

        species_dir_path = os.path.join(args.outdir, species.replace(" ", "_"))

        os.makedirs(species_dir_path, exist_ok=True)

        # Write sequence ids and refseq ids to respective files
        with open(os.path.join(species_dir_path, 'ids.txt'), 'w') as ids_file:
            ids_file.write("\n".join(group['id']))

        with open(os.path.join(species_dir_path, 'refseq_ids.txt'), 'w') as refseq_ids_file:
            refseqs_species = refseqs[refseqs['species'] == species]
            refseq_ids_file.write("\n".join(refseqs_species['id']))

        species_count.append((species, group.shape[0]))

    with open(os.path.join(args.outdir, 'species.tsv'), 'w') as f_species:
        f_species.write('\n'.join(f'{spec}\t{c}' for spec, c in species_count) + '\n')

    print('Writing fasta')

    for i, record in enumerate(SeqIO.parse(args.fasta, "fasta")):
        if i % 1000 == 0:
            print(f'Writing output for fasta record {i+1} of ({headers.shape[0]})')
        species = record.description.rsplit("|", 1)[1].strip()
        if species in species_found:
            with open(os.path.join(args.outdir, species.replace(" ", "_"), 'sequences.fasta'), 'a') as fasta_out:
                SeqIO.write(record, fasta_out, 'fasta')


if __name__ == '__main__':
    main()
