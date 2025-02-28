import yaml
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def main():
    fname_families = '/data/home/nils.petersen/dev/VirusDatasetCuration/data/out/organisms.yaml'
    fname_rvdb = '/data2/le-petersen/nils/base_datasets/rvdb/C-RVDBv29.0.fasta'
    fname_refseq = '/data2/le-petersen/nils/base_datasets/refseq/refseq_virus_20250227.fasta'
    fname_out = '/data2/le-petersen/nils/base_datasets/merged_datasets/refseq_virus_20250227_C-RVDBv29.0.fasta'

    with open(fname_families) as f_in:
        families = {
            species.upper(): family
            for species, family in yaml.safe_load(f_in).items()
        }

    with open(fname_out, 'w') as f_out:
        for seq in SeqIO.parse(fname_refseq, 'fasta'):
            SeqIO.write(seq, f_out, 'fasta')
        for seq in SeqIO.parse(fname_rvdb, 'fasta'):
            # >acc|GENBANK|HM118276.1|HIV-1 isolate BREPM2999_06 from Brazil envelope glycoprotein (env) gene, partial cds|Human immunodeficiency virus 1|VRL|25-JUL-2016
            if seq.description.strip().startswith('acc|REFSEQ'):
                continue
            try:
                _, _, genbank_id, therest =  map(str.strip, seq.description.split('|', 3))
                description, species, _, _ = therest.split('|')
            except ValueError:
                print(seq.description)
                continue
            family = families.get(species.upper(), '')
            SeqIO.write(
                SeqRecord(
                    id=seq.id,
                    description=f'{genbank_id}|{description}|{family}|{species}',
                    name=seq.name,
                    seq=seq.seq
                ),
                f_out,
                'fasta'
            )


if __name__ == '__main__':
    main()
