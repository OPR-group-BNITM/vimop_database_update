import yaml
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input_files',
        nargs='+',
        help='Input files from NCBI virus to extract families from.'
    )
    parser.add_argument(
        '--out',
        help='Yaml file to write output to',
        default='families.yaml'
    )
    args = parser.parse_args()

    families = {}

    for fname in args.input_files:
        for i, seq in enumerate(SeqIO.parse(fname, 'fasta')):
            if i % 10_000 == 0:
                print(i)
            _, family, organism = map(str.strip, seq.description.rsplit('|', 2))
            if not organism:
                print(f'WARGNING: No organism in {seq.description}')
                continue
            if organism in families and families[organism] != family:
                print(f'WARNING: Found different families {families[organism]} and {family} for organism {organism}')
                print(seq.description)
            else:
                families[organism] = family

    with open(args.out, 'w') as f_out:
        yaml.dump(families, f_out, default_flow_style=False, sort_keys=False)
    

if __name__ == '__main__':
    main()
