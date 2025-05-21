"""Add organism names to a yaml file with taxonomic groups.

The file is organised like this

curated:
  LASV:
    name: Mammarenavirus lassense
    taxa:
    - taxid: 3052310
    segments:
      S:
        refs:
        - NC_004296.1
        seqs:
        - NC_004296.1
        - KM822128.1
        - GU481068.1
    ...

And it will will add an entry with organism names for each taxid (here LASV) with organism names.

  LASV:
    name: Mammarenavirus lassense
    taxa:
    - taxid: 3052310
      organisms:
      - Mammarenavirus lassaense
      - Lassa virus Josiah
      - Lassa virus GA391
    segments:
      S:
        refs:
        - NC_004296.1
        seqs:
        - NC_004296.1
        - KM822128.1
        - GU481068.1
    ...

"""
import copy
import argparse
import yaml
from Bio import Entrez


Entrez.email = 'nils.petersen@bnitm.de'


def fetch_organism_names(taxon_id):
    # Use Entrez esearch and efetch to get child taxa
    handle = Entrez.esearch(db="taxonomy", term=f"txid{taxon_id}[Subtree]", retmax=10000)
    record = Entrez.read(handle)
    handle.close()

    # Retrieve taxon summaries
    id_list = record["IdList"]
    handle = Entrez.efetch(db="taxonomy", id=",".join(id_list), retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    return [
        str(record["ScientificName"])
        for record in records
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', help='yaml config file', required=True)
    parser.add_argument('--out', help='yaml out file', required=True)
    args = parser.parse_args()

    with open(args.config) as f_in:
        config = yaml.safe_load(f_in)

    out = copy.deepcopy(config)

    # category is e.g. curated, filters
    for category, category_dict in config.items():
        # identifier like DENV, ARENA
        for identifier, groupdict in category_dict.items():
            out[category][identifier]['taxa'] = [
                taxon | {'organisms': fetch_organism_names(taxon['taxid'])}
                for taxon in groupdict['taxa']
            ] 

    with open(args.out, "w") as f_out:
        yaml.dump(dict(out), f_out, default_flow_style=False, sort_keys=False)


if __name__ == '__main__':
    main()
