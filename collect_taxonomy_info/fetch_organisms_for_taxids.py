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

    out = {}

    # category is e.g. curated, filters
    for category, category_dict in config.items():
        # identifier like DENV, ARENA
        for identifier, groupdict in category_dict.items():
            for taxon in groupdict['taxa']:
                organisms = fetch_organism_names(taxon['taxid'])
                taxon_dict_out = taxon | {'organisms': organisms}
                
                out.setdefault(category, {}).setdefault(
                    identifier, {}).setdefault('taxa', []).append(taxon_dict_out)

    with open(args.out, "w") as f_out:
        yaml.dump(dict(out), f_out, default_flow_style=False, sort_keys=False)


if __name__ == '__main__':
    main()
