import argparse
from Bio import Entrez

# Set your email here
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

    # Print taxon IDs and names
    for record in records:
        tax_id = record["TaxId"]
        name = record["ScientificName"]
        print(f"{tax_id}\t{name}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--taxid', default="11266")
    args = parser.parse_args()

    fetch_organism_names(args.taxid)


if __name__ == '__main__':
    main()
