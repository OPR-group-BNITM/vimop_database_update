from Bio import Entrez
#import xml.etree.ElementTree as ET

# Set your email here
#Entrez.email = "your_email@example.com"

# Replace XXXX with your taxon ID
taxon_id = "11266"

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
