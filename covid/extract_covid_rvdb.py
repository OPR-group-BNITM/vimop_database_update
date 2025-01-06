from Bio import SeqIO

# Input and output file paths
input_fasta = "/Users/nils.petersen/dev/rvdb/C-RVDBv29.0.fasta"
output_fasta = "/Volumes/DataCurate/CurationDatasets/Covid/CovidRVDB/covid_rvdb_v29.fasta"

# Open the output FASTA file for writing
with open(output_fasta, "w") as output_handle:
    # Parse the input FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        # Check if the description contains "Severe acute respiratory syndrome coronavirus 2"
        if "severe acute respiratory syndrome coronavirus 2" in record.description.lower():
            # Extract the GenBank ID and reformat the header
            parts = record.description.split("|")
            genbank_id = parts[2] if len(parts) > 2 else "Unknown"
            description = parts[3] if len(parts) > 3 else ""
            # Create a new record with the updated description
            record.id = genbank_id
            record.description = f"|{description}|Coronaviridae|Severe acute respiratory syndrome coronavirus 2"

            # Write the updated record to the output file
            SeqIO.write(record, output_handle, "fasta")

print(f"Filtered sequences have been written to {output_fasta}")
