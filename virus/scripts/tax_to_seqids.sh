#!/usr/bin/env bash

set -euo pipefail
set -x

# Setup taxon kit
# - install: conda install taxonkit)
# - cd /Users/nils.petersen/.taxonkit
# - wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
# - tar -xzf taxdump.tar.gz

virus_taxid=10239

outdir="taxids"

mkdir -p "$outdir"

fname_seqs="$HOME/ViMOP_DB/virus/ALL.fasta"
fname_headers="${outdir}/headers.txt"
fname_virus_species="${outdir}/species.txt"
fname_virus_families="${outdir}/families.txt"
fname_taxids_species="${outdir}/species.taxids.tab"
fname_taxids_families="${outdir}/families.taxids.tab"
fname_seqid2taxid="${outdir}/seqid2taxid.tsv"

awk -F'|' '/^>/{ sub(/^>/, "", $1); print $1"|"$(NF-3)"|"$(NF-2) }' $fname_seqs > $fname_headers

awk -F'|' '{ f=$2; sub(/^[[:space:]]+|[[:space:]]+$/, "", f); if (!seen[f]++) print f }' $fname_headers > $fname_virus_families
awk -F'|' '{ f=$3; sub(/^[[:space:]]+|[[:space:]]+$/, "", f); if (f!="" && !seen[f]++) print f }' $fname_headers > $fname_virus_species

taxonkit name2taxid --threads 8 < $fname_virus_species > $fname_taxids_species
taxonkit name2taxid --threads 8 < $fname_virus_families > $fname_taxids_families

awk -v FS='|' -v OFS='\t' -v virus_taxid="$virus_taxid" \
    -v sp_map="$fname_taxids_species" -v fam_map="$fname_taxids_families" '
function trim(s){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }

BEGIN {
  # load species map (tab-separated: name \t taxid \t rank)
  while ((getline line < sp_map) > 0) {
    n = split(line, a, "\t")
    if (n >= 2) {
      name = trim(a[1]); tax = trim(a[2])
      if (tax != "" && tax != "NA") sp[name] = tax
    }
  }
  close(sp_map)

  # load family map (tab-separated)
  while ((getline line < fam_map) > 0) {
    n = split(line, a, "\t")
    if (n >= 2) {
      name = trim(a[1]); tax = trim(a[2])
      if (tax != "" && tax != "NA") fam[name] = tax
    }
  }
  close(fam_map)
}
# Process every line of headers.txt (format: seqid|family|species) — NO /^>/
{
  id      = trim($1)
  famname = trim($2)
  spname  = trim($3)

  tax = ""
  if (spname in sp)        tax = sp[spname]
  else if (famname in fam) tax = fam[famname]
  else                     tax = virus_taxid

  print id, tax
}' "$fname_headers" > "$fname_seqid2taxid"


rm -f $fname_headers $fname_virus_species $fname_virus_families $fname_taxids_species $fname_taxids_families
