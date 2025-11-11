#!/usr/bin/env bash
#SBATCH -c 1 # number of cores
#SBATCH --mem 16000 # memory pool for all cores
#SBATCH --mail-type=ALL

set -euo pipefail
set -x

source "$(/opt/conda/bin/conda info --base)/etc/profile.d/conda.sh"
conda activate datasets

source "${VIMUPDATE_SRC}/download_genomes/utils.sh"

virus_taxid="${VIMUPDATE_VIRUS_TAXID}"

OUT_DIR="${VIMUPDATE_GENOMES}/ncbi_virus"
mkdir -p "$OUT_DIR"

TMP_DIR="${VIMUPDATE_GENOMES}/tmp_virus_taxinfo"
mkdir -p "$TMP_DIR"

# SETUP taxonkit
WORKDIR=$(pwd)
mkdir -p $VIMUPDATE_TAXONKIT
cd "$VIMUPDATE_TAXONKIT"
wget "${VIMUPDATE_TAXDUMP_URL}"
tar -xzf taxdump.tar.gz
cd "$WORKDIR"

# get all virus taxids
echo "$virus_taxid" \
| taxonkit list --data-dir "$VIMUPDATE_TAXONKIT" \
| sed 's/.* //' \
> "${OUT_DIR}/virus_taxids.txt"

# helper: from a tax-summary TSV on stdin, print the Taxid column (by header name)
taxsummary_to_ids() {
  awk -F'\t' '
    NR==1 {
      for (i=1;i<=NF;i++) if ($i=="Taxid") {c=i; break}
      if (!c) { print "ERROR: Taxid column not found in header" > "/dev/stderr"; exit 1 }
      next
    }
    { print $c }
  '
}

# 1) Get all families to use as chunks
FAM_ALL="$TMP_DIR/families.taxsummary.tsv"
FAM_IDS="$TMP_DIR/families.ids.txt"

if [[ ! -s "$FAM_ALL" ]]
then
  run_with_retries 8 bash -c "
    datasets summary taxonomy taxon ${virus_taxid} --children --rank family --as-json-lines \
    | dataformat tsv taxonomy --template tax-summary
  " > "$FAM_ALL"
fi

taxsummary_to_ids < "$FAM_ALL" > "$FAM_IDS"

# 2) Summarize NCBI Virus entries chunk-by-chunk with retries
summarize_taxon() {
  local tx="$1"
  run_with_retries 8 bash -c "
    datasets summary virus genome taxon '$tx' --as-json-lines \
    | dataformat tsv virus-genome --fields accession,virus-tax-id,virus-name,sourcedb \
    | tail -n +2
  "
}

OUT_TSV="${TMP_DIR}/viruses.taxsummary.tsv"
{
  printf "accession\tvirus-tax-id\tvirus-name\tsourcedb"
  while read -r tx
  do
    [[ -n "$tx" ]] || continue
    echo "Summarizing family taxid $tx ..." >&2
    summarize_taxon "$tx"
  done < "$FAM_IDS"
} > "$OUT_TSV"

UNIQUE_VIRUS_TAX_IDS="$TMP_DIR/unique.virus.ids.txt"
VIRUS_TAXID_TO_SPECIES_AND_FAM="$TMP_DIR/virusid.to.family.tsv"

awk -F'\t' 'NR>1 && !seen[$2]++ { print $2 }' "$OUT_TSV" > "$UNIQUE_VIRUS_TAX_IDS"

{
  echo -e "virus-tax-id\tfamily\tfamily-tax-id\tspecies\tspecies-tax-id"
  taxonkit reformat --data-dir "$VIMUPDATE_TAXONKIT" -I 1 -f '{f}|{s}' -t "$UNIQUE_VIRUS_TAX_IDS" \
  | awk -F '\t' 'BEGIN{OFS="\t"}{
    split($2,n,"|");  # names: family|species
    split($3,i,"|");  # ids:   family|species
    print $1, n[1], i[1], n[2], i[2]
  }'
} > "$VIRUS_TAXID_TO_SPECIES_AND_FAM"

# merge

VIRUS_TABLE_MERGED="$OUT_DIR/viruses.tsv"

awk -F'\t' 'NR==FNR {
  # build map from virusid.to.family.tsv (skip header)
  if (FNR>1) m[$1]=$2 FS $3 FS $4 FS $5
  next
}
FNR==1 {
  # header from viruses.taxsummary.tsv + extra columns
  print $0 FS "family\tfamily-tax-id\tspecies\tspecies-tax-id"
  next
}
/^Summarizing/ { next }  # drop that summary line
{
  k=$2
  if (k in m) print $0 FS m[k]
  else        print $0 FS "NA\tNA\tNA\tNA"
}' "$VIRUS_TAXID_TO_SPECIES_AND_FAM" "$OUT_TSV" \
> $VIRUS_TABLE_MERGED

# Filter out covid

VIRUS_TABLE_NOCOVID="$OUT_DIR/viruses.nocovid.tsv"

awk -F"\t" -v covid_taxid="$VIMUPDATE_COVID_TAXID" ' $2 != covid_taxid { print $0 }' "$VIRUS_TABLE_MERGED" > "$VIRUS_TABLE_NOCOVID"

echo "Done. TSV: $VIRUS_TABLE_NOCOVID (No covid) and $VIRUS_TABLE_MERGED (including covid)"
