#!/usr/bin/env bash
set -euo pipefail

# archive_and_split.sh
#
# Usage:
#   archive_and_split.sh <datadir> <db_name: virus|centrifuge|contaminants> <output_dir> [url_prefix]

datadir=${1:-}
db_name=${2:-}
output_dir=${3:-}
url_prefix=${4:-}
chunk_size="200M"

# -------- Input validation --------
usage() {
  echo "Usage: $0 <datadir> <db_name: virus|centrifuge|contaminants> <output_dir> [url_prefix]" >&2
  exit 2
}

[[ -n "$datadir" && -n "$db_name" && -n "$output_dir" ]] || usage

if [[ ! -d "$datadir" ]]; then
  echo "ERROR: datadir '$datadir' is not a directory or not accessible." >&2
  exit 1
fi

case "$db_name" in
  virus|centrifuge|contaminants) ;;
  *)
    echo "ERROR: db_name must be one of: virus, centrifuge, contaminants (got '$db_name')." >&2
    exit 1
    ;;
esac

fn_config_yaml="${datadir}/${db_name}.yaml"
if [[ ! -f "$fn_config_yaml" ]]; then
  echo "ERROR: Expected YAML config not found: ${fn_config_yaml}" >&2
  exit 1
fi

if [[ -z "${url_prefix}" ]]; then
  echo "ERROR: url_prefix not provided (4th arg) and URL_PREFIX env var not set." >&2
  exit 1
fi

# Ensure required tools exist
for tool in grep awk sed sha256sum sort xz tar split; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: Required tool '$tool' not found in PATH."; exit 1; }
done

# -------- Get version from YAML --------
# Accepts lines like:
#   version: 1.2
#   version: '1.2'
#   version: "1.2"
version="$(
  grep -E '^[[:space:]]*version[[:space:]]*:' "$fn_config_yaml" \
  | head -n1 \
  | sed -E 's/^[[:space:]]*version[[:space:]]*:[[:space:]]*//; s/[[:space:]]*$//' \
  | tr -d '"'\'''
)"

if [[ -z "$version" ]]; then
  echo "ERROR: Could not extract 'version' from ${fn_config_yaml}." >&2
  exit 1
fi

# Check version format X.Y (e.g., 1.0 or 2.1)
if [[ ! "$version" =~ ^[0-9]+\.{1}[0-9]+$ ]]; then
  echo "ERROR: Version '$version' does not match required format X.Y (e.g., 1.0, 2.1)." >&2
  exit 1
fi

versioned_db="${db_name}.v${version}"

# -------- Compute checksum of the whole directory --------

checksum_whole_directory=$(
    find "$datadir/." -type f -exec sha256sum {} \; \
    | sort \
    | awk '{print $1}' \
    | sha256sum \
    | awk '{print $1}'
)

# -------- Compress the directory --------
mkdir -p "$output_dir/zipped_archives"
fname_zipped="$output_dir/zipped_archives/$versioned_db.tar.xz"

if [[ -f "$fname_zipped" ]]; then
  echo "ERROR: Zipped archive ${fname_zipped} already exists." >&2
  exit 1
fi

tar \
  --use-compress-program="xz -9e -T0" \
  -cvf "$fname_zipped" \
  -C "$(dirname "$datadir")" "$(basename "$datadir")"

# -------- Split the zipped file --------
split_dir="$output_dir/split/files"
mkdir -p "$split_dir"
split_prefix="${versioned_db}.tar.xz.part."
split -b "$chunk_size" "$fname_zipped" "$split_dir/$split_prefix"

# -------- Create the files YAML with checksums --------
fname_file_config="${versioned_db}.files.yaml"
cd "$split_dir"

checksum_zipped="$(sha256sum "$fname_zipped" | awk '{print $1}')"

{
  echo "${db_name}:"
  echo "  version: '${version}'"
  echo "  checksum_directory: '${checksum_whole_directory}'"
  echo "  checksum_zipped: '${checksum_zipped}'"
  echo "  files:"
  for fn in ${split_prefix}*; do
    [[ -f "$fn" ]] || { echo "ERROR: No split parts found with prefix '${split_prefix}'." >&2; exit 1; }
    checksum_part="$(sha256sum "$fn" | awk '{print $1}')"
    echo "  - name: '${fn}'"
    echo "    url: '${url_prefix}/${fn}'"
    echo "    checksum: '${checksum_part}'"
  done
} > "$fname_file_config"

echo "OK: Created archive, split parts, and ${split_dir}/${fname_file_config}"
