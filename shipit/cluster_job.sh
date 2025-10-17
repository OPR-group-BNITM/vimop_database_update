#!/usr/bin/env bash

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

db_name=virus
data_dir=/data2/le-petersen/shared_data/db_different_versions/VimopDatabaseUpdate/virus/data/output_databases/db2.3/virus
out_dir=/data2/le-petersen/shared_data/db_different_versions/

${script_dir}/archive_and_split.sh $data_dir $db_name $out_dir
