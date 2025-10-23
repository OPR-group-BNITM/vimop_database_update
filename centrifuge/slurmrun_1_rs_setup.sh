#!/usr/bin/env bash

set -euo pipefail

job_name=rs_setup
container=vimop_db_update_general

this_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
source "${this_script_dir}/../env.sh"
jobscript="${this_script_dir}/${job_name}.sh"

mkdir -p "${VIMOP_DB_UPDATE_SRC_SLURM_OUTPUT_DIR}"

sbatch \
    --container=$HOME/containers/${container} \
    --job-name=$job_name \
    -o "${VIMOP_DB_UPDATE_SRC_SLURM_OUTPUT_DIR}/slurm.${job_name}.%N.%j.out" \
    -e "${VIMOP_DB_UPDATE_SRC_SLURM_OUTPUT_DIR}/slurm.${job_name}.%N.%j.err" \
    --mail-user=$VIMOP_DB_UPDATE_USER_MAIL \
    ${jobscript}
