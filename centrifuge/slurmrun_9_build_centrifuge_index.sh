#!/usr/bin/env bash

set -euo pipefail

job_name=build_centrifuge_index
container=vimop_db_update_general

this_script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
source "${this_script_dir}/../env.sh"
jobscript="${this_script_dir}/${job_name}.sh"

mkdir -p "${VIMUPDATE_SLURM}"

sbatch \
    --container=$HOME/containers/${container} \
    --job-name=$job_name \
    -o "${VIMUPDATE_SLURM}/slurm.${job_name}.%N.%j.out" \
    -e "${VIMUPDATE_SLURM}/slurm.${job_name}.%N.%j.err" \
    --mail-user=$VIMUPDATE_MAIL \
    ${jobscript}
