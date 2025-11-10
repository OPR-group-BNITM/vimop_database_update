#!/usr/bin/env bash


container=vimop_db_update_general
source "env.sh"


# submit_job <job_name> <job_dir> [dep_jobid ...]
submit_job() {
  local job_name="$1"
  local job_dir="$2"
  shift 2
  local deps=("$@")

  local jobscript="${job_dir}/${job_name}.sh"

  # Basic checks
  if [ -z "$job_name" ] || [ -z "$job_dir" ]; then
    echo "usage: submit_job <job_name> <job_dir> [dep_jobid ...]" >&2
    return 2
  fi
  if [ ! -f "$jobscript" ]; then
    echo "error: job script not found: $jobscript" >&2
    return 2
  fi

  # Build dependency string if any deps were provided (afterok)
  local dep_arg=()
  if [ "${#deps[@]}" -gt 0 ]; then
    local depstr="afterok:${deps[0]}"
    for d in "${deps[@]:1}"; do
      depstr="${depstr}:$d"
    done
    dep_arg=(--dependency="$depstr")
  fi

  # Submit with --parsable so stdout is the JobID; echo it so caller can capture
  local jid
  jid=$(
    sbatch \
      --parsable \
      --container="$HOME/containers/${container}" \
      --container-security=seccomp=unconfined \
      --job-name="$job_name" \
      -o "${VIMUPDATE_SLURM}/slurm.${job_name}.%N.%j.out" \
      -e "${VIMUPDATE_SLURM}/slurm.${job_name}.%N.%j.err" \
      --mail-user="$VIMUPDATE_MAIL" \
      "${dep_arg[@]}" \
      "$jobscript"
  ) || return $?

  echo "$jid"
}


mkdir -p "${VIMUPDATE_SLURM}"

# Downloads
id_refseq_download=$(submit_job rs_setup download_genomes)
id_virus_taxinfo=$(submit_job virus_taxinfo download_genomes)
id_virus_sequences=$(submit_job virus_sequences download_genomes "$id_virus_taxinfo")
id_rvdb=$(submit_job rvdb download_genomes)
id_merge_ncbi_rvdb_covid=$(submit_job merge_ncbi_rvdb_covid download_genomes "$id_virus_sequences" "$id_rvdb")
id_rna=$(submit_job rna download_genomes)
id_contaminants=$(submit_job contaminants download_genomes)

# Database setup
id_virus_db=$(submit_job virus_db virus "$id_merge_ncbi_rvdb_covid")
id_host_db=$(submit_job contaminant_db host "$id_rna" "$id_refseq_download" "$id_contaminants")
id_centrifuge_db=$(submit_job build_centrifuge_index centrifuge "$id_virus_db" "$id_refseq_download")

# Package
submit_job ship shipit "$id_virus_db" "$id_host_db" "$id_centrifuge_db"
