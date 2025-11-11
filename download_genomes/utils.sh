#!/usr/bin/env bash


run_with_retries() {
  local max_tries=${1:-8}; shift
  local try=1
  local delay=2
  while true
  do
    if "$@"
    then
      return 0
    fi
    if (( try >= max_tries ))
    then
      echo "ERROR: command failed after ${max_tries} attempts: $*" >&2
      return 1
    fi
    echo "warn: transient failure (attempt $try). retrying in ${delay}s…" >&2
    sleep "$delay"
    delay=$(( delay * 2 ))
    try=$(( try + 1 ))
  done
}
