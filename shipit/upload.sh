#!/usr/bin/env bash

set -euo pipefail

#scp -r "${VIMUPDATE_DB}/split/files" "${VIMUPDATE_FILESHARE_PREFIX}/files"
scp -r "${VIMUPDATE_DB}/split/vimop_db.${VIMUPDATE_DB_VERSION}.yaml" "${VIMUPDATE_FILESHARE_PREFIX}/vimop_db.${VIMUPDATE_DB_VERSION}.yaml"
scp -r "${VIMUPDATE_DB}/split/vimop_db.${VIMUPDATE_DB_VERSION}.yaml" "${VIMUPDATE_FILESHARE_PREFIX}/latest.yaml"
