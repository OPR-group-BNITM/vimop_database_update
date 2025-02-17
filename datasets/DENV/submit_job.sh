#!/usr/bin/env bash

casename=DENV

/opt/gaia-tools/bin/queue_job.sh nextflow "../../shared/job.sh ALL $casename"
