#!/usr/bin/env bash

casename=LCMV

/opt/gaia-tools/bin/queue_job.sh nextflow "../../shared/job.sh ALL $casename"
