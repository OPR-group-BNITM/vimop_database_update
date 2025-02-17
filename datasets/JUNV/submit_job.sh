#!/usr/bin/env bash

casename=JUNV

/opt/gaia-tools/bin/queue_job.sh nextflow "../../shared/job.sh ALL $casename"
