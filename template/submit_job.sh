#!/usr/bin/env bash

casename=dummy # replace (e.g. with 'covid')

/opt/gaia-tools/bin/queue_job.sh curator "../../shared/job.sh ALL $casename"
