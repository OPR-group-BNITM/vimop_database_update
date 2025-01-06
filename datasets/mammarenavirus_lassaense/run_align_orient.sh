#!/usr/bin/env bash

set -x

nextflow ../../shared/workflow/align_to_refs.nf -c align_orient.nextflow.config -with-conda
