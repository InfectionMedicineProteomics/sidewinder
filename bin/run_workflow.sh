#!/usr/bin/env bash


snakemake -c 10 --use-singularity --use-conda "$@"
