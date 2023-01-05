#!/usr/bin/env bash

output_dir=$(grep -E '^output_dir: ' config/config.yml | cut -f2 -d' ' | tr -d "'")

snakemake --use-singularity --singularity-args "--bind $output_dir:/data" --use-conda "$@"
