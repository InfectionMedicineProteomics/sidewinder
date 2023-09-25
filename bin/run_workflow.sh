#!/usr/bin/env bash

#output_dir=$(grep -E '^output_dir: ' config/user_config.yml | cut -f2 -d' ' | tr -d "'")

# Removed from command:
# --singularity-args "--bind $output_dir:/data"

snakemake --use-singularity --use-conda --cores 150 "$@"  # Add -k to continue with independant jobs.
