#!/usr/bin/env bash

#output_dir=$(grep -E '^output_dir: ' config/user_config.yml | cut -f2 -d' ' | tr -d "'")

# Removed from command:
# --singularity-args "--bind $output_dir:/data"

# Add -k to continue with independant jobs.
snakemake --use-singularity --singularity-args "--bind /srv/data1" --use-conda --rerun-triggers mtime --cores 150 "$@"
