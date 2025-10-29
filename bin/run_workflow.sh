#!/usr/bin/env bash

export XDG_CACHE_HOME=/srv/data1/jo0348st/.cache

export TMPDIR=/srv/data1/jo0348st/.tmp

# Add -k to continue with independant jobs.
snakemake \
    --cores 46 \
    --use-conda \
    --use-singularity \
    --singularity-args "--nv --bind /srv/data1" \
    --resources nvidia_gpu=4 \
    --rerun-triggers mtime \
    "$@"
