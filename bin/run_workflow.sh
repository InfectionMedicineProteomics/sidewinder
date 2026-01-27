#!/usr/bin/env bash

export XDG_CACHE_HOME=/srv/data1/jo0348st/.cache

export TMPDIR=/srv/data1/jo0348st/.tmp

export CUDA_VILISBLE_DEVISES=0,1

# Add -k to continue with independant jobs.
snakemake \
    --cores 46 \
    --use-conda \
    --use-singularity \
    --singularity-args "--nv --bind /srv/data1" \
    --resources nvidia_gpu=2 \
    --rerun-triggers mtime \
    "$@"
