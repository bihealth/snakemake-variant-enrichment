#!/bin/bash

# Set a name for the job (-J or --job-name).
#SBATCH --job-name=enrichment

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=logs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=2

# Force allocation of the two cores on ONE node.
#SBATCH --nodes=1

# Set the total memory. Units can be given in T|G|M|K.
#SBATCH --mem=1G

# Optionally, set the partition to be used (-p or --partition).
#SBATCH --partition=critical

# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=7-00:00

set -x


export TMPDIR=/fast/users/${USER}/scratch/tmp
export LOGDIR=logs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
mkdir -p $LOGDIR

unset DRMAA_LIBRARY_PATH
eval "$($(which conda) shell.bash hook)"
conda activate snakemake-variant-enrichment

JOBS=100

# Note that Slurm DRMAA differs slightly from original Slurm syntax
# --mem-per-cpu doesn't accept units and the default unit here is MB
# -t only accepts HH:MM
snakemake \
    --jobscript bih_cluster_jobscript.sh \
    --use-conda \
    --cluster-config bih_cluster_config.json \
    --drmaa " \
        -p critical \
        -t {cluster.time} \
        --nodes=1 \
        --mem={cluster.memory} \
        --cpus-per-task={cluster.threads} \
        -o $LOGDIR/%x-%j.log" \
    -j $JOBS \
    -k \
    -p \
    $@

