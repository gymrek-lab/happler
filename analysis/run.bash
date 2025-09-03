#!/usr/bin/env bash
#SBATCH --export ALL
#SBATCH --partition condo
#SBATCH --account ddp268
#SBATCH --qos condo
#SBATCH --job-name happler-smk
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 0:30:00
#SBATCH --mem 8G
#SBATCH --output /dev/null
#SBATCH --signal=B:SIGUSR1@5


# An example bash script demonstrating how to run the entire snakemake pipeline
# This script creates a log file in the execution directory

# Before running this snakemake pipeline, remember to complete the config file
# with the required input info.
# Also, make sure that this script is executed from the directory that it lives in!

# Before we do anything, let's ensure that we get notified in case this job times out
trap "command -v slack &>/dev/null && slack \"TIMEOUT: snakemake job\"" SIGUSR1

out_path="out"
mkdir -p "$out_path"

# Clear leftover log files
echo ""> "$out_path/log"

# try to find and activate the snakemake conda env if we need it
if ! command -v 'snakemake' &>/dev/null && \
	command -v 'conda' &>/dev/null && \
   [ "$CONDA_DEFAULT_ENV" != "snakemake" ] && \
   conda info --envs | grep "$CONDA_ROOT/snakemake" &>/dev/null; then
        echo "Snakemake not detected. Attempting to switch to snakemake environment." >> "$out_path/log"
        eval "$(conda shell.bash hook)"
        conda activate snakemake
fi


# Check: Are we within a SLURM batch job or just an interactive node?
if [ "$ENVIRONMENT" = "BATCH" ]; then
    snakemake \
    --workflow-profile profile/slurm \
    --rerun-trigger {mtime,params,input} \
    --notemp \
    -k \
    -j 64 \
    -c 64 \
    "$@" &>"$out_path/log" &
else
    snakemake \
    --workflow-profile profile/default \
    --rerun-trigger {mtime,params,input} \
    --notemp \
    -k \
    -c 4 \
    "$@" &>"$out_path/log" &
fi

wait $!
exit_code=$?

# Send a slack message to notify us that the job completed
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
        slack "snakemake finished successfully" &>/dev/null
    else
        slack "snakemake happler-smk job failed" &>/dev/null
        slack "$(tail -n4 "$out_path/log")" &>/dev/null
    fi
fi
exit "$exit_code"
