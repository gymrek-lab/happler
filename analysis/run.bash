#!/usr/bin/env bash
#SBATCH --export ALL
#SBATCH --partition hotel
#SBATCH --qos hotel
#SBATCH --job-name happler-smk
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 8:00:00
#SBATCH --output /dev/null

# An example bash script demonstrating how to run the entire snakemake pipeline
# This script creates two separate log files in the output dir:
# 	1) log - the basic snakemake log of completed rules

# Before running this snakemake pipeline, remember to complete the config file
# with the required input info.
# Also, make sure that this script is executed from the directory that it lives in!

out_path="out"
mkdir -p "$out_path"

# clear leftover log files
echo ""> "$out_path/log"

# try to find and activate the snakemake conda env if we need it
if ! command -v 'snakemake' &>/dev/null && \
	command -v 'conda' &>/dev/null && \
   [ "$CONDA_DEFAULT_ENV" != "snakemake" ] && \
   conda info --envs | grep "$CONDA_ROOT/snakemake" &>/dev/null; then
        echo "Snakemake not detected. Attempting to switch to snakemake environment." >> "out/log"
        eval "$(conda shell.bash hook)"
        conda activate snakemake
fi


# check: are we being executed from within qsub?
if [ "$ENVIRONMENT" = "BATCH" ]; then
    snakemake \
    --workflow-profile profile/slurm \
    --rerun-trigger {mtime,params,input} \
    -k \
    -j 16 \
    -c 16 \
    "$@" &>"$out_path/log"
else
    snakemake \
    --workflow-profile profile/default \
    --notemp \
    --rerun-trigger {mtime,params,input} \
    -k \
    -c 4 \
    "$@" &>"$out_path/log"
fi

exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
        slack "snakemake finished successfully" &>/dev/null
    else
        slack "snakemake simulate_gwas job failed" &>/dev/null
        slack "$(tail -n4 "$out_path/log")" &>/dev/null
    fi
fi
exit "$exit_code"
