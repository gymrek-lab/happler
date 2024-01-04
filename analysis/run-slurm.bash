#!/usr/bin/env bash
#SBATCH --export=ALL
#SBATCH --partitin=hotel
#SBATCH --qos=hotel
#SBATCH --output=/dev/null
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=8:00:00

# An example bash script demonstrating how to run the entire snakemake pipeline
# This script creates two separate log files in the output dir:
# 	1) log - the basic snakemake log of completed rules

# Before running this snakemake pipeline, remember to complete the config file
# with the required input info.
# Also, make sure that this script is executed from the directory that it lives in!

out_path="out"
mkdir -p "$out_path/logs"

# clear leftover log files
echo ""> "$out_path/logs/log"

# try to find and activate the snakemake conda env if we need it
if ! command -v 'snakemake' &>/dev/null && \
	command -v 'conda' &>/dev/null && \
   [ "$CONDA_DEFAULT_ENV" != "snakemake" ] && \
   conda info --envs | grep "$CONDA_ROOT/snakemake" &>/dev/null; then
        echo "Snakemake not detected. Attempting to switch to snakemake environment." >> "out/logs/log"
        eval "$(conda shell.bash hook)"
        conda activate snakemake
fi


# check: are we being executed from within qsub?
if [ "$ENVIRONMENT" = "BATCH" ]; then
    snakemake \
    --cluster "sbatch --export=ALL --partition={resources.queue} --qos=hotel --output=/dev/null --time='00:{resources.runtime_min}:00' --nodes=1 --ntasks-per-node={threads}" \
    --cluster-cancel "scancel {cluster.jobid}" \
    --cluster-sync "srun --time {cluster.time} --mem {cluster.mem} --cpus-per-task {cluster.cpus}" \
    --default-resources 'runtime_min=30' 'queue="condo"' \
    --latency-wait 60 \
    --use-conda \
    --conda-frontend conda \
    --rerun-trigger {mtime,params,input} \
    -k \
    -j 12 \
    -c 12 \
    "$@" &>"$out_path/logs/log"
else
    snakemake \
    --latency-wait 60 \
    --use-conda \
    --conda-frontend conda \
    --notemp \
    --rerun-trigger {mtime,params,input} \
    -k \
    -c 4 \
    "$@" &>"$out_path/logs/log"
fi

exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
        slack "snakemake finished successfully" &>/dev/null
    else
        slack "snakemake simulate_gwas job failed" &>/dev/null
        slack "$(tail -n4 "$out_path/logs/log")" &>/dev/null
    fi
fi
exit "$exit_code"
