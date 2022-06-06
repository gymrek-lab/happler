#!/usr/bin/env bash
#PBS -V
#PBS -d .
#PBS -q hotel
#PBS -j oe
#PBS -o /dev/null
#PBS -N run.happler
#PBS -l nodes=1:ppn=1
#PBS -l walltime=1:00:00


# An example bash script demonstrating how to run the entire happler pipeline
# This script creates a log file in the output dir:
#   log - the basic happler log of completed rules

# Before running this happler pipeline, remember to complete the config file
# with the required input info.
# Also, make sure that this script is executed from the directory that it lives in!

out_path="out"

# try to find and activate the happler conda env if we need it
if ! command -v 'happler' &>/dev/null && \
    command -v 'conda' &>/dev/null && \
   [ "$CONDA_DEFAULT_ENV" != "happler" ] && \
   conda info --envs | grep "$CONDA_ROOT/happler" &>/dev/null; then
        echo "Snakemake not detected. Attempting to switch to happler environment." >> "out/logs/log"
        eval "$(conda shell.bash hook)"
        conda activate happler
fi

happler run \
--region '1:98001984-99001984' \
-S out/1_98001984-99001984/samples-10000.tsv \
--discard-multiallelic \
-o out/1_98001984-99001984/happler.haps \
--verbosity DEBUG \
out/1_98001984-99001984/merged.sorted.bcf \
<(zcat out/1_98001984-99001984/sim/str/0.15/phens.tsv.gz | tail -n+2 | cut -f1,3) &>"$out_path/logs/log"

exit_code="$?"
if command -v 'slack' &>/dev/null; then
    if [ "$exit_code" -eq 0 ]; then
        slack "snakemake finished successfully" &>/dev/null
    else
        slack "snakemake run_geuvadis job failed" &>/dev/null
        slack "$(tail -n4 "$out_path/logs/log")" &>/dev/null
    fi
fi
exit "$exit_code"
