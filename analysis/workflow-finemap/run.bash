#!/usr/bin/env bash
#SBATCH --export ALL
#SBATCH --partition ind-shared
#SBATCH --account ddp268
#SBATCH --job-name happler-smk
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --time 24:00:00
#SBATCH --mem 20G
#SBATCH --output /dev/null
#SBATCH --signal=B:SIGUSR1@5
#SBATCH --mail-type=END,FAIL


# An example bash script demonstrating how to run the entire snakemake pipeline
# This script creates a log file in the execution directory
echo "" > smk.log

# try to find and activate the snakemake conda env if we need it
if ! command -v 'snakemake' &>/dev/null && \
	command -v 'conda' &>/dev/null && \
   [ "$CONDA_DEFAULT_ENV" != "snakemake" ] && \
   conda info --envs | grep "$CONDA_ROOT/snakemake" &>/dev/null; then
        echo "Snakemake not detected. Attempting to switch to snakemake environment." >> "smk.log"
        eval "$(conda shell.bash hook)"
        conda activate snakemake
fi

snakemake \
--workflow-profile workflow-finemap/profile/default \
-s workflow-finemap/Snakefile \
-c 4 \
"$@" &>"smk.log" &

wait $!
exit_code=$?

exit "$exit_code"
