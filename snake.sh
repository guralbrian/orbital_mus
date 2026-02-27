#!/bin/bash
# Run Snakemake on the login node (lightweight scheduler) and let it submit
# compute-heavy rules as SLURM jobs. Use: nohup bash snake.sh &
#
# Snakemake itself uses negligible resources — only the individual rules
# run on compute nodes via the SLURM executor plugin.

eval "$(conda shell.bash hook)"
conda activate orbital_mus

snakemake -j 64 --latency-wait 120 \
    --executor slurm \
    --rerun-incomplete \
    --default-resources mem_mb=4000 runtime=120 nodes=1 tasks=1 cpus_per_task=4 slurm_account=rc_chrau_pi

snakemake --rulegraph | dot -Tsvg > dag.svg
