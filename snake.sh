#!/bin/bash
#SBATCH --job-name=orbital_mus
#SBATCH --time=240
#SBATCH --mem=1000
#SBATCH --output=.slurmlogs/%j.out
#SBATCH --error=.slurmlogs/%j.err
#SBATCH -N 1
#SBATCH -n 1

eval "$(conda shell.bash hook)"
conda activate orbital_mus

snakemake -j 4 --latency-wait 60 \
    --cluster "sbatch --mem={resources.mem_mb} -N 1 -n 4 -t 120 -o .slurmlogs/%j.out"

snakemake --rulegraph | dot -Tsvg > dag.svg
