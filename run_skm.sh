#!/bin/bash
#SBATCH -A sogenant_supergene_mgraminicola
#SBATCH --mem=50GB
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --partition=fast
#SBATCH --cpus-per-task=32
#SBATCH --job-name=snakemake_pipeline
#SBATCH -o snakemake_controller.log
#SBATCH -e snakemake_controller.err

module load snakemake/9.4.0
module load python/3.12

mkdir -p conda_cache/pkgs
export CONDA_PKGS_DIRS=conda_cache/pkgs


CONDARC=.condarc snakemake --use-conda --conda-frontend mamba --executor slurm --profile profiles/slurm

