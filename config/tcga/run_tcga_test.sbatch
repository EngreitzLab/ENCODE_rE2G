#!/bin/bash
#SBATCH --job-name=rE2G_tcga_test
#SBATCH --account=engreitz
#SBATCH --partition=engreitz
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/scratch/users/shethm/tcga_atac/rE2G_results/logs/%x_%j.log
set -euo pipefail
export PYTHONUNBUFFERED=1

cd "$OAK/Users/sheth/ENCODE_rE2G_tcga"
eval "$(conda shell.bash hook)"
conda activate /oak/stanford/groups/engreitz/Users/sheth/.conda/envs/run_snakemake
echo "host $(hostname), snakemake $(snakemake --version), mamba $(command -v mamba)"

snakemake --profile slurm --use-conda --configfile config/tcga/config_tcga_test.yaml
