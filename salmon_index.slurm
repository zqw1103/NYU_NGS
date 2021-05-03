#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=salmon_index
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

module purge

module load salmon/1.4.0

salmon index -t Homo_sapiens.GRCh38.cdna.all.normalized.fa \
-i Homo_sapiens.GRCh38.cdna.all.normalized.fa_index -k 31

echo _ESTATUS_ [ salmon index ]: $?
echo _END_ [ salmon_index.slurm ]: $(date)

