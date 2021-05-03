#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=norm_fasta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

module purge

module load gatk/4.1.9.0

gatk --java-options "-Xmx2G" NormalizeFasta \
INPUT=Homo_sapiens.GRCh38.cdna.all.fa \
OUTPUT=Homo_sapiens.GRCh38.cdna.all.normalized.fa \
TRUNCATE_SEQUENCE_NAMES_AT_WHITESPACE=true

echo _ESTATUS_ [ norm fasta ]: $?
echo _END_ [ normalizefasta.slurm ]: $(date)
