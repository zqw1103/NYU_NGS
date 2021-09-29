#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=filtered_SNPs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

module purge
echo script begin: $(date)

module load gatk/4.1.9.0

gatk SelectVariants \
-R /scratch/work/courses/BI7653/hw3.2021/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa \
-V snps.filtered.vcf.gz \
--exclude-filtered \
-O filtered.snps.vcf.gz

echo _ESTATUS_ [ subset snps ]: $?
echo _END_ [ subset_SNPs_temp.slurm ]: $(date)
