#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=subset_SNPs_VCF
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

module purge
echo script begin: $(date)

module load gatk/4.1.9.0

gatk SelectVariants \
-V /scratch/zw2946/ngs.week4/task1/cohort.vcf.gz \
-select-type SNP \
-O snps.cohort.vcf.gz

echo _ESTATUS_ [ subset snps ]: $?
echo _END_ [ subset_SNPs_temp.slurm ]: $(date)
