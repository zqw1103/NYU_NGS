#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=genotypeGVCFS
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

module purge
echo script begin: $(date)

module load gatk/4.1.9.0

gatk --java-options "-Xmx8G" GenotypeGVCFs \
-R /scratch/work/courses/BI7653/hw3.2021/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa \
-V /scratch/work/courses/BI7653/hw4.2021/cohort.g.vcf.gz \
-O cohort.vcf.gz \
--allow-old-rms-mapping-quality-annotation-data

echo _ESTATUS_ [ genotype gvcfs ]: $?
echo _END_ [ GenotypeGVCF_tem.slurm ]: $(date)
