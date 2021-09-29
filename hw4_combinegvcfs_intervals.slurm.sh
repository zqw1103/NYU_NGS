#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --mem=44GB
#SBATCH --job-name=intervals
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jmf11@nyu.edu


module purge

ref=/scratch/courses/BI7653/hw3.2019/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa

module load gatk/4.0.2.1

gatk --java-options "-Xmx42G" CombineGVCFs \
-R $ref \
-L 1:1-5000000 \
-L 2:1-5000000 \
-L 3:1-5000000 \
--variant /scratch/work/courses/BI7653/hw4.2021/gvcfs.list \
-O cohort.intervals.g.vcf.gz

echo _ESTATUS_ [ combine gvcfs intervals ]: $?
echo _END_ [ hw4_combinegvcfs_intervals.slurm ]: $(date)
