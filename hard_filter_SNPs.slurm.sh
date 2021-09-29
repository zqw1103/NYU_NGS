#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=hard_filter_SNPs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

module purge
echo script begin: $(date)

module load gatk/4.1.9.0

gatk VariantFiltration \
    -V /scratch/zw2946/ngs.week4/task2/snps.cohort.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O snps.filtered.vcf.gz

echo _ESTATUS_ [ subset snps ]: $?
echo _END_ [ subset_SNPs_temp.slurm ]: $(date)
