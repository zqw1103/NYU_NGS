#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=samtools.stats.slurm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

module purge

echo script begin: $(date)

module load samtools/intel/1.11

samtools stats /scratch/work/courses/BI7653/hw6.2021/CR2342.bam > CR2342_bam_samtools_stats.txt

echo _ESTATUS_ [ samtools.stats ]: $?

echo _END_ [ samtools_stats.slurm ]: $(date)
