#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=samtools_index
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu
#SBATCH --array=1-8

module purge

echo script begin: $(date)

# Path to 3-column (tab-delimited) table with sample name, fastq 1 file name, and fastq 2 file name

table=/scratch/zw2946/ngs.week8/task1/Aligned.sortedByCoord.out.bam.txt

line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
sample="$(printf "%s" "${line}" | basename $line Aligned.sortedByCoord.out.bam)"

cd $sample

module load samtools/intel/1.11

samtools index ${sample}Aligned.sortedByCoord.out.bam
