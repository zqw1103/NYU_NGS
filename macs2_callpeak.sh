#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=46GB
#SBATCH --job-name=macs2_callpeak
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu


module load macs2/intel/2.2.7.1

macs2 callpeak -t /scratch/zw2946/ngs.week11/task3/SRR7207011/SRR7207011.filtered.sorted.bam \
-c /scratch/zw2946/ngs.week11/task3/SRR7207089/SRR7207089.filtered.sorted.bam \
-f BAM \
-g hs \
-n SRR7207011 \
-B \
-q 0.01

macs2 callpeak -t /scratch/zw2946/ngs.week11/task3/SRR7207017/SRR7207017.filtered.sorted.bam \
-c /scratch/zw2946/ngs.week11/task3/SRR7207089/SRR7207089.filtered.sorted.bam \
-f BAM \
-g hs \
-n SRR7207017 \
-B \
-q 0.01

echo _ESTATUS_ [ macs2_callpeak ]: $?

echo _END_ [ macs2_callpeak.slurm ]: $(date)