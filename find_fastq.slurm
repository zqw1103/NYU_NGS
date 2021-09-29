#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=find_fastq
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

find $PWD -name \
'SRR*.fastq' > fastq_files.txt
