#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=find_processed_fq
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu

find $PWD -name "*.fq" > processed_fq_file.txt
