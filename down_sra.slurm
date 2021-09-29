#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=down_sra
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu
##SBATCH --array=1-3

module purge 

module load sra-tools/2.10.9 

fasterq-dump SRR7207011

fasterq-dump SRR7207017

fasterq-dump SRR7207089

echo _ESTATUS_ [ down_sra ]: $?


