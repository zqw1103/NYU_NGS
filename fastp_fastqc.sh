#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=fastp_array
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu
#SBATCH --array=1-3

module purge

table=/scratch/zw2946/ngs.week11/raw_data/fastq_files.txt

fqfilepath="$(head -n ${SLURM_ARRAY_TASK_ID} "${table}" | tail -n 1)"
fq="$(basename $fqfilepath)"
sample="$(basename $fq .fastq)"

mkdir "${sample}"
cd "${sample}"

echo Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample

fq_fastp=$(basename ${fq} .fastq).fq

module load fastp/intel/0.20.1

fastp -i $fqfilepath \
-o $fq_fastp \
--length_required 50 \
--html $sample.fastp.html \
--json $sample.fastp.json

echo _ESTATUS_ [ fastp for $sample ]: $?

module purge
module load fastqc/0.11.9

fastqc $fq_fastp

echo _ESTATUS_ [ fastqc for $sample ]: $?

echo _END_ [ fastp for $sample ]: $(date)