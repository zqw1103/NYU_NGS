#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=bwa_mem
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu
#SBATCH --array=1-3

module purge

ref=/scratch/work/courses/BI7653/hw3.2021/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa

table=/scratch/zw2946/ngs.week11/task2/processed_fq_file.txt

fqfilepath="$(head -n ${SLURM_ARRAY_TASK_ID} "${table}" | tail -n 1)"
fq="$(basename $fqfilepath)"
sample="$(basename $fq .fq)"

mkdir "${sample}"
cd "${sample}"

echo Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample

module load bwa/intel/0.7.17

bwa mem \
-M \
-t $SLURM_CPUS_PER_TASK \
"${ref}" \
$fqfilepath > $sample.sam

echo _ESTATUS_ [ bwa mem for $sample ]: $?

echo _END_ [ bwa_mem.slurm for $sample ]: $(date)
