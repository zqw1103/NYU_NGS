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
#SBATCH --array=1-6

module purge

table=/scratch/work/courses/BI7653/project.2021/project_fastqs.txt
line="$(head -n ${SLURM_ARRAY_TASK_ID} "${table}" | tail -n 1)"
fq="$(printf "%s" "${line}" | cut -f2)"
sample="$(basename "${fq}" .fastq.gz)"

echo The array index is: ${SLURM_ARRAY_TASK_ID} sample: $sample

fqdir=/scratch/work/courses/BI7653/project.2021/fastqs

mkdir "${sample}"
cd "${sample}"

fq_fastp=$(basename ${fq} .fastq.gz).sE.fastq.gz

module load fastp/intel/0.20.1

fastp -i $fqdir/$fq \
-o $fq_fastp \
-g \
--length_required 50 \
--html $sample.fastp.html \
--json $sample.fastp.json

echo _ESTATUS_ [ fastp for $sample ]: $?

module purge
module load fastqc/0.11.9

fastqc $fq_fastp

echo _ESTATUS_ [ fastqc for $sample ]: $?

echo _END_ [ fastp for $sample ]: $(date)
