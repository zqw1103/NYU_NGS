#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --mem=46GB
#SBATCH --job-name=convert_sort_index
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zw2946@nyu.edu
#SBATCH --array=1-3

module purge


# Path to fq file

table=/scratch/zw2946/ngs.week11/task2/processed_fq_file.txt


# Define sample variables for current array index

fqfilepath="$(head -n ${SLURM_ARRAY_TASK_ID} "${table}" | tail -n 1)"
fq="$(basename $fqfilepath)"
sample="$(basename $fq .fq)"

echo Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample

cd "${sample}"

module load samtools/intel/1.11

#sam to bam
samtools view -bh ${sample}.sam > ${sample}.bam


module load picard/2.17.11

java -Xmx44g -jar "${PICARD_JAR}" SortSam \
INPUT=${sample}.bam \
OUTPUT=${sample}.sorted.bam \
SORT_ORDER=coordinate \
TMP_DIR="${SLURM_JOBTMP}" \
MAX_RECORDS_IN_RAM=10000000 \
VALIDATION_STRINGENCY=LENIENT

echo _ESTATUS_ [ SortSam $sample ]: $?

module load samtools/intel/1.11

#index sorted bam
samtools index ${sample}.sorted.bam

echo _ESTATUS_ [ index $sample ]: $?

echo _END_ [ convert_sort_index for $sample ]: $(date)